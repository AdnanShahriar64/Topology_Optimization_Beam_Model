function [IEN, Type, beta] = Generate_IEN_Type(Nodes)

  X = Nodes(:,2:4);
  nEl = size(X,1) - 1;

  % ---------- IEN ----------
  IEN = [Nodes(1:end-1,1)'; Nodes(2:end,1)'];

  % ---------- unit direction (slope) of each element ----------
  d = diff(X);
  d = d ./ sqrt(sum(d.^2,2));

  % ---------- turn angle at each interior node (deg) ----------
  beta = zeros(1,nEl+1);
  for k = 2:nEl
      beta(k) = atan2d(norm(cross(d(k-1,:),d(k,:))), dot(d(k-1,:),d(k,:)));
  end

  % ---------- type: slope changes after a straight -> curved ----------
  kind = ones(1,nEl);
  for e = 2:nEl-1
      if kind(e-1) == 1 && beta(e) > 0.5
          kind(e) = 2;
      end
  end

  Type = [kind; 2*ones(1,nEl); ones(1,nEl); ones(1,nEl)];

end
