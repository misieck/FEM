function snapshot_index = find_stationary(V, max_v, test_length)

stable_snapshots = (find (max (V(:,3:end)) < max_v));
mod_stable = [stable_snapshots, stable_snapshots(end)+1:stable_snapshots(end)+test_length+1];
for i = 1:length(stable_snapshots)
    mod_stable(i) = mod_stable(i+test_length) - mod_stable(i);
end
stable_indises = find (mod_stable == test_length);
snapshot_index = stable_snapshots(stable_indises(1));
