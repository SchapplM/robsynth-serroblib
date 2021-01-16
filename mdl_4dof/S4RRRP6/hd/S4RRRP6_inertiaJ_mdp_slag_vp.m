% Calculate joint inertia matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:16
% EndTime: 2021-01-15 14:39:17
% DurationCPUTime: 0.21s
% Computational Cost: add. (179->87), mult. (357->121), div. (0->0), fcn. (273->4), ass. (0->34)
t47 = sin(qJ(2));
t68 = 0.2e1 * t47;
t46 = sin(qJ(3));
t67 = pkin(5) * t46;
t49 = cos(qJ(2));
t66 = pkin(5) * t49;
t48 = cos(qJ(3));
t65 = t46 * t48;
t64 = qJ(4) + pkin(6);
t63 = MDP(21) * pkin(3);
t62 = qJ(4) * t47;
t41 = t64 * t48;
t61 = t41 * MDP(19);
t42 = -t48 * pkin(3) - pkin(2);
t60 = t42 * MDP(21);
t59 = t46 * MDP(14);
t58 = t46 * MDP(19);
t57 = t48 * MDP(18);
t56 = t48 * t66;
t55 = MDP(12) * t65;
t54 = -MDP(20) * pkin(3) + MDP(13);
t53 = MDP(16) * t46 + MDP(17) * t48;
t52 = t46 * MDP(18) + t48 * MDP(19);
t39 = -t49 * pkin(2) - t47 * pkin(6) - pkin(1);
t34 = t56 + (t39 - t62) * t46;
t37 = t48 * t39;
t51 = (-t46 * t66 + t37) * MDP(16) - (t46 * t39 + t56) * MDP(17) - t34 * MDP(19);
t50 = -t57 + t58 + t60;
t45 = t48 ^ 2;
t43 = t46 ^ 2;
t40 = t64 * t46;
t38 = (pkin(3) * t46 + pkin(5)) * t47;
t33 = -t48 * t62 + t37 + (-pkin(3) - t67) * t49;
t1 = [MDP(1) + (t33 ^ 2 + t34 ^ 2 + t38 ^ 2) * MDP(21) + (t49 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + (-t48 * MDP(13) + MDP(5) + t59) * t68) * t49 + 0.2e1 * (-t33 * MDP(18) - t51) * t49 + ((-t33 * t48 - t34 * t46) * MDP(20) + t52 * t38) * t68 + (-0.2e1 * pkin(1) * MDP(10) + (t45 * MDP(11) + 0.2e1 * t53 * pkin(5) + MDP(4) - 0.2e1 * t55) * t47) * t47; (-t33 * t46 + t34 * t48) * MDP(20) + (-t33 * t40 + t34 * t41) * MDP(21) + t50 * t38 + (-pkin(5) * MDP(10) - t46 * MDP(13) - t48 * MDP(14) + t40 * MDP(18) + t53 * pkin(6) + MDP(7) + t61) * t49 + (MDP(6) - pkin(5) * MDP(9) + MDP(11) * t65 + (-t43 + t45) * MDP(12) + (-pkin(2) * t46 - pkin(5) * t48) * MDP(16) + (-pkin(2) * t48 + t67) * MDP(17) + (t40 * t48 - t41 * t46) * MDP(20) + t52 * t42) * t47; MDP(8) + t43 * MDP(11) + 0.2e1 * t55 + 0.2e1 * (t40 * t46 + t41 * t48) * MDP(20) + (t40 ^ 2 + t41 ^ 2) * MDP(21) + (-0.2e1 * t57 + 0.2e1 * t58 + t60) * t42 + 0.2e1 * (t48 * MDP(16) - t46 * MDP(17)) * pkin(2); t33 * t63 + t37 * MDP(18) + (-MDP(15) + (-0.2e1 * pkin(3) - t67) * MDP(18)) * t49 + (-t59 + (-MDP(18) * qJ(4) + t54) * t48) * t47 + t51; -t61 + (-MDP(17) * pkin(6) + MDP(14)) * t48 - (MDP(18) + t63) * t40 + (-MDP(16) * pkin(6) + t54) * t46; MDP(15) + (0.2e1 * MDP(18) + t63) * pkin(3); t38 * MDP(21) + t52 * t47; t50; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
