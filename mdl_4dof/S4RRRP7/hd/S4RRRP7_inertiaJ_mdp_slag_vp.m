% Calculate joint inertia matrix for
% S4RRRP7
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
%   see S4RRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP7_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:05
% EndTime: 2019-12-31 17:21:06
% DurationCPUTime: 0.28s
% Computational Cost: add. (181->89), mult. (356->128), div. (0->0), fcn. (261->4), ass. (0->31)
t55 = sin(qJ(2));
t74 = 0.2e1 * t55;
t54 = sin(qJ(3));
t73 = pkin(5) * t54;
t56 = cos(qJ(3));
t72 = pkin(5) * t56;
t57 = cos(qJ(2));
t71 = t54 * t57;
t48 = -t57 * pkin(2) - t55 * pkin(6) - pkin(1);
t43 = t54 * t48 + t57 * t72;
t51 = t54 ^ 2;
t53 = t56 ^ 2;
t70 = t51 + t53;
t69 = t54 * MDP(13);
t68 = t56 * MDP(12);
t67 = MDP(17) - MDP(20);
t66 = -MDP(21) * pkin(3) - MDP(18);
t65 = -2 * qJ(4) * MDP(20) - MDP(15);
t64 = -t56 * pkin(3) - t54 * qJ(4);
t63 = -pkin(3) * t54 + t56 * qJ(4);
t40 = -t57 * qJ(4) + t43;
t46 = t56 * t48;
t41 = -t46 + (pkin(3) + t73) * t57;
t62 = t40 * t56 + t41 * t54;
t61 = t56 * MDP(13) - t54 * MDP(14);
t60 = (-pkin(5) * t71 + t46) * MDP(16) - t43 * MDP(17);
t59 = t54 * MDP(18) - t56 * MDP(20);
t49 = pkin(6) * t71;
t47 = -pkin(2) + t64;
t44 = (pkin(5) - t63) * t55;
t1 = [MDP(1) + (t40 ^ 2 + t41 ^ 2 + t44 ^ 2) * MDP(21) + (t57 * MDP(15) + 0.2e1 * pkin(1) * MDP(9) + (MDP(5) - t61) * t74) * t57 + 0.2e1 * (t41 * MDP(18) - t40 * MDP(20) - t60) * t57 + ((-t40 * t54 + t41 * t56) * MDP(19) + t59 * t44) * t74 + (-0.2e1 * pkin(1) * MDP(10) + (t53 * MDP(11) - 0.2e1 * t54 * t68 + MDP(4) + 0.2e1 * (t54 * MDP(16) + t56 * MDP(17)) * pkin(5)) * t55) * t55; t49 * MDP(16) + (-t44 * t56 + t49) * MDP(18) + t62 * MDP(19) - t44 * t54 * MDP(20) + (t62 * pkin(6) + t44 * t47) * MDP(21) + (-pkin(5) * MDP(10) - t69 + MDP(7) + (t67 * pkin(6) - MDP(14)) * t56) * t57 + (MDP(6) - pkin(5) * MDP(9) + t56 * t54 * MDP(11) + (-t51 + t53) * MDP(12) + (-pkin(2) * t54 - t72) * MDP(16) + (-pkin(2) * t56 + t73) * MDP(17) + t59 * t47) * t55; MDP(8) + t51 * MDP(11) + (t70 * pkin(6) ^ 2 + t47 ^ 2) * MDP(21) + 0.2e1 * t70 * MDP(19) * pkin(6) + 0.2e1 * (pkin(2) * MDP(16) - t47 * MDP(18)) * t56 + 0.2e1 * (-pkin(2) * MDP(17) - t47 * MDP(20) + t68) * t54; t46 * MDP(18) + t43 * MDP(20) + (-t41 * pkin(3) + t40 * qJ(4)) * MDP(21) + ((-0.2e1 * pkin(3) - t73) * MDP(18) + t65) * t57 + (t64 * MDP(19) + t61) * t55 + t60; t69 + t56 * MDP(14) + t63 * MDP(19) + ((MDP(21) * qJ(4) - t67) * t56 + (-MDP(16) + t66) * t54) * pkin(6); 0.2e1 * pkin(3) * MDP(18) + (pkin(3) ^ 2 + (qJ(4) ^ 2)) * MDP(21) - t65; t56 * t55 * MDP(19) + t57 * MDP(18) + t41 * MDP(21); (MDP(21) * pkin(6) + MDP(19)) * t54; t66; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
