% Calculate joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR9_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:38
% EndTime: 2019-12-31 18:24:39
% DurationCPUTime: 0.13s
% Computational Cost: add. (135->67), mult. (242->92), div. (0->0), fcn. (171->6), ass. (0->39)
t58 = -pkin(3) - pkin(7);
t52 = sin(pkin(8));
t46 = t52 * pkin(1) + pkin(6);
t78 = pkin(4) + t46;
t57 = cos(qJ(3));
t44 = t78 * t57;
t77 = t44 * t57;
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t76 = t54 * t56;
t55 = sin(qJ(3));
t49 = t55 ^ 2;
t51 = t57 ^ 2;
t75 = t49 + t51;
t74 = MDP(15) * pkin(3);
t73 = t46 ^ 2 * MDP(15);
t72 = t54 * MDP(21);
t71 = t56 * MDP(22);
t70 = qJ(4) * MDP(15);
t69 = MDP(17) * t76;
t53 = cos(pkin(8));
t47 = -t53 * pkin(1) - pkin(2);
t68 = MDP(13) - t74;
t67 = MDP(10) - t68;
t66 = -MDP(11) + MDP(14) + t70;
t65 = -MDP(18) * t54 - MDP(19) * t56;
t64 = t56 * MDP(21) - t54 * MDP(22);
t63 = t71 + t72;
t62 = -t55 * qJ(4) + t47;
t61 = MDP(12) + t64;
t60 = (MDP(21) * t58 + MDP(18)) * t56 + (-MDP(22) * t58 - MDP(19)) * t54;
t50 = t56 ^ 2;
t48 = t54 ^ 2;
t43 = t78 * t55;
t42 = -t57 * pkin(3) + t62;
t41 = t58 * t57 + t62;
t40 = t56 * t41 + t54 * t43;
t39 = -t54 * t41 + t56 * t43;
t1 = [-0.2e1 * t47 * t57 * MDP(10) + MDP(1) + (t52 ^ 2 + t53 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * t57 * MDP(13) + MDP(15) * t42) * t42 + (t48 * MDP(16) + 0.2e1 * t69 + t73) * t51 + (MDP(20) + MDP(5) + t73) * t49 + 0.2e1 * (t47 * MDP(11) - t42 * MDP(14) + (MDP(6) + t65) * t57) * t55 + 0.2e1 * (t39 * t55 + t56 * t77) * MDP(21) + 0.2e1 * (-t40 * t55 - t54 * t77) * MDP(22) + 0.2e1 * t75 * MDP(12) * t46; 0; t75 * MDP(15) + MDP(4); t63 * t44 + (-pkin(3) * MDP(12) + MDP(7) + t60) * t55 + (MDP(8) - MDP(16) * t76 + (t48 - t50) * MDP(17) + t61 * qJ(4)) * t57 + (-t67 * t55 + t66 * t57) * t46; t67 * t57 + (t63 + t66) * t55; -0.2e1 * t69 + t50 * MDP(16) + MDP(9) + (-0.2e1 * MDP(13) + t74) * pkin(3) + (0.2e1 * MDP(14) + t70 + 0.2e1 * t71 + 0.2e1 * t72) * qJ(4); (MDP(15) * t46 + t61) * t55; -t57 * MDP(15); t68; MDP(15); t55 * MDP(20) + t39 * MDP(21) - t40 * MDP(22) + t65 * t57; -t64 * t57; t60; t64; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
