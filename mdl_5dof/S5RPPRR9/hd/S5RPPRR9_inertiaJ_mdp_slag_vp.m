% Calculate joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR9_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:49
% EndTime: 2019-12-31 18:02:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (143->68), mult. (237->101), div. (0->0), fcn. (186->6), ass. (0->32)
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t59 = -pkin(1) - pkin(2);
t46 = t54 * qJ(2) + t53 * t59;
t43 = -pkin(6) + t46;
t55 = sin(qJ(5));
t73 = t43 * t55;
t57 = cos(qJ(5));
t72 = t43 * t57;
t71 = t55 * t57;
t58 = cos(qJ(4));
t70 = t55 * t58;
t69 = t57 * t58;
t56 = sin(qJ(4));
t68 = t56 * MDP(16);
t67 = t58 * MDP(21);
t66 = MDP(18) * t71;
t44 = t53 * qJ(2) - t54 * t59;
t42 = pkin(3) + t44;
t65 = -MDP(19) * t57 + MDP(20) * t55;
t64 = (-t53 * t70 - t57 * t54) * MDP(22) - (t53 * t69 - t55 * t54) * MDP(23);
t63 = t57 * MDP(22) - t55 * MDP(23);
t62 = -MDP(22) * t55 - MDP(23) * t57;
t61 = MDP(15) + t63;
t60 = t55 * MDP(19) + t57 * MDP(20) + t62 * pkin(7);
t52 = t57 ^ 2;
t51 = t56 ^ 2;
t50 = t55 ^ 2;
t39 = t58 * pkin(4) + t56 * pkin(7) + t42;
t38 = t55 * t39 + t43 * t69;
t37 = t57 * t39 - t43 * t70;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t44 ^ 2 + t46 ^ 2) * MDP(9) - 0.2e1 * t42 * t68 + (t52 * MDP(17) + MDP(10) - 0.2e1 * t66) * t51 + (0.2e1 * t42 * MDP(15) + t67 + 0.2e1 * (MDP(11) + t65) * t56) * t58 + 0.2e1 * t44 * MDP(7) + 0.2e1 * t46 * MDP(8) + 0.2e1 * (t37 * t58 - t51 * t73) * MDP(22) + 0.2e1 * (-t38 * t58 - t51 * t72) * MDP(23); -pkin(1) * MDP(6) - MDP(4) + t64 * t58 + (t46 * MDP(9) + t62 * t51 + MDP(8)) * t53 + (-t58 * MDP(15) - t44 * MDP(9) - MDP(7) + t68) * t54; MDP(6) + (t53 ^ 2 + t54 ^ 2) * MDP(9); 0; 0; MDP(9); (-t43 * MDP(16) - MDP(13) + t60) * t58 + (-MDP(12) - t43 * MDP(15) - MDP(17) * t71 + (t50 - t52) * MDP(18) + (pkin(4) * t55 - t72) * MDP(22) + (pkin(4) * t57 + t73) * MDP(23)) * t56; (-t58 * MDP(16) - t61 * t56) * t53; t61 * t58 - t68; t50 * MDP(17) + 0.2e1 * pkin(4) * t63 + MDP(14) + 0.2e1 * t66; t37 * MDP(22) - t38 * MDP(23) + t65 * t56 + t67; t64; t62 * t56; t60; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
