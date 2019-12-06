% Calculate joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPPR1_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:07
% EndTime: 2019-12-05 17:29:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (201->64), mult. (382->96), div. (0->0), fcn. (343->8), ass. (0->36)
t71 = cos(pkin(7));
t61 = -pkin(1) * t71 - pkin(2);
t67 = sin(pkin(8));
t70 = cos(pkin(8));
t54 = -pkin(3) * t70 - qJ(4) * t67 + t61;
t69 = cos(pkin(9));
t52 = t69 * t54;
t68 = sin(pkin(7));
t60 = pkin(1) * t68 + qJ(3);
t66 = sin(pkin(9));
t84 = pkin(6) * t67;
t45 = -t69 * t84 + t52 + (-t60 * t66 - pkin(4)) * t70;
t83 = t60 * t70;
t48 = t66 * t54 + t69 * t83;
t46 = -t66 * t84 + t48;
t72 = sin(qJ(5));
t73 = cos(qJ(5));
t56 = t66 * t73 + t69 * t72;
t49 = t56 * t67;
t88 = -(t45 * t73 - t46 * t72) * MDP(18) + (t45 * t72 + t46 * t73) * MDP(19) + t49 * MDP(16);
t87 = 0.2e1 * t67;
t86 = -0.2e1 * t70;
t82 = t66 ^ 2 + t69 ^ 2;
t81 = t61 * MDP(8);
t47 = -t66 * t83 + t52;
t79 = t47 * t69 + t48 * t66;
t55 = -t66 * t72 + t69 * t73;
t78 = t69 * MDP(10) + t66 * MDP(9);
t50 = t55 * t67;
t76 = t49 * MDP(18) + t50 * MDP(19);
t75 = MDP(18) * t55 - MDP(19) * t56;
t65 = t70 ^ 2;
t63 = t67 ^ 2;
t59 = t60 ^ 2;
t58 = t63 * t59;
t1 = [MDP(1) + (t59 * t65 + t58) * MDP(8) + (t47 ^ 2 + t48 ^ 2 + t58) * MDP(12) + t65 * MDP(17) + (t68 ^ 2 + t71 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(5) * t86 + MDP(6) * t87 + t81) * t61 + (MDP(13) * t50 - 0.2e1 * t49 * MDP(14) + MDP(15) * t86) * t50 + t76 * (pkin(4) * t66 + t60) * t87 - 0.2e1 * t79 * MDP(11) * t67 + 0.2e1 * (t65 * MDP(7) + (MDP(7) + t78) * t63) * t60 + 0.2e1 * (t48 * MDP(10) - t47 * MDP(9) + t88) * t70; (-t47 * t66 + t48 * t69 - t83) * t67 * MDP(12); MDP(4) + (t63 + t65) * MDP(8) + (t82 * t63 + t65) * MDP(12); t81 + t79 * MDP(12) + (-t82 * MDP(11) + MDP(6)) * t67 + (MDP(10) * t66 - MDP(9) * t69 - MDP(5) - t75) * t70; 0; t82 * MDP(12) + MDP(8); (MDP(12) * t60 + t78) * t67 + t76; -t70 * MDP(12); 0; MDP(12); t50 * MDP(15) - MDP(17) * t70 - t88; -t76; t75; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
