% Calculate joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR9_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:05
% EndTime: 2019-12-31 17:10:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (195->80), mult. (416->126), div. (0->0), fcn. (378->6), ass. (0->43)
t95 = (MDP(14) * qJ(3));
t94 = 2 * MDP(13);
t93 = -2 * MDP(16);
t92 = 2 * MDP(21);
t68 = sin(pkin(7));
t91 = pkin(5) * t68;
t73 = cos(qJ(2));
t90 = pkin(5) * t73;
t71 = sin(qJ(2));
t89 = pkin(6) * t71;
t88 = pkin(6) + qJ(3);
t59 = -t73 * pkin(2) - t71 * qJ(3) - pkin(1);
t69 = cos(pkin(7));
t51 = t68 * t59 + t69 * t90;
t86 = pkin(2) * MDP(14);
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t57 = t72 * t68 + t70 * t69;
t52 = t57 * t71;
t85 = t52 * MDP(18);
t56 = t70 * t68 - t72 * t69;
t53 = t56 * t71;
t84 = t53 * MDP(15);
t83 = t53 * MDP(17);
t82 = t56 * MDP(20);
t81 = t68 * MDP(12);
t80 = t69 * MDP(11);
t79 = t73 * MDP(19);
t77 = MDP(11) * t68 + MDP(12) * t69;
t76 = -t80 + t81 - t86;
t60 = t88 * t68;
t61 = t88 * t69;
t75 = t57 * MDP(17) - t56 * MDP(18) + (-t72 * t60 - t70 * t61) * MDP(20) - (-t70 * t60 + t72 * t61) * MDP(21);
t67 = t71 ^ 2;
t64 = -t69 * pkin(3) - pkin(2);
t58 = (pkin(3) * t68 + pkin(5)) * t71;
t55 = t69 * t59;
t50 = -t68 * t90 + t55;
t47 = -t68 * t89 + t51;
t46 = -t69 * t89 + t55 + (-pkin(3) - t91) * t73;
t45 = t70 * t46 + t72 * t47;
t44 = t72 * t46 - t70 * t47;
t1 = [MDP(1) + t67 * MDP(4) + (t67 * pkin(5) ^ 2 + t50 ^ 2 + t51 ^ 2) * MDP(14) - (t52 * t93 - t84) * t53 + (0.2e1 * pkin(1) * MDP(9) + t79 + 0.2e1 * t83 + 0.2e1 * t85) * t73 + 0.2e1 * (-t50 * t73 + t67 * t91) * MDP(11) + 0.2e1 * (t67 * pkin(5) * t69 + t51 * t73) * MDP(12) + 0.2e1 * (-t44 * t73 + t58 * t52) * MDP(20) + (t45 * t73 - t58 * t53) * t92 + (-0.2e1 * pkin(1) * MDP(10) + 0.2e1 * t73 * MDP(5) + (-t50 * t69 - t51 * t68) * t94) * t71; -t57 * t84 + (-t57 * t52 + t53 * t56) * MDP(16) + (t64 * t52 + t58 * t56) * MDP(20) + (-t64 * t53 + t58 * t57) * MDP(21) + (-pkin(5) * MDP(10) + t77 * qJ(3) + MDP(7) - t75) * t73 + (MDP(6) - t77 * pkin(2) + (-MDP(9) + t76) * pkin(5)) * t71 + (MDP(13) + t95) * (-t50 * t68 + t51 * t69); 0.2e1 * t64 * t82 + MDP(8) + (0.2e1 * t80 - 0.2e1 * t81 + t86) * pkin(2) + (MDP(15) * t57 + t56 * t93 + t64 * t92) * t57 + (t94 + t95) * (t68 ^ 2 + t69 ^ 2) * qJ(3); t52 * MDP(20) - t53 * MDP(21) + (pkin(5) * MDP(14) + t77) * t71; t57 * MDP(21) + t76 + t82; MDP(14); t44 * MDP(20) - t45 * MDP(21) - t79 - t83 - t85; t75; 0; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
