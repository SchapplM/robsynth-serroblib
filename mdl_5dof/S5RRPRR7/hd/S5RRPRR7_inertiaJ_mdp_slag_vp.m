% Calculate joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR7_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:41
% EndTime: 2019-12-31 20:15:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (184->68), mult. (286->81), div. (0->0), fcn. (234->6), ass. (0->40)
t69 = sin(qJ(5));
t70 = sin(qJ(4));
t72 = cos(qJ(5));
t73 = cos(qJ(4));
t55 = t69 * t73 + t72 * t70;
t56 = -t69 * t70 + t72 * t73;
t98 = t55 * MDP(22) + t56 * MDP(23);
t97 = -2 * MDP(7);
t96 = 2 * MDP(8);
t74 = cos(qJ(2));
t65 = -t74 * pkin(1) - pkin(2);
t61 = -pkin(7) + t65;
t95 = -pkin(8) + t61;
t75 = -pkin(2) - pkin(7);
t94 = -pkin(8) + t75;
t93 = MDP(9) * pkin(2);
t92 = t56 * MDP(22) - t55 * MDP(23);
t91 = t56 * MDP(19) - t55 * MDP(20);
t89 = t65 * MDP(9);
t85 = t70 * MDP(15);
t84 = t73 * MDP(15);
t83 = t73 * MDP(16);
t71 = sin(qJ(2));
t62 = t71 * pkin(1) + qJ(3);
t48 = t95 * t70;
t49 = t95 * t73;
t82 = (-t69 * t48 + t72 * t49) * MDP(22) + (-t72 * t48 - t69 * t49) * MDP(23) + t91;
t57 = t94 * t70;
t58 = t94 * t73;
t81 = (-t69 * t57 + t72 * t58) * MDP(22) + (-t72 * t57 - t69 * t58) * MDP(23) + t91;
t80 = MDP(4) + (MDP(10) * t73 - 0.2e1 * MDP(11) * t70) * t73 + (MDP(17) * t56 - 0.2e1 * MDP(18) * t55) * t56;
t79 = -t70 * MDP(16) + t84;
t78 = t96 + 0.2e1 * t83 + 0.2e1 * t85;
t77 = (MDP(22) * t72 - MDP(23) * t69) * pkin(4);
t76 = 0.2e1 * t98;
t68 = t70 * pkin(4);
t67 = t73 * MDP(12);
t64 = qJ(3) + t68;
t59 = t62 + t68;
t1 = [MDP(1) + ((2 * MDP(7)) + t89) * t65 + t59 * t76 + 0.2e1 * (t74 * MDP(5) - t71 * MDP(6)) * pkin(1) + (t62 * MDP(9) + t78) * t62 + t80; pkin(2) * t97 + qJ(3) * t96 + (-t65 * pkin(2) + t62 * qJ(3)) * MDP(9) + ((MDP(5) - MDP(7)) * t74 + (-MDP(6) + MDP(8)) * t71) * pkin(1) + t80 + t98 * (t59 + t64) + (t83 + t85) * (qJ(3) + t62); t64 * t76 + (t97 + t93) * pkin(2) + (MDP(9) * qJ(3) + t78) * qJ(3) + t80; MDP(7) + t89; MDP(7) - t93; MDP(9); -t70 * MDP(13) + t61 * t79 + t67 + t82; t75 * t84 + t67 + (-MDP(16) * t75 - MDP(13)) * t70 + t81; t79 + t92; MDP(14) + MDP(21) + 0.2e1 * t77; t82; t81; t92; MDP(21) + t77; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
