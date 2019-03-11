% Calculate joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPPR8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:00
% EndTime: 2019-03-09 03:00:01
% DurationCPUTime: 0.26s
% Computational Cost: add. (232->102), mult. (333->123), div. (0->0), fcn. (221->4), ass. (0->47)
t67 = sin(qJ(6));
t69 = cos(qJ(6));
t76 = t67 * MDP(27) + t69 * MDP(28);
t98 = t76 - MDP(15) + MDP(20);
t97 = 2 * qJ(2);
t96 = 2 * MDP(14);
t95 = 2 * MDP(19);
t71 = (-pkin(3) - pkin(4));
t68 = sin(qJ(3));
t72 = -pkin(1) - pkin(7);
t89 = -qJ(5) - t72;
t52 = t89 * t68;
t94 = t52 * t68;
t93 = t67 * t69;
t62 = t68 ^ 2;
t70 = cos(qJ(3));
t64 = t70 ^ 2;
t92 = t62 + t64;
t91 = t68 * qJ(4);
t90 = t70 * qJ(4) - qJ(2);
t88 = t72 ^ 2 * MDP(17);
t60 = -pkin(8) + t71;
t84 = MDP(16) - MDP(13);
t83 = MDP(17) + MDP(21);
t82 = 0.2e1 * t70;
t81 = MDP(23) * t93;
t80 = pkin(3) * MDP(17) + MDP(14);
t79 = t71 * MDP(21) + MDP(19);
t78 = MDP(24) * t69 - MDP(25) * t67;
t77 = t69 * MDP(27) - t67 * MDP(28);
t75 = MDP(18) + t77;
t74 = -t67 * MDP(24) - t69 * MDP(25) - t76 * t60;
t73 = qJ(4) ^ 2;
t66 = qJ(4) + pkin(5);
t63 = t69 ^ 2;
t61 = t67 ^ 2;
t57 = t70 * pkin(3) + t91;
t56 = t68 * pkin(3) - t90;
t55 = -t71 * t70 + t91;
t54 = t89 * t70;
t51 = t92 * t72;
t50 = t71 * t68 + t90;
t49 = t70 * pkin(5) + t60 * t68 + t90;
t48 = -t54 * t70 - t94;
t47 = t67 * t49 + t69 * t54;
t46 = t69 * t49 - t67 * t54;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + (MDP(5) * t97) + ((pkin(1) ^ 2 + qJ(2) ^ 2) * MDP(6)) + t56 ^ 2 * MDP(17) + (t50 ^ 2 + t52 ^ 2 + t54 ^ 2) * MDP(21) + (MDP(26) + MDP(7) + t88) * t64 + (t63 * MDP(22) - 0.2e1 * t81 + t88) * t62 + ((qJ(2) * MDP(13)) - t56 * MDP(16) + t50 * MDP(18)) * t82 + ((MDP(12) * t97) + t56 * t96 + t50 * t95 + (-MDP(8) + t78) * t82) * t68 - 0.2e1 * t51 * MDP(15) + 0.2e1 * t48 * MDP(20) + 0.2e1 * (t46 * t70 - t67 * t94) * MDP(27) + 0.2e1 * (-t47 * t70 - t69 * t94) * MDP(28); t51 * MDP(17) + t48 * MDP(21) - (pkin(1) * MDP(6)) + t98 * t92 + MDP(4); t83 * t92 + MDP(6); -t57 * MDP(15) + t55 * MDP(20) + t79 * t54 - (qJ(4) * MDP(21) + t75) * t52 + (MDP(9) + (MDP(12) + t80) * t72 + t74) * t70 + (-MDP(10) - MDP(22) * t93 + (t61 - t63) * MDP(23) + t76 * t66 + (qJ(4) * MDP(17) + t84) * t72) * t68; t57 * MDP(17) + t55 * MDP(21) + (MDP(12) + MDP(14) - MDP(19)) * t70 + (t75 + t84) * t68; MDP(11) + (pkin(3) * t96) + ((pkin(3) ^ 2) + t73) * MDP(17) + (t71 * t95) + ((t71 ^ 2) + t73) * MDP(21) + t61 * MDP(22) + 0.2e1 * t81 + 0.2e1 * t77 * t66 + 0.2e1 * (MDP(16) + MDP(18)) * qJ(4); t54 * MDP(21) + (-MDP(17) * t72 - t98) * t70; -t83 * t70; t79 - t80; t83; t68 * MDP(19) + t50 * MDP(21) + t75 * t70; 0; 0; 0; MDP(21); t70 * MDP(26) + t46 * MDP(27) - t47 * MDP(28) + t78 * t68; t76 * t70; t74; -t76; t77; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
