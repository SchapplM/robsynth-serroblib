% Calculate joint inertia matrix for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:40
% EndTime: 2019-03-08 19:37:41
% DurationCPUTime: 0.25s
% Computational Cost: add. (236->98), mult. (509->142), div. (0->0), fcn. (501->10), ass. (0->56)
t80 = sin(pkin(11));
t81 = sin(pkin(6));
t82 = cos(pkin(11));
t86 = sin(qJ(2));
t89 = cos(qJ(2));
t65 = (t80 * t86 - t82 * t89) * t81;
t118 = t65 ^ 2;
t90 = -pkin(4) - pkin(9);
t74 = pkin(2) * t80 + pkin(8);
t117 = pkin(5) + t74;
t116 = MDP(5) * pkin(2);
t84 = sin(qJ(6));
t87 = cos(qJ(6));
t115 = t84 * t87;
t88 = cos(qJ(4));
t114 = t84 * t88;
t113 = t87 * t88;
t85 = sin(qJ(4));
t77 = t85 ^ 2;
t79 = t88 ^ 2;
t112 = t77 + t79;
t111 = MDP(16) * pkin(4);
t110 = MDP(16) * t74;
t109 = MDP(22) * t84;
t108 = MDP(23) * t87;
t75 = -pkin(2) * t82 - pkin(3);
t95 = -qJ(5) * t85 + t75;
t69 = -pkin(4) * t88 + t95;
t107 = t69 * MDP(16);
t106 = t74 ^ 2 * MDP(16);
t105 = MDP(16) * qJ(5);
t104 = MDP(12) - MDP(15);
t103 = MDP(18) * t115;
t102 = MDP(14) - t111;
t100 = MDP(11) - t102;
t99 = -t104 + t105;
t98 = -MDP(19) * t84 - MDP(20) * t87;
t97 = MDP(22) * t87 - t84 * MDP(23);
t96 = t108 + t109;
t94 = MDP(13) + t97;
t93 = t96 + t99;
t92 = (MDP(22) * t90 + MDP(19)) * t87 + (-MDP(23) * t90 - MDP(20)) * t84;
t83 = cos(pkin(6));
t78 = t87 ^ 2;
t76 = t84 ^ 2;
t71 = t117 * t88;
t70 = t117 * t85;
t68 = t90 * t88 + t95;
t67 = (t80 * t89 + t82 * t86) * t81;
t63 = t67 * t88 + t83 * t85;
t62 = t67 * t85 - t83 * t88;
t61 = t68 * t87 + t70 * t84;
t60 = -t68 * t84 + t70 * t87;
t59 = t62 * t84 + t65 * t87;
t58 = t62 * t87 - t65 * t84;
t1 = [MDP(1) + (t67 ^ 2 + t83 ^ 2 + t118) * MDP(5) + (t62 ^ 2 + t63 ^ 2 + t118) * MDP(16); t67 * t80 * t116 + (t63 * t113 + t58 * t85) * MDP(22) + (-t63 * t114 - t59 * t85) * MDP(23) + (t89 * MDP(3) - t86 * MDP(4)) * t81 + (-t82 * t116 + t107 + (-MDP(11) + MDP(14)) * t88 + t104 * t85) * t65 + (MDP(13) + t110) * (t62 * t85 + t63 * t88); -0.2e1 * t75 * t88 * MDP(11) + MDP(2) + (t80 ^ 2 + t82 ^ 2) * MDP(5) * pkin(2) ^ 2 + (0.2e1 * t88 * MDP(14) + t107) * t69 + (t76 * MDP(17) + 0.2e1 * t103 + t106) * t79 + (MDP(21) + MDP(6) + t106) * t77 + 0.2e1 * (t75 * MDP(12) - t69 * MDP(15) + (MDP(7) + t98) * t88) * t85 + 0.2e1 * (t71 * t113 + t60 * t85) * MDP(22) + 0.2e1 * (-t71 * t114 - t61 * t85) * MDP(23) + 0.2e1 * t112 * MDP(13) * t74; t83 * MDP(5) + (-t62 * t88 + t63 * t85) * MDP(16); 0; t112 * MDP(16) + MDP(5); -t100 * t62 + t93 * t63; t96 * t71 + (-pkin(4) * MDP(13) + MDP(8) + t92) * t85 + (MDP(9) - MDP(17) * t115 + (t76 - t78) * MDP(18) + t94 * qJ(5)) * t88 + (-t100 * t85 + t99 * t88) * t74; t100 * t88 + t93 * t85; -0.2e1 * t103 + t78 * MDP(17) + MDP(10) + (-0.2e1 * MDP(14) + t111) * pkin(4) + (0.2e1 * MDP(15) + t105 + 0.2e1 * t108 + 0.2e1 * t109) * qJ(5); t62 * MDP(16); (t94 + t110) * t85; -t88 * MDP(16); t102; MDP(16); t58 * MDP(22) - t59 * MDP(23); MDP(21) * t85 + t60 * MDP(22) - t61 * MDP(23) + t98 * t88; -t97 * t88; t92; t97; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
