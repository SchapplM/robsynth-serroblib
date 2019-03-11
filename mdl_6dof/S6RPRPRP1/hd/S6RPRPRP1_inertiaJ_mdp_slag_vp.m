% Calculate joint inertia matrix for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6RPRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:14
% EndTime: 2019-03-09 03:03:14
% DurationCPUTime: 0.35s
% Computational Cost: add. (558->116), mult. (974->177), div. (0->0), fcn. (1001->8), ass. (0->54)
t92 = sin(pkin(9));
t86 = pkin(1) * t92 + pkin(7);
t122 = qJ(4) + t86;
t121 = cos(qJ(3));
t91 = sin(pkin(10));
t93 = cos(pkin(10));
t96 = sin(qJ(3));
t79 = -t93 * t121 + t91 * t96;
t77 = t79 ^ 2;
t108 = t122 * t96;
t76 = t122 * t121;
t70 = -t91 * t108 + t93 * t76;
t97 = cos(qJ(5));
t120 = t70 * t97;
t81 = t91 * t121 + t93 * t96;
t95 = sin(qJ(5));
t89 = t95 ^ 2;
t119 = t81 * t89;
t90 = t97 ^ 2;
t118 = t89 + t90;
t117 = MDP(22) * pkin(5);
t116 = qJ(6) * t81;
t85 = pkin(3) * t91 + pkin(8);
t115 = qJ(6) + t85;
t114 = MDP(17) * t95;
t73 = t90 * t81;
t113 = (-t73 - t119) * MDP(21);
t87 = -pkin(3) * t93 - pkin(4);
t82 = -pkin(5) * t97 + t87;
t112 = t82 * MDP(22);
t111 = t95 * MDP(20);
t110 = t97 * MDP(20);
t109 = t95 * t97 * MDP(15);
t94 = cos(pkin(9));
t88 = -t94 * pkin(1) - pkin(2);
t83 = -t121 * pkin(3) + t88;
t67 = t79 * pkin(4) - t81 * pkin(8) + t83;
t63 = t97 * t67 - t70 * t95;
t68 = t93 * t108 + t76 * t91;
t107 = t121 * MDP(10);
t106 = -MDP(21) * pkin(5) + MDP(16);
t105 = MDP(19) + t117;
t61 = pkin(5) * t79 - t97 * t116 + t63;
t62 = t120 + (t67 - t116) * t95;
t104 = t61 * t97 + t62 * t95;
t74 = t115 * t95;
t75 = t115 * t97;
t103 = t74 * t95 + t75 * t97;
t102 = -t79 * t85 + t81 * t87;
t101 = t63 * MDP(19) - (t67 * t95 + t120) * MDP(20);
t100 = MDP(19) * t97 - t111;
t78 = t81 ^ 2;
t65 = pkin(5) * t81 * t95 + t68;
t1 = [MDP(1) - 0.2e1 * t88 * t107 + (t68 ^ 2 + t70 ^ 2 + t83 ^ 2) * MDP(13) + t77 * MDP(18) + (t61 ^ 2 + t62 ^ 2 + t65 ^ 2) * MDP(22) + (t92 ^ 2 + t94 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t90 * MDP(14) - 0.2e1 * t109) * t78 + (0.2e1 * t88 * MDP(11) + MDP(5) * t96 + 0.2e1 * t121 * MDP(6)) * t96 + 0.2e1 * (-t70 * MDP(12) + t101) * t79 + 0.2e1 * ((MDP(16) * t97 - t114) * t79 - t104 * MDP(21) + (t95 * MDP(19) + MDP(12) + t110) * t68) * t81; (t68 * t79 + t70 * t81) * MDP(13) + (t65 * t79 + (-t61 * t95 + t62 * t97) * t81) * MDP(22); MDP(4) + (t78 + t77) * MDP(13) + (t118 * t78 + t77) * MDP(22); t96 * MDP(7) + t121 * MDP(8) + (t73 - t119) * MDP(15) + (-t61 * t74 + t62 * t75 + t65 * t82) * MDP(22) + (-t96 * MDP(10) - t121 * MDP(11)) * t86 + (t79 * MDP(17) - t68 * MDP(19) + t102 * MDP(20) + (t74 * t81 + t62) * MDP(21)) * t97 + (t97 * t81 * MDP(14) + t79 * MDP(16) + t102 * MDP(19) + t68 * MDP(20) + (-t75 * t81 - t61) * MDP(21)) * t95 + ((-t79 * t91 - t81 * t93) * MDP(12) + (-t68 * t93 + t70 * t91) * MDP(13)) * pkin(3); t107 - t96 * MDP(11) - t113 + t103 * MDP(22) * t81 + (-t100 + t112) * t79 + (-t79 * t93 + t81 * t91) * MDP(13) * pkin(3); MDP(9) + t89 * MDP(14) + 0.2e1 * t109 + 0.2e1 * t103 * MDP(21) + (t74 ^ 2 + t75 ^ 2 + t82 ^ 2) * MDP(22) + (t91 ^ 2 + t93 ^ 2) * MDP(13) * pkin(3) ^ 2 - 0.2e1 * t100 * t87; t83 * MDP(13) + t104 * MDP(22) + t100 * t79 + t113; 0; (-t74 * t97 + t75 * t95) * MDP(22); t118 * MDP(22) + MDP(13); t61 * t117 + MDP(18) * t79 + (t106 * t97 - t114) * t81 + t101; (-t105 * t95 - t110) * t81; -t74 * t117 + (-MDP(20) * t85 + MDP(17)) * t97 + (-MDP(19) * t85 + t106) * t95; t105 * t97 - t111; MDP(22) * pkin(5) ^ 2 + MDP(18); t65 * MDP(22); t79 * MDP(22); t112; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
