% Calculate joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRP8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:25
% EndTime: 2019-03-09 02:16:27
% DurationCPUTime: 0.44s
% Computational Cost: add. (650->129), mult. (1076->178), div. (0->0), fcn. (1114->6), ass. (0->59)
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t104 = sin(qJ(4));
t135 = cos(qJ(4));
t85 = t104 * t100 - t135 * t101;
t139 = t85 ^ 2;
t86 = t135 * t100 + t104 * t101;
t82 = t86 ^ 2;
t145 = t82 + t139;
t130 = MDP(28) * t86;
t144 = MDP(23) + MDP(25);
t126 = MDP(24) - MDP(27);
t103 = sin(qJ(5));
t105 = cos(qJ(5));
t140 = -t126 * t103 + t144 * t105 + MDP(16);
t137 = pkin(5) * t86;
t136 = pkin(8) * t86;
t102 = -pkin(1) - qJ(3);
t134 = -pkin(7) + t102;
t92 = t100 * pkin(3) + qJ(2);
t77 = pkin(4) * t86 + pkin(8) * t85 + t92;
t88 = t134 * t100;
t89 = t134 * t101;
t80 = t104 * t89 + t135 * t88;
t73 = t103 * t77 + t105 * t80;
t91 = t100 ^ 2 + t101 ^ 2;
t98 = t103 ^ 2;
t99 = t105 ^ 2;
t133 = t98 + t99;
t132 = qJ(6) * t86;
t131 = t105 * t85;
t117 = pkin(5) * t105 + qJ(6) * t103;
t90 = -pkin(4) - t117;
t129 = t90 * MDP(28);
t128 = MDP(19) * t105;
t125 = t103 * t80 - t105 * t77;
t124 = t133 * MDP(26);
t123 = t133 * MDP(28);
t122 = -MDP(28) * pkin(5) - MDP(25);
t121 = MDP(28) * pkin(8) + MDP(26);
t120 = pkin(4) * t85 - t136;
t119 = t85 * t90 + t136;
t118 = MDP(23) - t122;
t116 = -pkin(5) * t103 + qJ(6) * t105;
t115 = MDP(28) * qJ(6) - t126;
t114 = -MDP(17) + t124;
t113 = -t125 * MDP(23) - t73 * MDP(24);
t79 = t104 * t88 - t135 * t89;
t74 = t116 * t85 + t79;
t112 = -t79 * MDP(23) - t74 * MDP(25);
t111 = t79 * MDP(24) - t74 * MDP(27);
t110 = t100 * MDP(7) + t101 * MDP(8);
t109 = -t105 * MDP(20) + t103 * MDP(21);
t108 = -t118 * t103 + t115 * t105;
t106 = qJ(2) ^ 2;
t84 = t91 * t102;
t71 = t125 - t137;
t70 = t132 + t73;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t106) * MDP(6) + (t91 * t102 ^ 2 + t106) * MDP(10) + t82 * MDP(22) + (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) * MDP(28) + (MDP(18) * t99 - 0.2e1 * t103 * t128 + MDP(11)) * t139 - 0.2e1 * t92 * t85 * MDP(17) - 0.2e1 * t84 * MDP(9) + 0.2e1 * ((-t71 * MDP(26) - t111) * t105 + (t70 * MDP(26) + t112) * t103) * t85 + 0.2e1 * (MDP(5) + t110) * qJ(2) + 0.2e1 * ((MDP(12) + t109) * t85 + t92 * MDP(16) - MDP(25) * t71 + MDP(27) * t70 + t113) * t86; MDP(4) - pkin(1) * MDP(6) - t91 * MDP(9) + t84 * MDP(10) + t74 * t85 * MDP(28) + (-t126 * t145 + t70 * t130) * t105 + (t71 * t130 - t144 * t145) * t103; MDP(6) + t91 * MDP(10) + (t133 * t82 + t139) * MDP(28); qJ(2) * MDP(10) + (t103 * t70 - t105 * t71) * MDP(28) + t114 * t85 + t140 * t86 + t110; 0; MDP(10) + t123; t74 * t129 - t86 * MDP(14) - t79 * MDP(16) - t80 * MDP(17) + (-MDP(13) + (t98 - t99) * MDP(19)) * t85 + (t86 * MDP(21) + t120 * MDP(24) + t119 * MDP(27) + t121 * t70 + t112) * t105 + (-MDP(18) * t131 + t86 * MDP(20) + t120 * MDP(23) - t119 * MDP(25) + t121 * t71 + t111) * t103; (pkin(8) * t123 + t114) * t86 + (t129 - t140) * t85; 0; MDP(15) + t98 * MDP(18) + (t133 * pkin(8) ^ 2 + t90 ^ 2) * MDP(28) + 0.2e1 * pkin(8) * t124 + 0.2e1 * (MDP(23) * pkin(4) - MDP(25) * t90) * t105 + 0.2e1 * (-MDP(24) * pkin(4) - MDP(27) * t90 + t128) * t103; t86 * MDP(22) + (-t125 + 0.2e1 * t137) * MDP(25) + (0.2e1 * t132 + t73) * MDP(27) + (-pkin(5) * t71 + qJ(6) * t70) * MDP(28) + (t117 * MDP(26) + t109) * t85 + t113; t108 * t86; t115 * t103 + t118 * t105; t103 * MDP(20) + t105 * MDP(21) + t116 * MDP(26) + t108 * pkin(8); MDP(22) + 0.2e1 * pkin(5) * MDP(25) + 0.2e1 * qJ(6) * MDP(27) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(28); -MDP(25) * t86 - MDP(26) * t131 + MDP(28) * t71; t103 * t130; -t105 * MDP(28); t121 * t103; t122; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
