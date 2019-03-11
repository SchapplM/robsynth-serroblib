% Calculate joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP10_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:34
% EndTime: 2019-03-09 06:32:37
% DurationCPUTime: 1.09s
% Computational Cost: add. (861->200), mult. (1538->277), div. (0->0), fcn. (1502->6), ass. (0->71)
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t112 = t131 * t132 - t134 * t135;
t113 = t131 * t135 + t134 * t132;
t154 = MDP(26) + MDP(28);
t179 = MDP(27) - MDP(30);
t169 = -pkin(9) - pkin(8);
t116 = t169 * t135;
t152 = t169 * t132;
t94 = -t131 * t116 - t134 * t152;
t95 = -t134 * t116 + t131 * t152;
t143 = t113 * MDP(23) - t112 * MDP(24) - t154 * t94 - t179 * t95;
t176 = t132 * MDP(19) + t135 * MDP(20);
t182 = t132 * MDP(16) + t135 * MDP(17) - pkin(8) * t176 + t143;
t136 = cos(qJ(3));
t105 = t113 * t136;
t163 = t135 * t136;
t165 = t132 * t136;
t107 = -t131 * t165 + t134 * t163;
t178 = t107 * MDP(23) - t105 * MDP(24);
t133 = sin(qJ(3));
t115 = t133 * pkin(3) - t136 * pkin(8) + qJ(2);
t110 = t135 * t115;
t137 = -pkin(1) - pkin(7);
t164 = t132 * t137;
t87 = -pkin(9) * t163 + t110 + (pkin(4) - t164) * t133;
t162 = t135 * t137;
t151 = t133 * t162;
t89 = t151 + (-pkin(9) * t136 + t115) * t132;
t82 = t131 * t87 + t134 * t89;
t175 = t133 * MDP(25) - t82 * MDP(27) + t178;
t174 = -2 * MDP(22);
t173 = 0.2e1 * MDP(27);
t172 = 2 * MDP(28);
t171 = 2 * MDP(29);
t170 = 2 * MDP(30);
t168 = (pkin(1) * MDP(6));
t167 = t133 * pkin(5);
t166 = t132 * t135;
t160 = t107 * MDP(21);
t126 = t131 * pkin(4);
t119 = t126 + qJ(6);
t159 = t119 * MDP(30);
t158 = t131 * MDP(27);
t155 = MDP(18) + MDP(25);
t125 = t133 * qJ(6);
t79 = t125 + t82;
t123 = -t135 * pkin(4) - pkin(3);
t149 = MDP(15) * t166;
t148 = t131 * t89 - t134 * t87;
t147 = pkin(5) * t172 + MDP(25);
t111 = pkin(4) * t165 - t136 * t137;
t104 = t113 * t133;
t106 = t112 * t133;
t146 = -t154 * t104 + t179 * t106;
t142 = t135 * MDP(16) - t132 * MDP(17);
t141 = t135 * MDP(19) - t132 * MDP(20);
t139 = (t134 * MDP(26) - t158) * pkin(4);
t130 = t136 ^ 2;
t129 = t135 ^ 2;
t128 = t133 ^ 2;
t127 = t132 ^ 2;
t121 = t134 * pkin(4) + pkin(5);
t97 = t132 * t115 + t151;
t96 = -t133 * t164 + t110;
t86 = t112 * pkin(5) - t113 * qJ(6) + t123;
t83 = t105 * pkin(5) - t107 * qJ(6) + t111;
t80 = t148 - t167;
t1 = [(t79 ^ 2 + t80 ^ 2 + t83 ^ 2) * MDP(31) + MDP(1) + t155 * t128 + (t105 * t174 + t160) * t107 + ((-2 * MDP(4) + t168) * pkin(1)) + (0.2e1 * t136 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t129 * MDP(14) + MDP(7) - 0.2e1 * t149) * t130 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t142) * t136 + t178) * t133 + (-t79 * t105 + t80 * t107) * t171 + (t83 * t105 - t80 * t133) * t172 + (-t83 * t107 + t79 * t133) * t170 + 0.2e1 * (t111 * t105 - t133 * t148) * MDP(26) + (t111 * t107 - t82 * t133) * t173 + 0.2e1 * (-t130 * t164 + t96 * t133) * MDP(19) + 0.2e1 * (-t130 * t162 - t97 * t133) * MDP(20); MDP(4) - t168 + (t104 * t107 + t106 * t105) * MDP(29) + (t80 * t104 - t79 * t106 - t83 * t136) * MDP(31) + t176 * (-t128 - t130) + t154 * (-t104 * t133 - t136 * t105) - t179 * (-t106 * t133 + t136 * t107); MDP(6) + (t104 ^ 2 + t106 ^ 2 + t130) * MDP(31); t113 * t160 + (-t113 * t105 - t107 * t112) * MDP(22) + (t123 * t105 + t111 * t112) * MDP(26) + (t123 * t107 + t111 * t113) * MDP(27) + (t86 * t105 + t83 * t112) * MDP(28) + (-t95 * t105 + t94 * t107 - t79 * t112 + t80 * t113) * MDP(29) + (-t86 * t107 - t83 * t113) * MDP(30) + (t79 * t95 + t80 * t94 + t83 * t86) * MDP(31) + (-t137 * MDP(13) - MDP(10) + t182) * t133 + (MDP(9) + t137 * MDP(12) + MDP(14) * t166 + (-t127 + t129) * MDP(15) + (-pkin(3) * t132 + t162) * MDP(19) + (-pkin(3) * t135 - t164) * MDP(20)) * t136; -t133 * MDP(13) + (t104 * t113 + t106 * t112) * MDP(29) + (t104 * t94 - t106 * t95) * MDP(31) + (-t86 * MDP(31) - t154 * t112 - t113 * t179 + MDP(12) + t141) * t136; MDP(11) + t127 * MDP(14) + 0.2e1 * t149 + (t86 ^ 2 + t94 ^ 2 + t95 ^ 2) * MDP(31) + (MDP(21) * t113 - 0.2e1 * t86 * MDP(30) + t112 * t174 + t123 * t173 + t94 * t171) * t113 + 0.2e1 * t141 * pkin(3) + 0.2e1 * (t123 * MDP(26) + t86 * MDP(28) - t95 * MDP(29)) * t112; t96 * MDP(19) - t97 * MDP(20) + (-t119 * t105 - t121 * t107) * MDP(29) + t79 * MDP(30) + (t79 * t119 - t80 * t121) * MDP(31) + t142 * t136 + (MDP(18) + (pkin(5) + t121) * MDP(28) + t159 + t139) * t133 - t154 * t148 + t175; (-t104 * t121 - t106 * t119) * MDP(31) - t176 * t133 + t146; (-t119 * t112 - t121 * t113) * MDP(29) + (t95 * t119 - t94 * t121) * MDP(31) + t182; (t119 ^ 2 + t121 ^ 2) * MDP(31) + 0.2e1 * t139 + t121 * t172 + 0.2e1 * t159 + t155; -t148 * MDP(26) + (-t148 + 0.2e1 * t167) * MDP(28) + (-pkin(5) * t107 - t105 * qJ(6)) * MDP(29) + (0.2e1 * t125 + t82) * MDP(30) + (-t80 * pkin(5) + t79 * qJ(6)) * MDP(31) + t175; (-t104 * pkin(5) - t106 * qJ(6)) * MDP(31) + t146; (-pkin(5) * t113 - t112 * qJ(6)) * MDP(29) + (-t94 * pkin(5) + t95 * qJ(6)) * MDP(31) + t143; (0.2e1 * qJ(6) + t126) * MDP(30) + (t121 * pkin(5) + t119 * qJ(6)) * MDP(31) + (t154 * t134 - t158) * pkin(4) + t147; qJ(6) * t170 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(31) + t147; -t133 * MDP(28) + t107 * MDP(29) + t80 * MDP(31); t104 * MDP(31); t113 * MDP(29) + t94 * MDP(31); -t121 * MDP(31) - MDP(28); -MDP(31) * pkin(5) - MDP(28); MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
