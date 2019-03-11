% Calculate joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:28
% EndTime: 2019-03-08 21:16:30
% DurationCPUTime: 0.63s
% Computational Cost: add. (495->154), mult. (1058->218), div. (0->0), fcn. (1081->10), ass. (0->76)
t165 = MDP(15) + MDP(19);
t198 = t165 * qJ(4);
t186 = pkin(8) * MDP(15);
t137 = sin(pkin(11));
t139 = cos(pkin(11));
t166 = MDP(13) - MDP(18);
t167 = MDP(12) + MDP(16);
t195 = t167 * t137 + t166 * t139;
t194 = MDP(14) + MDP(17);
t164 = qJ(5) * t137 + pkin(3);
t190 = pkin(4) + pkin(5);
t117 = t139 * t190 + t164;
t193 = 0.2e1 * t117;
t141 = sin(qJ(3));
t192 = 0.2e1 * t141;
t191 = -2 * MDP(21);
t144 = cos(qJ(3));
t189 = pkin(8) * t144;
t188 = -pkin(9) + qJ(4);
t187 = pkin(3) * MDP(15);
t185 = cos(pkin(6));
t183 = t137 * t141;
t143 = cos(qJ(6));
t182 = t137 * t143;
t138 = sin(pkin(6));
t142 = sin(qJ(2));
t181 = t138 * t142;
t145 = cos(qJ(2));
t180 = t138 * t145;
t179 = t139 * t141;
t123 = -pkin(3) * t144 - qJ(4) * t141 - pkin(2);
t110 = t137 * t123 + t139 * t189;
t177 = MDP(24) * t144;
t140 = sin(qJ(6));
t112 = t140 * t179 - t141 * t182;
t176 = t112 * MDP(23);
t119 = t137 * t140 + t139 * t143;
t113 = t119 * t141;
t175 = t113 * MDP(22);
t174 = t119 * MDP(25);
t120 = -t139 * t140 + t182;
t173 = t120 * MDP(20);
t121 = -pkin(4) * t139 - t164;
t172 = t121 * MDP(19);
t171 = t137 * MDP(13);
t170 = t137 * MDP(18);
t169 = t139 * MDP(12);
t168 = t139 * MDP(16);
t163 = qJ(5) * t139 - pkin(8);
t128 = t137 * t189;
t109 = t123 * t139 - t128;
t106 = -qJ(5) * t144 + t110;
t133 = t144 * pkin(4);
t95 = pkin(5) * t144 + t128 + t133 + (-pkin(9) * t141 - t123) * t139;
t97 = pkin(9) * t183 + t106;
t156 = (-t140 * t97 + t143 * t95) * MDP(25) - (t140 * t95 + t143 * t97) * MDP(26);
t107 = -t109 + t133;
t154 = t106 * t139 + t107 * t137;
t153 = -t109 * t137 + t110 * t139;
t152 = t137 * MDP(12) + MDP(13) * t139;
t151 = t137 * MDP(16) - MDP(18) * t139;
t150 = t112 * MDP(25) + t113 * MDP(26);
t149 = MDP(25) * t143 - MDP(26) * t140;
t124 = t188 * t137;
t125 = t188 * t139;
t148 = t120 * MDP(22) - t119 * MDP(23) + (t124 * t143 - t125 * t140) * MDP(25) - (t124 * t140 + t125 * t143) * MDP(26);
t147 = -t120 * MDP(26) + t137 * t166 - t139 * t167 + t172 - t174 - t187;
t116 = t141 * t185 + t144 * t181;
t115 = t141 * t181 - t144 * t185;
t111 = (pkin(4) * t137 - t163) * t141;
t105 = (-t137 * t190 + t163) * t141;
t102 = t116 * t139 - t137 * t180;
t100 = t116 * t137 + t139 * t180;
t93 = t100 * t140 + t102 * t143;
t92 = t100 * t143 - t102 * t140;
t1 = [MDP(1) + t165 * (t100 ^ 2 + t102 ^ 2 + t115 ^ 2); (-t100 * t109 + t102 * t110) * MDP(15) + (t100 * t107 + t102 * t106 + t111 * t115) * MDP(19) + (-t112 * t115 + t144 * t92) * MDP(25) + (-t113 * t115 - t144 * t93) * MDP(26) + t166 * (t102 * t144 + t115 * t179) + (-t142 * MDP(4) + (t144 * MDP(10) + MDP(3)) * t145) * t138 + (-MDP(11) * t180 + t115 * t186 + t194 * (t100 * t139 - t102 * t137)) * t141 + t167 * (t100 * t144 + t115 * t183); MDP(2) + (t109 ^ 2 + t110 ^ 2) * MDP(15) + (t106 ^ 2 + t107 ^ 2 + t111 ^ 2) * MDP(19) + (t113 * MDP(20) + t112 * t191) * t113 + (0.2e1 * pkin(2) * MDP(10) + MDP(6) * t192 + 0.2e1 * t175 - 0.2e1 * t176 + t177) * t144 + 0.2e1 * t150 * t105 + 0.2e1 * (-t109 * MDP(12) + t110 * MDP(13) + t107 * MDP(16) - t106 * MDP(18) + t156) * t144 + ((-t109 * t139 - t110 * t137) * MDP(14) + (-t106 * t137 + t107 * t139) * MDP(17) + t151 * t111) * t192 + (-0.2e1 * pkin(2) * MDP(11) + (MDP(5) + (0.2e1 * t152 + t186) * pkin(8)) * t141) * t141; -t116 * MDP(11) + (-MDP(10) + t147) * t115 + (t194 + t198) * (t100 * t137 + t102 * t139); t153 * MDP(14) + t154 * MDP(17) + t113 * t173 + (-t112 * t120 - t113 * t119) * MDP(21) + (t105 * t119 + t112 * t117) * MDP(25) + (t105 * t120 + t117 * t113) * MDP(26) + (-t168 - t170 + t172) * t111 + (-pkin(8) * MDP(11) + MDP(8) + t148) * t144 + (MDP(7) + t151 * t121 - t152 * pkin(3) + (-MDP(10) - t169 + t171 - t187) * pkin(8)) * t141 + (t153 * MDP(15) + t154 * MDP(19) + t195 * t144) * qJ(4); MDP(9) + t174 * t193 + (-0.2e1 * t168 - 0.2e1 * t170 + t172) * t121 + (0.2e1 * t169 - 0.2e1 * t171 + t187) * pkin(3) + (MDP(26) * t193 + t119 * t191 + t173) * t120 + (0.2e1 * t194 + t198) * (t137 ^ 2 + t139 ^ 2) * qJ(4); t165 * t115; t111 * MDP(19) + (t186 + t195) * t141 - t150; t147; t165; t100 * MDP(19); MDP(17) * t179 + t107 * MDP(19) + (MDP(16) + t149) * t144; (MDP(19) * qJ(4) + MDP(17)) * t137; 0; MDP(19); MDP(25) * t92 - MDP(26) * t93; t156 + t175 - t176 + t177; t148; 0; t149; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
