% Calculate joint inertia matrix for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP13_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:35
% EndTime: 2019-03-09 13:01:37
% DurationCPUTime: 1.22s
% Computational Cost: add. (1062->234), mult. (2273->322), div. (0->0), fcn. (2290->8), ass. (0->97)
t147 = sin(pkin(6));
t154 = cos(qJ(2));
t191 = t147 * t154;
t151 = sin(qJ(2));
t192 = t147 * t151;
t135 = pkin(8) * t192;
t148 = cos(pkin(6));
t201 = pkin(1) * t154;
t170 = -pkin(2) - t201;
t113 = pkin(3) * t192 + t135 + (-pkin(9) + t170) * t148;
t155 = -pkin(2) - pkin(9);
t168 = -qJ(3) * t151 - pkin(1);
t119 = (t154 * t155 + t168) * t147;
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t108 = t153 * t113 - t150 * t119;
t109 = t150 * t113 + t153 * t119;
t126 = t148 * t153 - t150 * t191;
t208 = t126 * MDP(17) + t108 * MDP(20) - t109 * MDP(21);
t207 = 0.2e1 * t148;
t206 = 0.2e1 * t150;
t204 = 2 * MDP(27);
t203 = 2 * MDP(28);
t202 = 2 * MDP(29);
t200 = -qJ(6) - pkin(10);
t199 = (MDP(30) * pkin(5));
t198 = pkin(2) * MDP(14);
t197 = qJ(6) * t153;
t106 = -pkin(4) * t192 - t108;
t149 = sin(qJ(5));
t196 = t106 * t149;
t152 = cos(qJ(5));
t195 = t106 * t152;
t116 = t126 * t149 - t152 * t192;
t194 = t116 * t152;
t117 = t126 * t152 + t149 * t192;
t193 = t117 * t149;
t190 = t149 * t155;
t189 = t152 * t155;
t128 = t148 * t151 * pkin(1) + pkin(8) * t191;
t142 = t149 ^ 2;
t145 = t152 ^ 2;
t188 = t142 + t145;
t143 = t150 ^ 2;
t146 = t153 ^ 2;
t187 = -t143 - t146;
t105 = t116 * pkin(5) + t106;
t186 = t105 * MDP(30);
t131 = t150 * pkin(4) - t153 * pkin(10) + qJ(3);
t129 = t152 * t131;
t112 = -t152 * t197 + t129 + (pkin(5) - t190) * t150;
t185 = t112 * MDP(30);
t184 = t116 * MDP(25);
t183 = t117 * MDP(22);
t182 = t117 * MDP(24);
t125 = t148 * t150 + t153 * t191;
t181 = t125 * MDP(26);
t180 = t126 * MDP(16);
t178 = t126 * MDP(21);
t139 = -t152 * pkin(5) - pkin(4);
t177 = t139 * MDP(30);
t176 = t149 * MDP(25);
t175 = t152 * MDP(22);
t174 = t152 * MDP(28);
t173 = t153 * MDP(30);
t172 = t155 * MDP(21);
t171 = t150 * t189;
t140 = t148 * qJ(3);
t122 = -t140 - t128;
t169 = t149 * t152 * MDP(23);
t167 = -MDP(29) * pkin(5) + MDP(24);
t107 = pkin(10) * t192 + t109;
t118 = pkin(3) * t191 - t122;
t111 = t125 * pkin(4) - t126 * pkin(10) + t118;
t103 = -t149 * t107 + t152 * t111;
t166 = t155 * MDP(20) + MDP(17);
t101 = t125 * pkin(5) - t117 * qJ(6) + t103;
t104 = t152 * t107 + t149 * t111;
t102 = -t116 * qJ(6) + t104;
t165 = -t101 * t149 + t102 * t152;
t132 = t200 * t149;
t133 = t200 * t152;
t164 = -t132 * t149 - t133 * t152;
t163 = -t128 * MDP(10) + (t148 * t201 - t135) * MDP(9);
t120 = -t150 * t190 + t129;
t121 = t149 * t131 + t171;
t161 = t120 * MDP(27) - t121 * MDP(28);
t160 = t152 * MDP(27) - t149 * MDP(28);
t159 = -t149 * MDP(27) - t174;
t158 = t152 * MDP(24) - MDP(16) - t176;
t157 = t103 * MDP(27) - t104 * MDP(28) + t181 + t182 - t184;
t156 = t149 * MDP(24) + t152 * MDP(25) + pkin(10) * t159 - MDP(18);
t130 = (pkin(5) * t149 - t155) * t153;
t124 = t148 * t170 + t135;
t123 = (-pkin(2) * t154 + t168) * t147;
t115 = t171 + (t131 - t197) * t149;
t1 = [t148 ^ 2 * MDP(8) + t126 ^ 2 * MDP(15) + (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) * MDP(14) + (t101 ^ 2 + t102 ^ 2 + t105 ^ 2) * MDP(30) + MDP(1) + (-0.2e1 * t116 * MDP(23) + t183) * t117 + (-0.2e1 * MDP(18) * t192 - 0.2e1 * t180 + t181 + 0.2e1 * t182 - 0.2e1 * t184) * t125 + (-t101 * t117 - t102 * t116) * t202 + (t103 * t125 + t106 * t116) * t204 + (-t104 * t125 + t106 * t117) * t203 + (t124 * MDP(12) - t122 * MDP(13) + t163) * t207 + 0.2e1 * (t125 * MDP(20) + t178) * t118 + ((t151 * MDP(6) + t154 * MDP(7)) * t207 + 0.2e1 * (-t122 * MDP(11) + t123 * MDP(12)) * t154 + ((MDP(19) + MDP(4)) * t151 ^ 2 + 0.2e1 * (-MDP(10) * t151 + MDP(9) * t154) * pkin(1)) * t147 + 0.2e1 * (t124 * MDP(11) - t123 * MDP(13) + MDP(5) * t191 + t208) * t151) * t147; t135 * MDP(12) + (0.2e1 * t140 + t128) * MDP(13) + (-t124 * pkin(2) - t122 * qJ(3)) * MDP(14) + qJ(3) * t178 + (-t112 * t117 - t115 * t116) * MDP(29) + (t101 * t112 + t102 * t115 + t105 * t130) * MDP(30) + (MDP(8) + (-0.2e1 * pkin(2) - t201) * MDP(12)) * t148 + (qJ(3) * MDP(20) + t161) * t125 + (t118 * MDP(20) + t157 - t180) * t150 + ((qJ(3) * MDP(11) + MDP(7)) * t154 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(18) - t172) * t150) * t151) * t147 + (t126 * MDP(15) + t118 * MDP(21) + (-t193 - t194) * MDP(23) + (-t116 * t155 + t196) * MDP(27) + (-t117 * t155 + t195) * MDP(28) + (-t101 * t152 - t102 * t149) * MDP(29) + t117 * t175 + t166 * t192 + t158 * t125) * t153 + t163; MDP(8) + t143 * MDP(26) + (t112 ^ 2 + t115 ^ 2 + t130 ^ 2) * MDP(30) + (-0.2e1 * MDP(12) + t198) * pkin(2) + (MDP(14) * qJ(3) + MDP(20) * t206 + 0.2e1 * MDP(13)) * qJ(3) + (t145 * MDP(22) + MDP(15) - 0.2e1 * t169) * t146 + (t120 * t150 - t146 * t190) * t204 + (-t121 * t150 - t146 * t189) * t203 + (0.2e1 * qJ(3) * MDP(21) + (-t112 * t152 - t115 * t149) * t202 + t158 * t206) * t153; MDP(11) * t192 + t148 * MDP(12) + t124 * MDP(14) + (MDP(20) * t192 - t116 * MDP(27) - t117 * MDP(28) - t186) * t153 + (-MDP(21) * t192 + (t193 - t194) * MDP(29) + t165 * MDP(30) + t159 * t125) * t150; -t130 * t173 - t198 + MDP(12) + (t115 * t150 * MDP(30) + MDP(28) * t187) * t152 + (MDP(27) * t187 - t150 * t185) * t149; MDP(14) + (t143 * t188 + t146) * MDP(30); MDP(19) * t192 + t149 * t183 + (-t149 * t116 + t117 * t152) * MDP(23) + (-pkin(4) * t116 - t195) * MDP(27) + (-pkin(4) * t117 + t196) * MDP(28) + (t133 * t116 - t132 * t117 + t165) * MDP(29) + (t101 * t132 - t102 * t133 + t105 * t139) * MDP(30) + t156 * t125 + t208; (-t112 * t149 + t115 * t152) * MDP(29) + (t112 * t132 - t115 * t133 + t130 * t139) * MDP(30) + (t156 - t172) * t150 + (t149 * t175 + (-t142 + t145) * MDP(23) + (-pkin(4) * t149 + t189) * MDP(27) + (-pkin(4) * t152 - t190) * MDP(28) + (-t132 * t152 + t133 * t149) * MDP(29) + t166) * t153; (MDP(20) + t160 - t177) * t153 + (MDP(29) * t188 + MDP(30) * t164 - MDP(21)) * t150; MDP(19) + t142 * MDP(22) + 0.2e1 * t169 + t164 * t202 + (t132 ^ 2 + t133 ^ 2 + t139 ^ 2) * MDP(30) + 0.2e1 * t160 * pkin(4); (-t117 * MDP(29) + t101 * MDP(30)) * pkin(5) + t157; pkin(5) * t185 + t150 * MDP(26) + (t152 * t167 - t176) * t153 + t161; (-t174 + (-MDP(27) - t199) * t149) * t150; t132 * t199 + (-MDP(28) * pkin(10) + MDP(25)) * t152 + (-MDP(27) * pkin(10) + t167) * t149; MDP(30) * pkin(5) ^ 2 + MDP(26); t186; t130 * MDP(30); -t173; t177; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
