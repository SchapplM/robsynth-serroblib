% Calculate joint inertia matrix for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP11_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:06
% EndTime: 2019-03-09 17:49:09
% DurationCPUTime: 1.15s
% Computational Cost: add. (1210->241), mult. (2584->334), div. (0->0), fcn. (2645->8), ass. (0->103)
t222 = pkin(4) + pkin(9);
t163 = sin(pkin(6));
t170 = cos(qJ(2));
t210 = t163 * t170;
t221 = 2 * MDP(16);
t220 = 2 * MDP(18);
t219 = 2 * MDP(19);
t218 = 2 * MDP(20);
t217 = pkin(3) + pkin(10);
t167 = sin(qJ(2));
t216 = pkin(1) * t167;
t215 = pkin(1) * t170;
t214 = MDP(29) * pkin(5);
t213 = MDP(30) * pkin(5);
t212 = (pkin(3) * MDP(21));
t211 = t163 * t167;
t164 = cos(pkin(6));
t209 = t164 * MDP(8);
t166 = sin(qJ(3));
t169 = cos(qJ(3));
t140 = -t164 * t169 + t166 * t211;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t128 = t140 * t168 + t165 * t210;
t208 = t165 * t128;
t149 = t222 * t166;
t207 = t165 * t149;
t206 = t165 * t168;
t205 = t167 * MDP(6);
t129 = -t140 * t165 + t168 * t210;
t204 = t168 * t129;
t203 = -qJ(6) - t217;
t190 = pkin(8) * t210;
t135 = t190 + (pkin(9) + t216) * t164;
t136 = (-pkin(2) * t170 - pkin(9) * t167 - pkin(1)) * t163;
t122 = t169 * t135 + t166 * t136;
t150 = t222 * t169;
t202 = MDP(25) * t168;
t121 = -t166 * t135 + t169 * t136;
t153 = pkin(3) * t210;
t120 = -t121 + t153;
t201 = t120 * MDP(21);
t200 = t128 * MDP(25);
t199 = t129 * MDP(22);
t198 = t140 * MDP(12);
t197 = t165 * MDP(27);
t196 = t165 * MDP(28);
t195 = t168 * MDP(28);
t194 = t170 * MDP(15);
t193 = pkin(9) ^ 2 * MDP(21);
t192 = MDP(11) + MDP(26);
t191 = MDP(17) - MDP(20);
t189 = qJ(4) * t210;
t188 = MDP(23) * t206;
t187 = -t166 * qJ(4) - pkin(2);
t186 = MDP(19) - t212;
t144 = -t169 * t217 + t187;
t185 = qJ(6) * t169 - t144;
t141 = t164 * t166 + t169 * t211;
t115 = t141 * pkin(4) + pkin(10) * t210 + t120;
t152 = pkin(8) * t211;
t134 = t152 + (-pkin(2) - t215) * t164;
t176 = -t141 * qJ(4) + t134;
t116 = t140 * t217 + t176;
t111 = t168 * t115 - t165 * t116;
t184 = -MDP(27) * t217 + MDP(24);
t109 = t141 * pkin(5) + t129 * qJ(6) + t111;
t112 = t165 * t115 + t168 * t116;
t110 = t128 * qJ(6) + t112;
t183 = t109 * t168 + t110 * t165;
t145 = t168 * t149;
t123 = t166 * pkin(5) + t165 * t185 + t145;
t124 = -t168 * t185 + t207;
t182 = t123 * t168 + t124 * t165;
t181 = (MDP(28) * t217 - MDP(25)) * t165;
t125 = -t165 * t144 + t145;
t126 = t168 * t144 + t207;
t180 = t125 * MDP(27) - t126 * MDP(28);
t179 = t168 * MDP(27) - t196;
t119 = t189 - t122;
t178 = -MDP(24) * t165 + MDP(12) - t202;
t177 = MDP(18) + t179;
t117 = -t140 * pkin(4) - t119;
t175 = pkin(9) * MDP(21) + t177;
t174 = -t129 * MDP(24) + t111 * MDP(27) - t112 * MDP(28) + t200;
t173 = -pkin(3) * MDP(18) + t168 * t184 + MDP(13) + t181;
t162 = t169 ^ 2;
t161 = t168 ^ 2;
t160 = t166 ^ 2;
t159 = t165 ^ 2;
t158 = t163 ^ 2;
t155 = t165 * pkin(5) + qJ(4);
t151 = t159 + t161;
t148 = -t169 * pkin(3) + t187;
t147 = t203 * t168;
t146 = t203 * t165;
t143 = t164 * t216 + t190;
t142 = t164 * t215 - t152;
t139 = t168 * t169 * pkin(5) + t150;
t127 = t146 * t165 + t147 * t168;
t118 = t140 * pkin(3) + t176;
t113 = -t128 * pkin(5) + t117;
t1 = [(t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(30) + (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) * MDP(21) + t158 * t167 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t163 * t205 + t209) * t164 + t192 * t141 ^ 2 + (-0.2e1 * t128 * MDP(23) - 0.2e1 * t141 * MDP(24) + t199) * t129 + (0.2e1 * MDP(5) * t167 + t194) * t158 * t170 + (-t118 * t141 + t119 * t210) * t218 + (-t121 * t210 + t134 * t140) * t221 + (-t118 * t140 - t120 * t210) * t219 + 0.2e1 * (t122 * t210 + t134 * t141) * MDP(17) + 0.2e1 * (t142 * t164 + t158 * t215) * MDP(9) + 0.2e1 * (-t143 * t164 - t158 * t216) * MDP(10) + (t119 * t140 + t120 * t141) * t220 + 0.2e1 * (t111 * t141 - t117 * t128) * MDP(27) + 0.2e1 * (-t112 * t141 - t117 * t129) * MDP(28) + 0.2e1 * (t109 * t129 + t110 * t128) * MDP(29) + 0.2e1 * (-t198 + t200) * t141 + 0.2e1 * (-t141 * MDP(13) + t140 * MDP(14) + MDP(7) * t164) * t210; (t109 * t123 + t110 * t124 + t113 * t139) * MDP(30) + t209 + t142 * MDP(9) - t143 * MDP(10) + (t125 * t141 - t150 * t128) * MDP(27) + (-t126 * t141 - t150 * t129) * MDP(28) + (t123 * t129 + t124 * t128) * MDP(29) + (t170 * MDP(7) + t205) * t163 + (-t140 * MDP(16) - t141 * MDP(17)) * pkin(2) + (-t140 * MDP(19) - t141 * MDP(20) + t118 * MDP(21)) * t148 + (-MDP(13) * t210 - t198 + t134 * MDP(17) + t120 * MDP(18) - t118 * MDP(20) + t192 * t141 + (t141 * MDP(18) + t201 + (MDP(16) - MDP(19)) * t210) * pkin(9) + t174) * t166 + (-MDP(14) * t210 + t165 * t199 - t134 * MDP(16) - t119 * MDP(18) + t118 * MDP(19) + (t204 - t208) * MDP(23) + (t109 * t165 - t110 * t168) * MDP(29) + t179 * t117 + t178 * t141 + (-t140 * MDP(18) - t119 * MDP(21) + t191 * t210) * pkin(9)) * t169; MDP(8) + pkin(2) * t169 * t221 + (t123 ^ 2 + t124 ^ 2 + t139 ^ 2) * MDP(30) + (MDP(21) * t148 + t169 * t219) * t148 + (t159 * MDP(22) + 0.2e1 * t188 + t193) * t162 + (t192 + t193) * t160 + 0.2e1 * ((t123 * t165 - t124 * t168) * MDP(29) + t179 * t150) * t169 + (t160 + t162) * pkin(9) * t220 + 0.2e1 * (-pkin(2) * MDP(17) - t148 * MDP(20) + t178 * t169 + t180) * t166; -t163 * t194 + t121 * MDP(16) - t122 * MDP(17) + (-t121 + 0.2e1 * t153) * MDP(19) + (-0.2e1 * t189 + t122) * MDP(20) + (-t120 * pkin(3) - t119 * qJ(4)) * MDP(21) - t168 * t199 + (t168 * t128 + t129 * t165) * MDP(23) + (-qJ(4) * t128 + t117 * t165) * MDP(27) + (-qJ(4) * t129 + t117 * t168) * MDP(28) + (t146 * t128 + t147 * t129 - t183) * MDP(29) + (t109 * t147 + t110 * t146 + t113 * t155) * MDP(30) + (-qJ(4) * MDP(18) - MDP(14)) * t140 + t173 * t141; -t182 * MDP(29) + (t123 * t147 + t124 * t146 + t139 * t155) * MDP(30) + (t195 + t197) * t150 + ((-MDP(16) + t186) * pkin(9) + t173) * t166 + (MDP(14) - MDP(22) * t206 + (t159 - t161) * MDP(23) + (-t146 * t168 + t147 * t165) * MDP(29) - t191 * pkin(9) + t175 * qJ(4)) * t169; MDP(15) + t161 * MDP(22) - 0.2e1 * t188 - 0.2e1 * t127 * MDP(29) + (t146 ^ 2 + t147 ^ 2 + t155 ^ 2) * MDP(30) + (-2 * MDP(19) + t212) * pkin(3) + (MDP(21) * qJ(4) + 0.2e1 * t195 + 0.2e1 * t197 + t218) * qJ(4); -MDP(19) * t210 + t201 + (t204 + t208) * MDP(29) + t183 * MDP(30) + t177 * t141; MDP(30) * t182 + t166 * t175; -t151 * MDP(29) + t127 * MDP(30) + t186; t151 * MDP(30) + MDP(21); t141 * MDP(26) + (t129 * MDP(29) + t109 * MDP(30)) * pkin(5) + t174; t123 * t213 + t166 * MDP(26) + (-t202 + (-MDP(24) + t214) * t165) * t169 + t180; t147 * t213 + t181 + (t184 - t214) * t168; -t196 + (MDP(27) + t213) * t168; MDP(30) * pkin(5) ^ 2 + MDP(26); t113 * MDP(30); t139 * MDP(30); t155 * MDP(30); 0; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
