% Calculate joint inertia matrix for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR13_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:25:24
% EndTime: 2019-03-09 04:25:28
% DurationCPUTime: 1.23s
% Computational Cost: add. (1702->226), mult. (4532->333), div. (0->0), fcn. (5130->12), ass. (0->108)
t176 = cos(pkin(7));
t173 = sin(pkin(7));
t235 = cos(pkin(6));
t199 = t235 * t173;
t174 = sin(pkin(6));
t175 = cos(pkin(12));
t227 = t174 * t175;
t245 = t176 * t227 + t199;
t177 = sin(qJ(6));
t180 = cos(qJ(6));
t244 = t177 * MDP(31) + t180 * MDP(32);
t172 = sin(pkin(12));
t201 = pkin(1) * t235;
t156 = qJ(2) * t227 + t172 * t201;
t140 = pkin(9) * t245 + t156;
t164 = t175 * t201;
t230 = t172 * t174;
t144 = t235 * pkin(2) + t164 + (-pkin(9) * t176 - qJ(2)) * t230;
t149 = (-pkin(9) * t172 * t173 - pkin(2) * t175 - pkin(1)) * t174;
t179 = sin(qJ(3));
t182 = cos(qJ(3));
t128 = -t179 * t140 + (t144 * t176 + t149 * t173) * t182;
t243 = 2 * MDP(17);
t242 = 2 * MDP(25);
t241 = 0.2e1 * MDP(31);
t240 = 0.2e1 * MDP(32);
t239 = pkin(3) + pkin(10);
t167 = t174 ^ 2;
t238 = pkin(1) * t167;
t153 = t173 * t227 - t235 * t176;
t237 = t153 * pkin(3);
t236 = pkin(3) * MDP(18);
t152 = qJ(4) * t153;
t226 = t176 * t179;
t143 = t179 * t199 + (t172 * t182 + t175 * t226) * t174;
t122 = t143 * pkin(4) + t239 * t153 - t128;
t142 = t179 * t230 - t182 * t245;
t132 = -t173 * t144 + t176 * t149;
t190 = -t143 * qJ(4) + t132;
t123 = t239 * t142 + t190;
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t119 = t181 * t122 - t178 * t123;
t117 = -t143 * pkin(5) - t119;
t234 = t117 * t177;
t233 = t117 * t180;
t133 = -t142 * t181 - t153 * t178;
t232 = t133 * t178;
t228 = t173 * t182;
t157 = t178 * t176 + t181 * t228;
t231 = t157 * t181;
t229 = t173 * t179;
t225 = t177 * t239;
t224 = t180 * t239;
t222 = qJ(4) * MDP(18);
t127 = -t128 + t237;
t221 = t127 * MDP(18);
t134 = t142 * t178 - t153 * t181;
t130 = t134 * t177 - t143 * t180;
t220 = t130 * MDP(29);
t131 = t134 * t180 + t143 * t177;
t219 = t131 * MDP(26);
t218 = t131 * MDP(28);
t217 = t133 * MDP(30);
t216 = t134 * MDP(20);
t215 = t134 * MDP(21);
t214 = t134 * MDP(25);
t213 = t142 * MDP(11);
t212 = t153 * MDP(12);
t158 = t181 * t176 - t178 * t228;
t211 = t158 * MDP(25);
t209 = t178 * MDP(24);
t208 = t178 * MDP(25);
t207 = t180 * MDP(26);
t205 = t239 * MDP(25);
t204 = MDP(13) - MDP(16);
t203 = MDP(14) - MDP(17);
t129 = t182 * t140 + t144 * t226 + t149 * t229;
t200 = t177 * t180 * MDP(27);
t198 = MDP(16) - t236;
t197 = -MDP(24) * t239 + MDP(21);
t196 = -t152 + t129;
t120 = t178 * t122 + t181 * t123;
t194 = MDP(28) * t180 - MDP(29) * t177;
t160 = t178 * pkin(5) - t181 * pkin(11) + qJ(4);
t150 = t180 * t160 + t178 * t225;
t151 = t177 * t160 - t178 * t224;
t193 = t150 * MDP(31) - t151 * MDP(32);
t192 = t180 * MDP(31) - t177 * MDP(32);
t124 = -t142 * pkin(4) + t196;
t189 = -MDP(20) + t194;
t188 = MDP(24) + t192;
t186 = t177 * MDP(28) + t180 * MDP(29) - pkin(11) * t244;
t118 = t143 * pkin(11) + t120;
t121 = t133 * pkin(5) - t134 * pkin(11) + t124;
t115 = -t177 * t118 + t180 * t121;
t116 = t180 * t118 + t177 * t121;
t185 = t115 * MDP(31) - t116 * MDP(32) + t217 + t218 - t220;
t184 = -MDP(22) + t186;
t171 = t181 ^ 2;
t170 = t180 ^ 2;
t169 = t178 ^ 2;
t168 = t177 ^ 2;
t155 = -qJ(2) * t230 + t164;
t147 = t180 * t158 + t177 * t229;
t146 = -t177 * t158 + t180 * t229;
t125 = t142 * pkin(3) + t190;
t1 = [(t125 ^ 2 + t127 ^ 2 + t196 ^ 2) * MDP(18) + t134 ^ 2 * MDP(19) + MDP(1) + (t167 * pkin(1) ^ 2 + t155 ^ 2 + t156 ^ 2) * MDP(7) + (t212 + 0.2e1 * t213) * t153 + (MDP(8) + MDP(23)) * t143 ^ 2 + (-0.2e1 * t130 * MDP(27) + t219) * t131 + 0.2e1 * (-t153 * MDP(10) - t142 * MDP(9) + t215) * t143 + (-0.2e1 * t143 * MDP(22) - 0.2e1 * t216 + t217 + 0.2e1 * t218 - 0.2e1 * t220) * t133 + 0.2e1 * (-t128 * t153 + t132 * t142) * MDP(13) + 0.2e1 * (-t125 * t142 - t127 * t153) * MDP(16) + 0.2e1 * (t119 * t143 + t124 * t133) * MDP(24) + 0.2e1 * (t127 * t143 - t142 * t196) * MDP(15) + (-t120 * t143 + t124 * t134) * t242 + 0.2e1 * (t129 * t153 + t132 * t143) * MDP(14) + (-t125 * t143 - t153 * t196) * t243 + (-t116 * t133 + t117 * t131) * t240 + (t115 * t133 + t117 * t130) * t241 + 0.2e1 * (t155 * t235 + t175 * t238) * MDP(4) + 0.2e1 * (-t156 * t235 - t172 * t238) * MDP(5) + 0.2e1 * (-t155 * t172 + t156 * t175) * MDP(6) * t174; (t157 * t130 + t146 * t133) * MDP(31) + (t157 * t131 - t147 * t133) * MDP(32) + (-t157 * MDP(24) - t211) * t143 + (-t175 * MDP(4) + t172 * MDP(5) - pkin(1) * MDP(7)) * t174 + (t125 * MDP(18) + t204 * t142 + t203 * t143) * t176 + ((-t143 * MDP(15) - t221) * t182 + (-t142 * MDP(15) + MDP(18) * t196 + t133 * MDP(24) + t214) * t179 + (t203 * t179 - t204 * t182) * t153) * t173; MDP(7) + (t176 ^ 2 + (t179 ^ 2 + t182 ^ 2) * t173 ^ 2) * MDP(18); t143 * MDP(10) - t213 - t212 + t128 * MDP(13) - t129 * MDP(14) + (-t143 * pkin(3) - qJ(4) * t142) * MDP(15) + (-t128 + 0.2e1 * t237) * MDP(16) + (t196 - t152) * MDP(17) + (-t127 * pkin(3) + qJ(4) * t196) * MDP(18) + qJ(4) * t214 + (qJ(4) * MDP(24) + t193) * t133 + (-t216 + t124 * MDP(24) + (-MDP(22) + t205) * t143 + t185) * t178 + (t134 * MDP(19) + t124 * MDP(25) + t131 * t207 + (-t130 * t180 - t131 * t177) * MDP(27) + (t130 * t239 + t234) * MDP(31) + (t131 * t239 + t233) * MDP(32) + t197 * t143 + t189 * t133) * t181; (t146 * t178 + t177 * t231) * MDP(31) + (-t147 * t178 + t180 * t231) * MDP(32) + ((MDP(13) - t198) * t182 + (t181 * MDP(25) - t203 + t209 + t222) * t179) * t173; t169 * MDP(30) + MDP(12) + (-0.2e1 * MDP(16) + t236) * pkin(3) + (t243 + 0.2e1 * t209 + t222) * qJ(4) + (t170 * MDP(26) + MDP(19) - 0.2e1 * t200) * t171 + (t150 * t178 + t171 * t225) * t241 + (-t151 * t178 + t171 * t224) * t240 + (qJ(4) * t242 + 0.2e1 * t178 * t189) * t181; -t153 * MDP(16) + t221 + (-t181 * t130 - t177 * t232) * MDP(31) + (-t181 * t131 - t180 * t232) * MDP(32) + (t181 * MDP(24) + MDP(15) - t208) * t143; -MDP(18) * t228; t198 + t244 * (-t169 - t171); MDP(18); t215 + t143 * MDP(23) + t119 * MDP(24) - t120 * MDP(25) + t177 * t219 + (-t177 * t130 + t131 * t180) * MDP(27) + (-pkin(5) * t130 - t233) * MDP(31) + (-pkin(5) * t131 + t234) * MDP(32) + t184 * t133; -t188 * t157 - t211; (t184 + t205) * t178 + (t177 * t207 + (-t168 + t170) * MDP(27) + (-pkin(5) * t177 - t224) * MDP(31) + (-pkin(5) * t180 + t225) * MDP(32) + t197) * t181; t188 * t181 - t208; t168 * MDP(26) + 0.2e1 * pkin(5) * t192 + MDP(23) + 0.2e1 * t200; t185; t146 * MDP(31) - t147 * MDP(32); t178 * MDP(30) + t194 * t181 + t193; -t244 * t178; t186; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
