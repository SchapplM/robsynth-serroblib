% Calculate joint inertia matrix for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:31
% EndTime: 2019-03-10 01:31:36
% DurationCPUTime: 1.59s
% Computational Cost: add. (2014->258), mult. (3831->342), div. (0->0), fcn. (4041->8), ass. (0->102)
t184 = sin(qJ(3));
t188 = cos(qJ(3));
t229 = pkin(8) + pkin(9);
t165 = t229 * t184;
t166 = t229 * t188;
t183 = sin(qJ(4));
t187 = cos(qJ(4));
t143 = -t165 * t187 - t166 * t183;
t144 = -t165 * t183 + t166 * t187;
t161 = t183 * t184 - t187 * t188;
t162 = t183 * t188 + t184 * t187;
t132 = -pkin(10) * t161 + t144;
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t194 = -pkin(10) * t162 + t143;
t122 = t182 * t132 - t186 * t194;
t123 = t186 * t132 + t182 * t194;
t136 = t186 * t161 + t162 * t182;
t137 = -t161 * t182 + t162 * t186;
t208 = MDP(30) + MDP(32);
t241 = MDP(31) - MDP(34);
t196 = t137 * MDP(27) - t136 * MDP(28) - t208 * t122 - t241 * t123;
t193 = t162 * MDP(20) - t161 * MDP(21) + t143 * MDP(23) - t144 * MDP(24) + t196;
t249 = t193 - (MDP(16) * t184 + MDP(17) * t188) * pkin(8) + t184 * MDP(13) + t188 * MDP(14);
t211 = t182 * MDP(31);
t246 = (MDP(30) * t186 - t211) * pkin(4);
t210 = t187 * MDP(23);
t245 = (-t183 * MDP(24) + t210) * pkin(3);
t185 = sin(qJ(2));
t152 = t162 * t185;
t153 = t161 * t185;
t129 = t186 * t152 - t153 * t182;
t130 = -t152 * t182 - t153 * t186;
t219 = t130 * MDP(27) - t129 * MDP(28);
t243 = t153 * MDP(20) + t152 * MDP(21);
t189 = cos(qJ(2));
t164 = -pkin(2) * t189 - pkin(8) * t185 - pkin(1);
t160 = t188 * t164;
t222 = pkin(9) * t185;
t225 = pkin(7) * t184;
t138 = -t188 * t222 + t160 + (-pkin(3) - t225) * t189;
t223 = pkin(7) * t189;
t207 = t188 * t223;
t140 = t207 + (t164 - t222) * t184;
t125 = t187 * t138 - t140 * t183;
t115 = -pkin(4) * t189 + pkin(10) * t153 + t125;
t126 = t138 * t183 + t187 * t140;
t124 = -pkin(10) * t152 + t126;
t110 = t186 * t115 - t182 * t124;
t242 = t110 * MDP(30) + t219;
t111 = t182 * t115 + t186 * t124;
t240 = t125 * MDP(23) - t126 * MDP(24) + t110 * MDP(32) - t241 * t111 + t242 - t243;
t237 = -2 * MDP(19);
t236 = 0.2e1 * MDP(23);
t235 = 0.2e1 * MDP(24);
t234 = -2 * MDP(26);
t233 = 0.2e1 * MDP(31);
t232 = 0.2e1 * MDP(32);
t231 = 2 * MDP(33);
t230 = 2 * MDP(34);
t228 = pkin(3) * t183;
t227 = pkin(4) * t186;
t226 = pkin(5) * t189;
t224 = pkin(7) * t188;
t221 = qJ(6) * t189;
t220 = t184 * t188;
t174 = pkin(3) * t187 + pkin(4);
t156 = t182 * t174 + t186 * t228;
t163 = (pkin(3) * t184 + pkin(7)) * t185;
t218 = MDP(18) * t162;
t217 = MDP(25) * t137;
t168 = t186 * t174;
t202 = t182 * t228 - t168;
t213 = t202 * MDP(30);
t212 = t156 * MDP(31);
t209 = MDP(22) + MDP(29);
t190 = 0.2e1 * qJ(6);
t206 = t190 + t156;
t175 = -pkin(3) * t188 - pkin(2);
t205 = MDP(12) * t220;
t204 = t168 + t227;
t203 = MDP(15) + t209;
t139 = pkin(4) * t152 + t163;
t148 = pkin(4) * t161 + t175;
t201 = MDP(13) * t188 - MDP(14) * t184;
t195 = MDP(29) - t212 - t213;
t191 = 0.2e1 * pkin(5);
t180 = t188 ^ 2;
t179 = t185 ^ 2;
t178 = t184 ^ 2;
t176 = t182 * pkin(4);
t172 = pkin(5) + t227;
t171 = t176 + qJ(6);
t154 = -pkin(5) + t202;
t151 = qJ(6) + t156;
t147 = t164 * t184 + t207;
t146 = -t184 * t223 + t160;
t121 = pkin(5) * t136 - qJ(6) * t137 + t148;
t112 = pkin(5) * t129 - qJ(6) * t130 + t139;
t109 = -t110 + t226;
t108 = t111 - t221;
t1 = [(t108 ^ 2 + t109 ^ 2 + t112 ^ 2) * MDP(35) - 0.2e1 * pkin(1) * t185 * MDP(10) + MDP(1) - (-MDP(18) * t153 + t152 * t237) * t153 + (MDP(25) * t130 + t129 * t234) * t130 + t203 * t189 ^ 2 + (MDP(11) * t180 + MDP(4) - 0.2e1 * t205) * t179 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t201) * t185 + t243 - t219) * t189 + 0.2e1 * (t147 * t189 + t179 * t224) * MDP(17) + 0.2e1 * (-t146 * t189 + t179 * t225) * MDP(16) + (-t125 * t189 + t152 * t163) * t236 + (t126 * t189 - t153 * t163) * t235 + (-t108 * t189 - t112 * t130) * t230 + 0.2e1 * (-t110 * t189 + t129 * t139) * MDP(30) + (t111 * t189 + t130 * t139) * t233 + (t109 * t189 + t112 * t129) * t232 + (-t108 * t129 + t109 * t130) * t231; -t153 * t218 + t130 * t217 + (t108 * t123 + t109 * t122 + t112 * t121) * MDP(35) + (-t129 * t137 - t130 * t136) * MDP(26) + (-t108 * t136 + t109 * t137 + t122 * t130 - t123 * t129) * MDP(33) + (-t152 * t162 + t153 * t161) * MDP(19) + (t112 * t136 + t121 * t129) * MDP(32) + (-t112 * t137 - t121 * t130) * MDP(34) + (t129 * t148 + t136 * t139) * MDP(30) + (t130 * t148 + t137 * t139) * MDP(31) + (t152 * t175 + t161 * t163) * MDP(23) + (-t153 * t175 + t162 * t163) * MDP(24) + (-pkin(7) * MDP(9) + MDP(11) * t220 + MDP(6) + (-t178 + t180) * MDP(12) + (-pkin(2) * t184 - t224) * MDP(16) + (-pkin(2) * t188 + t225) * MDP(17)) * t185 + (-pkin(7) * MDP(10) + MDP(7) - t249) * t189; MDP(8) + t178 * MDP(11) + 0.2e1 * t205 + t175 * t161 * t236 + (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) * MDP(35) + (-0.2e1 * MDP(34) * t121 + t122 * t231 + t136 * t234 + t148 * t233 + t217) * t137 + (t161 * t237 + t175 * t235 + t218) * t162 + 0.2e1 * (MDP(16) * t188 - MDP(17) * t184) * pkin(2) + 0.2e1 * (MDP(30) * t148 + MDP(32) * t121 - MDP(33) * t123) * t136; (-MDP(15) - MDP(22) + (-pkin(5) + t154) * MDP(32) + (-qJ(6) - t151) * MDP(34) - t245 - t195) * t189 + t201 * t185 + t146 * MDP(16) - t147 * MDP(17) + (-t129 * t151 + t130 * t154) * MDP(33) + (t108 * t151 + t109 * t154) * MDP(35) + t240; (-t136 * t151 + t137 * t154) * MDP(33) + (t122 * t154 + t123 * t151) * MDP(35) + t249; (t151 ^ 2 + t154 ^ 2) * MDP(35) + 0.2e1 * t245 - 0.2e1 * t213 - 0.2e1 * t212 - 0.2e1 * MDP(32) * t154 + t151 * t230 + t203; (-t129 * t171 - t130 * t172) * MDP(33) + (t108 * t171 - t109 * t172) * MDP(35) + ((-pkin(5) - t172) * MDP(32) + (-qJ(6) - t171) * MDP(34) - t246 - t209) * t189 + t240; (-t136 * t171 - t137 * t172) * MDP(33) + (-t122 * t172 + t123 * t171) * MDP(35) + t193; t204 * MDP(30) + (-t176 - t156) * MDP(31) + (t191 + t204) * MDP(32) + (t176 + t206) * MDP(34) + (t151 * t171 - t154 * t172) * MDP(35) + (t210 + (-t208 * t182 - MDP(24)) * t183) * pkin(3) + t209; (t171 ^ 2 + t172 ^ 2) * MDP(35) + 0.2e1 * t246 + t172 * t232 + t171 * t230 + t209; -t189 * MDP(29) - t111 * MDP(31) + (t110 - 0.2e1 * t226) * MDP(32) + (-pkin(5) * t130 - qJ(6) * t129) * MDP(33) + (t111 - 0.2e1 * t221) * MDP(34) + (-pkin(5) * t109 + qJ(6) * t108) * MDP(35) + t242; (-pkin(5) * t137 - qJ(6) * t136) * MDP(33) + (-pkin(5) * t122 + qJ(6) * t123) * MDP(35) + t196; (t191 - t202) * MDP(32) + t206 * MDP(34) + (-pkin(5) * t154 + qJ(6) * t151) * MDP(35) + t195; MDP(29) + t191 * MDP(32) + (t190 + t176) * MDP(34) + (pkin(5) * t172 + qJ(6) * t171) * MDP(35) + (t208 * t186 - t211) * pkin(4); MDP(29) + pkin(5) * t232 + qJ(6) * t230 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35); MDP(32) * t189 + MDP(33) * t130 + MDP(35) * t109; MDP(33) * t137 + MDP(35) * t122; MDP(35) * t154 - MDP(32); -MDP(35) * t172 - MDP(32); -MDP(35) * pkin(5) - MDP(32); MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
