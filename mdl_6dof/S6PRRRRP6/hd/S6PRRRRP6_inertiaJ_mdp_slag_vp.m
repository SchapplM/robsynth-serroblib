% Calculate joint inertia matrix for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:29
% EndTime: 2019-03-09 00:34:33
% DurationCPUTime: 1.32s
% Computational Cost: add. (1453->271), mult. (3471->401), div. (0->0), fcn. (3855->12), ass. (0->107)
t180 = sin(pkin(7));
t255 = 0.2e1 * t180;
t254 = MDP(29) * pkin(11) + MDP(27);
t253 = MDP(24) + MDP(26);
t216 = MDP(25) - MDP(28);
t252 = 2 * MDP(26);
t251 = 2 * MDP(27);
t186 = sin(qJ(3));
t250 = pkin(2) * t186;
t190 = cos(qJ(3));
t249 = pkin(2) * t190;
t182 = cos(pkin(7));
t185 = sin(qJ(4));
t189 = cos(qJ(4));
t237 = t180 * t186;
t160 = -t189 * t182 + t185 * t237;
t248 = pkin(5) * t160;
t184 = sin(qJ(5));
t247 = pkin(10) * t184;
t188 = cos(qJ(5));
t246 = pkin(10) * t188;
t244 = pkin(3) * MDP(17);
t243 = qJ(6) * t160;
t181 = sin(pkin(6));
t183 = cos(pkin(6));
t187 = sin(qJ(2));
t191 = cos(qJ(2));
t234 = t182 * t191;
t146 = t183 * t237 + (t186 * t234 + t187 * t190) * t181;
t159 = -t180 * t181 * t191 + t182 * t183;
t137 = t146 * t185 - t159 * t189;
t242 = t137 * t185;
t236 = t180 * t190;
t215 = pkin(9) * t236;
t156 = t215 + (pkin(10) + t250) * t182;
t157 = (-pkin(3) * t190 - pkin(10) * t186 - pkin(2)) * t180;
t142 = -t185 * t156 + t157 * t189;
t140 = pkin(4) * t236 - t142;
t241 = t140 * t184;
t240 = t140 * t188;
t161 = t182 * t185 + t189 * t237;
t147 = t161 * t184 + t188 * t236;
t239 = t147 * t188;
t148 = t161 * t188 - t184 * t236;
t238 = t148 * t184;
t235 = t182 * MDP(9);
t233 = t184 * t189;
t232 = t186 * MDP(7);
t171 = pkin(9) * t237;
t155 = t171 + (-pkin(3) - t249) * t182;
t139 = pkin(4) * t160 - pkin(11) * t161 + t155;
t143 = t156 * t189 + t157 * t185;
t141 = -pkin(11) * t236 + t143;
t128 = t184 * t139 + t188 * t141;
t168 = -pkin(4) * t189 - pkin(11) * t185 - pkin(3);
t153 = t184 * t168 + t189 * t246;
t177 = t184 ^ 2;
t179 = t188 ^ 2;
t231 = t177 + t179;
t230 = MDP(16) * t190;
t229 = MDP(27) * t185;
t228 = t147 * MDP(22);
t227 = t148 * MDP(19);
t226 = t148 * MDP(21);
t225 = t160 * MDP(23);
t224 = t161 * MDP(13);
t223 = t161 * MDP(14);
t206 = -pkin(5) * t188 - qJ(6) * t184;
t167 = -pkin(4) + t206;
t222 = t167 * MDP(29);
t221 = t184 * MDP(21);
t220 = t185 * MDP(18);
t219 = t188 * MDP(19);
t218 = t188 * MDP(20);
t214 = pkin(10) * MDP(17) - MDP(14);
t213 = pkin(10) * MDP(18) - MDP(15);
t212 = -MDP(29) * pkin(5) - MDP(26);
t209 = -t188 * t139 + t141 * t184;
t208 = -0.2e1 * qJ(6) * MDP(28) - MDP(23);
t207 = -MDP(24) + t212;
t205 = -pkin(5) * t184 + qJ(6) * t188;
t204 = MDP(29) * qJ(6) - t216;
t125 = t128 + t243;
t126 = t209 - t248;
t203 = t125 * t188 + t126 * t184;
t149 = -qJ(6) * t189 + t153;
t165 = t188 * t168;
t150 = -t165 + (pkin(5) + t247) * t189;
t200 = t149 * t188 + t150 * t184;
t199 = t188 * MDP(21) - t184 * MDP(22);
t198 = t188 * MDP(22) + t221;
t197 = (-pkin(10) * t233 + t165) * MDP(24) - t153 * MDP(25);
t196 = t184 * MDP(26) - t188 * MDP(28);
t195 = -MDP(13) + t199;
t194 = t150 * MDP(26) - t149 * MDP(28) - t197;
t193 = -MDP(24) * t209 - t128 * MDP(25) + t225 + t226 - t228;
t176 = t180 ^ 2;
t173 = pkin(11) * t233;
t163 = t182 * t250 + t215;
t162 = t182 * t249 - t171;
t158 = (pkin(10) - t205) * t185;
t145 = -t183 * t236 + (t186 * t187 - t190 * t234) * t181;
t138 = t146 * t189 + t159 * t185;
t131 = t138 * t188 + t145 * t184;
t130 = t138 * t184 - t145 * t188;
t129 = pkin(5) * t147 - qJ(6) * t148 + t140;
t1 = [MDP(1) + (t130 ^ 2 + t131 ^ 2 + t137 ^ 2) * MDP(29); (-t145 * t182 - t159 * t236) * MDP(10) + (-t146 * t182 + t159 * t237) * MDP(11) + (t137 * t236 + t145 * t160) * MDP(17) + (t138 * t236 + t145 * t161) * MDP(18) + (t130 * t148 - t131 * t147) * MDP(27) + (t125 * t131 + t126 * t130 + t129 * t137) * MDP(29) + (MDP(3) * t191 - MDP(4) * t187) * t181 + t253 * (-t130 * t160 + t137 * t147) - t216 * (t131 * t160 - t137 * t148); t161 ^ 2 * MDP(12) + (t125 ^ 2 + t126 ^ 2 + t129 ^ 2) * MDP(29) + t176 * t186 ^ 2 * MDP(5) + MDP(2) + (t232 * t255 + t235) * t182 + (-0.2e1 * t147 * MDP(20) + t227) * t148 + ((MDP(8) * t182 - t223) * t255 + (0.2e1 * MDP(6) * t186 + t230) * t176) * t190 + (0.2e1 * MDP(15) * t236 - 0.2e1 * t224 + t225 + 0.2e1 * t226 - 0.2e1 * t228) * t160 + 0.2e1 * (t162 * t182 + t176 * t249) * MDP(10) + 0.2e1 * (-t142 * t236 + t155 * t160) * MDP(17) + 0.2e1 * (t143 * t236 + t155 * t161) * MDP(18) + 0.2e1 * (-t163 * t182 - t176 * t250) * MDP(11) + 0.2e1 * (t140 * t147 - t160 * t209) * MDP(24) + 0.2e1 * (t125 * t160 - t129 * t148) * MDP(28) + (-t126 * t160 + t129 * t147) * t252 + 0.2e1 * (-t128 * t160 + t140 * t148) * MDP(25) + (-t125 * t147 + t126 * t148) * t251; -t146 * MDP(11) + (t130 * t150 + t131 * t149 + t137 * t158) * MDP(29) + (t130 * t188 - t131 * t184) * t229 + t216 * (t131 * t189 + t188 * t242) + (-t189 * MDP(17) - MDP(10) + t220) * t145 + t253 * (t130 * t189 + t184 * t242); t235 + t162 * MDP(10) - t163 * MDP(11) - pkin(3) * t161 * MDP(18) + (-t147 * t149 + t148 * t150) * MDP(27) + (t125 * t149 + t126 * t150) * MDP(29) + (t190 * MDP(8) + t232) * t180 + (t147 * MDP(26) - t148 * MDP(28) + t129 * MDP(29)) * t158 + (-t194 - t244) * t160 + (-t155 * MDP(17) + t126 * MDP(26) - t125 * MDP(28) + t213 * t236 - t193 + t224) * t189 + (t161 * MDP(12) + t155 * MDP(18) + t148 * t219 + (-t238 - t239) * MDP(20) + (pkin(10) * t147 + t241) * MDP(24) + (pkin(10) * t148 + t240) * MDP(25) + (-t125 * t184 + t126 * t188) * MDP(27) + t214 * t236 + t196 * t129 + t195 * t160) * t185; MDP(9) - 0.2e1 * pkin(3) * t220 + (t149 ^ 2 + t150 ^ 2 + t158 ^ 2) * MDP(29) + 0.2e1 * ((-t149 * t184 + t150 * t188) * MDP(27) + t196 * t158) * t185 + (t189 * MDP(23) - 0.2e1 * t195 * t185 + 0.2e1 * t194 + 0.2e1 * t244) * t189 + (t179 * MDP(19) - 0.2e1 * t184 * t218 + MDP(12) + 0.2e1 * (t184 * MDP(24) + t188 * MDP(25)) * pkin(10)) * t185 ^ 2; -t138 * MDP(18) + (t216 * t184 - t188 * t253 - MDP(17) + t222) * t137 + t254 * (t130 * t184 + t131 * t188); t223 - t180 * t230 + t142 * MDP(17) - t143 * MDP(18) + t184 * t227 + (-t147 * t184 + t148 * t188) * MDP(20) + (-pkin(4) * t147 - t240) * MDP(24) + (-pkin(4) * t148 + t241) * MDP(25) + (-t129 * t188 + t147 * t167) * MDP(26) + t203 * MDP(27) + (-t129 * t184 - t148 * t167) * MDP(28) + t129 * t222 + ((t238 - t239) * MDP(27) + t203 * MDP(29)) * pkin(11) + (-MDP(15) + (-t253 * t184 - t216 * t188) * pkin(11) + t198) * t160; t173 * MDP(24) + (-t158 * t188 + t173) * MDP(26) + t200 * MDP(27) - t158 * t184 * MDP(28) + (t200 * pkin(11) + t158 * t167) * MDP(29) + (-t221 + (t216 * pkin(11) - MDP(22)) * t188 - t213) * t189 + (t184 * t219 + (-t177 + t179) * MDP(20) + (-pkin(4) * t184 - t246) * MDP(24) + (-pkin(4) * t188 + t247) * MDP(25) + t196 * t167 - t214) * t185; MDP(16) + t177 * MDP(19) + (t231 * pkin(11) ^ 2 + t167 ^ 2) * MDP(29) + t231 * pkin(11) * t251 + 0.2e1 * (MDP(24) * pkin(4) - MDP(26) * t167) * t188 + 0.2e1 * (-MDP(25) * pkin(4) - MDP(28) * t167 + t218) * t184; t207 * t130 + t204 * t131; (-t209 + 0.2e1 * t248) * MDP(26) + (-pkin(5) * t148 - qJ(6) * t147) * MDP(27) + (t128 + 0.2e1 * t243) * MDP(28) + (-pkin(5) * t126 + qJ(6) * t125) * MDP(29) + t193; t165 * MDP(26) + t153 * MDP(28) + (-pkin(5) * t150 + qJ(6) * t149) * MDP(29) + ((-0.2e1 * pkin(5) - t247) * MDP(26) + t208) * t189 + (t206 * MDP(27) + t199) * t185 + t197; t205 * MDP(27) + (t207 * t184 + t204 * t188) * pkin(11) + t198; pkin(5) * t252 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29) - t208; t130 * MDP(29); -t160 * MDP(26) + t148 * MDP(27) + t126 * MDP(29); t189 * MDP(26) + MDP(29) * t150 + t188 * t229; t254 * t184; t212; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
