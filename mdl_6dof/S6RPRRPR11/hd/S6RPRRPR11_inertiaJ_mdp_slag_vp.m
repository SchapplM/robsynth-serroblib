% Calculate joint inertia matrix for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR11_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:37
% EndTime: 2019-03-09 05:44:40
% DurationCPUTime: 1.50s
% Computational Cost: add. (3645->281), mult. (9642->436), div. (0->0), fcn. (11127->14), ass. (0->124)
t251 = (MDP(25) * qJ(5));
t275 = MDP(24) + t251;
t205 = sin(pkin(13));
t209 = cos(pkin(13));
t213 = sin(qJ(6));
t216 = cos(qJ(6));
t185 = t213 * t205 - t216 * t209;
t186 = t216 * t205 + t213 * t209;
t262 = pkin(11) + qJ(5);
t190 = t262 * t205;
t191 = t262 * t209;
t223 = t186 * MDP(28) - t185 * MDP(29) + (-t216 * t190 - t213 * t191) * MDP(31) - (-t213 * t190 + t216 * t191) * MDP(32);
t230 = t205 * MDP(22) + t209 * MDP(23);
t274 = -qJ(5) * t230 - MDP(18) + t223;
t206 = sin(pkin(12));
t208 = sin(pkin(6));
t210 = cos(pkin(12));
t255 = t208 * t210;
t212 = cos(pkin(6));
t266 = pkin(1) * t212;
t180 = qJ(2) * t255 + t206 * t266;
t207 = sin(pkin(7));
t211 = cos(pkin(7));
t254 = t210 * t211;
t235 = t208 * t254;
t161 = (t207 * t212 + t235) * pkin(9) + t180;
t215 = sin(qJ(3));
t218 = cos(qJ(3));
t196 = t210 * t266;
t258 = t206 * t208;
t165 = t212 * pkin(2) + t196 + (-pkin(9) * t211 - qJ(2)) * t258;
t172 = (-pkin(9) * t206 * t207 - pkin(2) * t210 - pkin(1)) * t208;
t233 = t165 * t211 + t172 * t207;
t147 = -t215 * t161 + t218 * t233;
t273 = 0.2e1 * MDP(22);
t272 = 0.2e1 * MDP(23);
t271 = 2 * MDP(24);
t270 = -2 * MDP(27);
t269 = 0.2e1 * MDP(31);
t268 = 0.2e1 * MDP(32);
t202 = t208 ^ 2;
t267 = pkin(1) * t202;
t265 = pkin(10) * t205;
t217 = cos(qJ(4));
t264 = pkin(10) * t217;
t214 = sin(qJ(4));
t263 = pkin(11) * t214;
t261 = pkin(3) * MDP(20);
t260 = pkin(4) * MDP(25);
t259 = pkin(10) * MDP(21);
t257 = t207 * t215;
t256 = t207 * t218;
t156 = -t207 * t165 + t211 * t172;
t163 = -t212 * t256 + t215 * t258 - t218 * t235;
t164 = t212 * t257 + (t206 * t218 + t215 * t254) * t208;
t143 = t163 * pkin(3) - t164 * pkin(10) + t156;
t148 = t218 * t161 + t215 * t233;
t178 = -t207 * t255 + t212 * t211;
t146 = t178 * pkin(10) + t148;
t137 = t214 * t143 + t217 * t146;
t134 = t163 * qJ(5) + t137;
t145 = -t178 * pkin(3) - t147;
t157 = t164 * t214 - t178 * t217;
t158 = t164 * t217 + t178 * t214;
t140 = t157 * pkin(4) - t158 * qJ(5) + t145;
t131 = t209 * t134 + t205 * t140;
t189 = -t217 * pkin(4) - t214 * qJ(5) - pkin(3);
t174 = t205 * t189 + t209 * t264;
t252 = MDP(21) * t214;
t136 = t217 * t143 - t214 * t146;
t135 = -t163 * pkin(4) - t136;
t250 = t135 * MDP(25);
t150 = t158 * t205 - t163 * t209;
t151 = t158 * t209 + t163 * t205;
t141 = t216 * t150 + t213 * t151;
t249 = t141 * MDP(29);
t142 = -t213 * t150 + t216 * t151;
t248 = t142 * MDP(26);
t247 = t142 * MDP(28);
t246 = t157 * MDP(30);
t245 = t158 * MDP(16);
t244 = t158 * MDP(17);
t243 = t163 * MDP(19);
t176 = t186 * t214;
t242 = t176 * MDP(29);
t177 = t185 * t214;
t241 = t177 * MDP(26);
t240 = t177 * MDP(28);
t239 = t185 * MDP(31);
t238 = t205 * MDP(23);
t237 = t209 * MDP(22);
t236 = t217 * MDP(30);
t130 = -t205 * t134 + t209 * t140;
t234 = -t130 * t205 + t131 * t209;
t182 = t214 * t211 + t217 * t257;
t166 = -t205 * t182 - t209 * t256;
t167 = t209 * t182 - t205 * t256;
t154 = t216 * t166 - t213 * t167;
t155 = t213 * t166 + t216 * t167;
t229 = t154 * MDP(31) - t155 * MDP(32);
t228 = t176 * MDP(31) - t177 * MDP(32);
t227 = -t237 + t238 - t260;
t226 = pkin(10) * MDP(25) + t230;
t225 = t150 * MDP(22) + t151 * MDP(23) + t250;
t184 = t209 * t189;
t162 = -t209 * t263 + t184 + (-pkin(5) - t265) * t217;
t168 = -t205 * t263 + t174;
t152 = t216 * t162 - t213 * t168;
t153 = t213 * t162 + t216 * t168;
t224 = t152 * MDP(31) - t153 * MDP(32) - t240 - t242;
t221 = t186 * MDP(32) + t227 + t239;
t128 = t157 * pkin(5) - t151 * pkin(11) + t130;
t129 = -t150 * pkin(11) + t131;
t126 = t216 * t128 - t213 * t129;
t127 = t213 * t128 + t216 * t129;
t220 = t126 * MDP(31) - t127 * MDP(32) + t246 + t247 - t249;
t204 = t214 ^ 2;
t200 = -t209 * pkin(5) - pkin(4);
t188 = (pkin(5) * t205 + pkin(10)) * t214;
t181 = -t217 * t211 + t214 * t257;
t179 = -qJ(2) * t258 + t196;
t173 = -t205 * t264 + t184;
t132 = t150 * pkin(5) + t135;
t1 = [MDP(1) + (t130 ^ 2 + t131 ^ 2 + t135 ^ 2) * MDP(25) + t158 ^ 2 * MDP(15) + t178 ^ 2 * MDP(12) + (t202 * pkin(1) ^ 2 + t179 ^ 2 + t180 ^ 2) * MDP(7) + (0.2e1 * t178 * MDP(10) + MDP(8) * t164) * t164 + (t141 * t270 + t248) * t142 + (-0.2e1 * t178 * MDP(11) - 0.2e1 * t164 * MDP(9) + t243 + 0.2e1 * t244) * t163 + (-0.2e1 * t163 * MDP(18) - 0.2e1 * t245 + t246 + 0.2e1 * t247 - 0.2e1 * t249) * t157 + 0.2e1 * (t147 * t178 + t156 * t163) * MDP(13) + 0.2e1 * (-t148 * t178 + t156 * t164) * MDP(14) + (-t131 * t157 + t135 * t151) * t272 + (-t127 * t157 + t132 * t142) * t268 + 0.2e1 * (t136 * t163 + t145 * t157) * MDP(20) + 0.2e1 * (-t137 * t163 + t145 * t158) * MDP(21) + (-t130 * t151 - t131 * t150) * t271 + (t130 * t157 + t135 * t150) * t273 + (t126 * t157 + t132 * t141) * t269 + 0.2e1 * (-t180 * t212 - t206 * t267) * MDP(5) + 0.2e1 * (t179 * t212 + t210 * t267) * MDP(4) + 0.2e1 * (-t179 * t206 + t180 * t210) * MDP(6) * t208; (t211 * t163 + t178 * t256) * MDP(13) + (t211 * t164 - t178 * t257) * MDP(14) + (-t157 * t256 - t181 * t163) * MDP(20) + (-t158 * t256 - t182 * t163) * MDP(21) + (t181 * t150 + t166 * t157) * MDP(22) + (t181 * t151 - t167 * t157) * MDP(23) + (-t167 * t150 - t166 * t151) * MDP(24) + (t130 * t166 + t131 * t167 + t135 * t181) * MDP(25) + (t181 * t141 + t154 * t157) * MDP(31) + (t181 * t142 - t155 * t157) * MDP(32) + (-t210 * MDP(4) + t206 * MDP(5) - pkin(1) * MDP(7)) * t208; MDP(7) + (t166 ^ 2 + t167 ^ 2 + t181 ^ 2) * MDP(25); t164 * MDP(10) - t163 * MDP(11) + t178 * MDP(12) + t147 * MDP(13) - t148 * MDP(14) - pkin(3) * t158 * MDP(21) + (-t174 * t150 - t173 * t151) * MDP(24) + (t130 * t173 + t131 * t174) * MDP(25) - t142 * t241 + (t177 * t141 - t142 * t176) * MDP(27) + (t132 * t176 + t188 * t141) * MDP(31) + (-t132 * t177 + t188 * t142) * MDP(32) + (t173 * MDP(22) - t174 * MDP(23) + t224 - t261) * t157 + (t245 - t145 * MDP(20) - t130 * MDP(22) + t131 * MDP(23) + (MDP(18) - t259) * t163 - t220) * t217 + (t158 * MDP(15) - t157 * MDP(16) + t163 * MDP(17) + t145 * MDP(21) + (-t130 * t209 - t131 * t205) * MDP(24) + t230 * t135 + (-t163 * MDP(20) + t225) * pkin(10)) * t214; (t166 * t173 + t167 * t174) * MDP(25) + t228 * t181 + (-t166 * MDP(22) + t167 * MDP(23) - t229) * t217 + ((-t166 * t209 - t167 * t205) * MDP(24) + t226 * t181) * t214 + (-t215 * MDP(14) + (MDP(20) * t217 + MDP(13) - t252) * t218) * t207; MDP(12) + t204 * MDP(15) - 0.2e1 * pkin(3) * t252 + (t204 * pkin(10) ^ 2 + t173 ^ 2 + t174 ^ 2) * MDP(25) - (t176 * t270 - t241) * t177 + (0.2e1 * t214 * MDP(16) + t236 + 0.2e1 * t240 + 0.2e1 * t242 + 0.2e1 * t261) * t217 + (-t173 * t217 + t204 * t265) * t273 + (t204 * pkin(10) * t209 + t174 * t217) * t272 + (-t152 * t217 + t188 * t176) * t269 + (t153 * t217 - t188 * t177) * t268 + (-t173 * t209 - t174 * t205) * t214 * t271; t244 + t243 + t136 * MDP(20) - t137 * MDP(21) + (-pkin(4) * t150 - t135 * t209) * MDP(22) + (-pkin(4) * t151 + t135 * t205) * MDP(23) + t234 * MDP(24) - pkin(4) * t250 + t186 * t248 + (-t186 * t141 - t142 * t185) * MDP(27) + (t132 * t185 + t200 * t141) * MDP(31) + (t132 * t186 + t200 * t142) * MDP(32) + ((-t150 * t209 + t151 * t205) * MDP(24) + t234 * MDP(25)) * qJ(5) + t274 * t157; -t182 * MDP(21) + (-MDP(20) + t221) * t181 + t275 * (-t166 * t205 + t167 * t209); -t186 * t241 + (-t186 * t176 + t177 * t185) * MDP(27) + (t200 * t176 + t188 * t185) * MDP(31) + (-t200 * t177 + t188 * t186) * MDP(32) + (-t259 - t274) * t217 + (MDP(17) - t230 * pkin(4) + (-MDP(20) + t227) * pkin(10)) * t214 + t275 * (-t173 * t205 + t174 * t209); 0.2e1 * t200 * t239 + MDP(19) + (0.2e1 * t237 - 0.2e1 * t238 + t260) * pkin(4) + (MDP(26) * t186 + t185 * t270 + t200 * t268) * t186 + (t271 + t251) * (t205 ^ 2 + t209 ^ 2) * qJ(5); t141 * MDP(31) + t142 * MDP(32) + t225; t181 * MDP(25); t214 * t226 + t228; t221; MDP(25); t220; t229; t224 - t236; t223; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
