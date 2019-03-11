% Calculate joint inertia matrix for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR7_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:30
% EndTime: 2019-03-09 18:39:35
% DurationCPUTime: 1.33s
% Computational Cost: add. (2575->243), mult. (5717->359), div. (0->0), fcn. (6697->12), ass. (0->120)
t213 = sin(qJ(6));
t217 = cos(qJ(6));
t256 = t213 * MDP(29) + t217 * MDP(30);
t215 = sin(qJ(3));
t218 = cos(qJ(3));
t281 = t215 * MDP(13) + t218 * MDP(14) + (-MDP(16) * t215 - MDP(17) * t218) * pkin(9);
t212 = cos(pkin(6));
t210 = sin(pkin(6));
t216 = sin(qJ(2));
t259 = t210 * t216;
t187 = -t212 * t218 + t215 * t259;
t188 = t212 * t215 + t218 * t259;
t209 = sin(pkin(12));
t211 = cos(pkin(12));
t166 = -t187 * t209 + t188 * t211;
t214 = sin(qJ(5));
t233 = t187 * t211 + t188 * t209;
t272 = cos(qJ(5));
t151 = t166 * t214 + t272 * t233;
t152 = t272 * t166 - t214 * t233;
t280 = t152 * MDP(22) - t151 * MDP(23);
t228 = MDP(32) * t217 - MDP(33) * t213;
t267 = -qJ(4) - pkin(9);
t195 = t267 * t215;
t196 = t267 * t218;
t176 = t209 * t195 - t211 * t196;
t232 = t209 * t215 - t211 * t218;
t164 = -t232 * pkin(10) + t176;
t175 = t211 * t195 + t196 * t209;
t191 = t209 * t218 + t211 * t215;
t226 = -pkin(10) * t191 + t175;
t147 = t164 * t214 - t272 * t226;
t148 = t272 * t164 + t214 * t226;
t173 = t191 * t214 + t272 * t232;
t174 = t272 * t191 - t214 * t232;
t278 = -t174 * MDP(22) + t173 * MDP(23) + t147 * MDP(25) + t148 * MDP(26);
t277 = 0.2e1 * MDP(16);
t276 = 2 * MDP(18);
t275 = 0.2e1 * MDP(25);
t274 = 0.2e1 * MDP(32);
t273 = 0.2e1 * MDP(33);
t271 = pkin(1) * t216;
t219 = cos(qJ(2));
t270 = pkin(1) * t219;
t269 = pkin(3) * t209;
t266 = MDP(18) * pkin(3);
t265 = MDP(19) * pkin(3);
t258 = t210 * t219;
t239 = pkin(8) * t258;
t180 = t239 + (pkin(9) + t271) * t212;
t181 = (-pkin(2) * t219 - pkin(9) * t216 - pkin(1)) * t210;
t160 = -t180 * t215 + t218 * t181;
t157 = -pkin(3) * t258 - qJ(4) * t188 + t160;
t161 = t180 * t218 + t181 * t215;
t159 = -qJ(4) * t187 + t161;
t136 = t211 * t157 - t159 * t209;
t132 = -pkin(4) * t258 - pkin(10) * t166 + t136;
t137 = t209 * t157 + t211 * t159;
t135 = -t233 * pkin(10) + t137;
t236 = -t272 * t132 + t135 * t214;
t125 = pkin(5) * t258 + t236;
t264 = t125 * t217;
t263 = t147 * t217;
t262 = t151 * t213;
t261 = t151 * t217;
t260 = t174 * t213;
t257 = t212 * t216;
t139 = t152 * t213 + t217 * t258;
t251 = t139 * MDP(30);
t140 = t152 * t217 - t213 * t258;
t250 = t140 * MDP(27);
t249 = t140 * MDP(29);
t248 = t151 * MDP(31);
t247 = t152 * MDP(21);
t246 = t152 * MDP(26);
t245 = t173 * MDP(31);
t244 = t174 * MDP(26);
t200 = pkin(3) * t211 + pkin(4);
t184 = t272 * t200 - t214 * t269;
t243 = t184 * MDP(25);
t185 = t214 * t200 + t272 * t269;
t242 = t185 * MDP(26);
t241 = t215 * MDP(11);
t240 = t217 * MDP(27);
t201 = -t218 * pkin(3) - pkin(2);
t238 = t213 * t217 * MDP(28);
t206 = t213 ^ 2;
t237 = t206 * MDP(27) + MDP(24) + 0.2e1 * t238;
t235 = -pkin(5) * t174 - pkin(11) * t173;
t182 = -pkin(5) - t184;
t183 = pkin(11) + t185;
t234 = -t173 * t183 + t174 * t182;
t231 = t188 * MDP(13) - t187 * MDP(14);
t229 = t217 * MDP(29) - t213 * MDP(30);
t227 = t213 * MDP(32) + t217 * MDP(33);
t198 = pkin(8) * t259;
t179 = t198 + (-pkin(2) - t270) * t212;
t128 = t214 * t132 + t272 * t135;
t225 = -MDP(21) + t229;
t224 = MDP(25) + t228;
t207 = t217 ^ 2;
t223 = (-t206 + t207) * t174 * MDP(28) + t240 * t260 - t278 + t256 * t173;
t169 = t187 * pkin(3) + t179;
t222 = -MDP(24) * t258 + (-t139 * t213 + t140 * t217) * MDP(28) + t213 * t250 + t280 + t256 * t151;
t178 = t232 * pkin(4) + t201;
t126 = -pkin(11) * t258 + t128;
t153 = t233 * pkin(4) + t169;
t129 = t151 * pkin(5) - t152 * pkin(11) + t153;
t122 = -t126 * t213 + t129 * t217;
t123 = t126 * t217 + t129 * t213;
t221 = t122 * MDP(32) - t123 * MDP(33) + t248 + t249 - t251;
t205 = t210 ^ 2;
t190 = pkin(1) * t257 + t239;
t189 = t212 * t270 - t198;
t146 = t173 * pkin(5) - t174 * pkin(11) + t178;
t141 = t147 * t213;
t134 = t146 * t213 + t148 * t217;
t133 = t146 * t217 - t148 * t213;
t124 = t125 * t213;
t1 = [t212 ^ 2 * MDP(8) + (t136 ^ 2 + t137 ^ 2 + t169 ^ 2) * MDP(19) + t152 ^ 2 * MDP(20) + MDP(1) + (t188 * MDP(11) - 0.2e1 * t187 * MDP(12)) * t188 + (-0.2e1 * t139 * MDP(28) + t250) * t140 + ((MDP(4) * t216 + 0.2e1 * MDP(5) * t219) * t216 + (MDP(15) + MDP(24)) * t219 ^ 2) * t205 + (-0.2e1 * t247 + t248 + 0.2e1 * t249 - 0.2e1 * t251) * t151 + 0.2e1 * (MDP(6) * t257 + (MDP(7) * t212 - t231 - t280) * t219) * t210 + (t122 * t151 + t125 * t139) * t274 + (-t123 * t151 + t125 * t140) * t273 + (-t160 * t258 + t179 * t187) * t277 + 0.2e1 * (t161 * t258 + t179 * t188) * MDP(17) + (t151 * t153 + t236 * t258) * t275 + 0.2e1 * (t128 * t258 + t152 * t153) * MDP(26) + 0.2e1 * (t189 * t212 + t205 * t270) * MDP(9) + (-t136 * t166 - t137 * t233) * t276 + 0.2e1 * (-t190 * t212 - t205 * t271) * MDP(10); (-t134 * t151 + t140 * t147) * MDP(33) + (-t136 * t191 - t137 * t232 - t175 * t166 - t176 * t233) * MDP(18) + (-t187 * t215 + t188 * t218) * MDP(12) + (-pkin(2) * t188 + t179 * t215) * MDP(17) + (-pkin(2) * t187 - t179 * t218) * MDP(16) + (t133 * t151 + t139 * t147) * MDP(32) + t212 * MDP(8) + (t136 * t175 + t137 * t176 + t169 * t201) * MDP(19) + t189 * MDP(9) - t190 * MDP(10) + t188 * t241 + (t151 * MDP(25) + t246) * t178 + (t153 * MDP(25) + t221 - t247) * t173 + ((-t139 * t217 - t140 * t213) * MDP(28) + t153 * MDP(26) + t152 * MDP(20) + t140 * t240 + t227 * t125 + t225 * t151) * t174 + (MDP(6) * t216 + (MDP(7) + t278 - t281) * t219) * t210; MDP(8) + pkin(2) * t218 * t277 + (t175 ^ 2 + t176 ^ 2 + t201 ^ 2) * MDP(19) + 0.2e1 * t178 * t244 + (0.2e1 * t218 * MDP(12) - 0.2e1 * pkin(2) * MDP(17) + t241) * t215 + (MDP(27) * t207 + MDP(20) - 0.2e1 * t238) * t174 ^ 2 + (0.2e1 * t225 * t174 + t178 * t275 + t245) * t173 + (-t175 * t191 - t176 * t232) * t276 + (t133 * t173 + t147 * t260) * t274 + (-t134 * t173 + t174 * t263) * t273; (t139 * t182 - t183 * t262 - t264) * MDP(32) + (t140 * t182 - t183 * t261 + t124) * MDP(33) + (-t184 * t258 - t236) * MDP(25) + (t185 * t258 - t128) * MDP(26) + (t136 * t211 + t137 * t209) * t265 - MDP(15) * t258 + (-t211 * t166 - t209 * t233) * t266 + t160 * MDP(16) - t161 * MDP(17) + t222 + t231; (-t211 * t191 - t209 * t232) * t266 + (t234 * t217 + t141) * MDP(33) + (t175 * t211 + t176 * t209) * t265 + (t234 * t213 - t263) * MDP(32) + t223 + t281; MDP(15) + (t209 ^ 2 + t211 ^ 2) * MDP(19) * pkin(3) ^ 2 - 0.2e1 * t228 * t182 + 0.2e1 * t243 - 0.2e1 * t242 + t237; MDP(19) * t169 + t224 * t151 + t246; MDP(19) * t201 + t224 * t173 + t244; 0; MDP(19); -t236 * MDP(25) - t128 * MDP(26) + (-pkin(5) * t139 - pkin(11) * t262 - t264) * MDP(32) + (-pkin(5) * t140 - pkin(11) * t261 + t124) * MDP(33) + t222; (t235 * t213 - t263) * MDP(32) + (t235 * t217 + t141) * MDP(33) + t223; t237 - t242 + t243 + t228 * (pkin(5) - t182); 0; 0.2e1 * pkin(5) * t228 + t237; t221; t133 * MDP(32) - t134 * MDP(33) + t229 * t174 + t245; -t227 * t183 + t256; t228; -t227 * pkin(11) + t256; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
