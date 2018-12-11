% Calculate Gravitation load on the joints for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR14_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x6] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:32:39
% EndTime: 2018-12-10 18:32:49
% DurationCPUTime: 2.24s
% Computational Cost: add. (4781->194), mult. (4878->310), div. (0->0), fcn. (4867->30), ass. (0->107)
t237 = sin(qJ(2));
t238 = sin(qJ(1));
t243 = cos(qJ(1));
t273 = pkin(6) + qJ(2);
t264 = cos(t273) / 0.2e1;
t274 = pkin(6) - qJ(2);
t269 = cos(t274);
t253 = t269 / 0.2e1 + t264;
t201 = t238 * t237 - t243 * t253;
t262 = sin(t273) / 0.2e1;
t267 = sin(t274);
t213 = t262 - t267 / 0.2e1;
t242 = cos(qJ(2));
t202 = t243 * t213 + t238 * t242;
t225 = pkin(7) - pkin(14);
t217 = cos(t225) / 0.2e1;
t224 = pkin(7) + pkin(14);
t223 = cos(t224);
t209 = t217 + t223 / 0.2e1;
t226 = sin(pkin(14));
t229 = sin(pkin(6));
t216 = sin(t224) / 0.2e1;
t222 = sin(t225);
t270 = t216 + t222 / 0.2e1;
t265 = t229 * t270;
t185 = t201 * t209 + t202 * t226 + t243 * t265;
t208 = t216 - t222 / 0.2e1;
t210 = t217 - t223 / 0.2e1;
t230 = cos(pkin(14));
t277 = t229 * t243;
t186 = t201 * t208 - t202 * t230 + t210 * t277;
t228 = sin(pkin(7));
t232 = cos(pkin(7));
t199 = -t201 * t228 + t232 * t277;
t271 = pkin(8) + qJ(4);
t261 = sin(t271) / 0.2e1;
t272 = pkin(8) - qJ(4);
t266 = sin(t272);
t211 = t261 - t266 / 0.2e1;
t263 = cos(t271) / 0.2e1;
t268 = cos(t272);
t214 = t263 - t268 / 0.2e1;
t241 = cos(qJ(4));
t162 = t185 * t211 + t186 * t241 - t199 * t214;
t227 = sin(pkin(8));
t231 = cos(pkin(8));
t177 = -t185 * t227 + t199 * t231;
t235 = sin(qJ(5));
t240 = cos(qJ(5));
t151 = t162 * t240 + t177 * t235;
t236 = sin(qJ(4));
t251 = t261 + t266 / 0.2e1;
t252 = t268 / 0.2e1 + t263;
t161 = -t185 * t252 + t186 * t236 - t199 * t251;
t234 = sin(qJ(6));
t239 = cos(qJ(6));
t293 = t151 * t234 - t161 * t239;
t292 = t151 * t239 + t161 * t234;
t289 = t162 * t235 - t177 * t240;
t204 = -t243 * t237 - t238 * t253;
t278 = t229 * t238;
t259 = -t204 * t228 + t232 * t278;
t286 = -g(1) * t199 - g(2) * t259;
t280 = t228 * t214;
t279 = t228 * t231;
t276 = t234 * t240;
t275 = t239 * t240;
t212 = t262 + t267 / 0.2e1;
t233 = cos(pkin(6));
t260 = t212 * t228 - t233 * t232;
t206 = t238 * t213 - t243 * t242;
t215 = t264 - t269 / 0.2e1;
t250 = t212 * t209 + t215 * t226 + t233 * t270;
t248 = t228 * t251;
t247 = -t204 * t209 - t206 * t226 - t238 * t265;
t194 = t212 * t208 + t233 * t210 - t215 * t230;
t246 = t194 * t241 + t250 * t211 + t260 * t214;
t187 = t204 * t208 - t206 * t230 + t210 * t278;
t244 = t187 * t241 - t247 * t211 - t259 * t214;
t198 = t215 * t208 + t212 * t230;
t197 = t215 * t209 - t212 * t226;
t193 = t204 * t230 + t206 * t208;
t192 = -t204 * t226 + t206 * t209;
t191 = -t201 * t230 - t202 * t208;
t190 = t201 * t226 - t202 * t209;
t189 = -t197 * t227 - t215 * t279;
t183 = -t250 * t227 - t260 * t231;
t180 = -t192 * t227 - t206 * t279;
t179 = -t190 * t227 + t202 * t279;
t178 = t247 * t227 + t259 * t231;
t175 = t197 * t211 + t198 * t241 + t215 * t280;
t174 = -t197 * t252 + t198 * t236 + t215 * t248;
t171 = t194 * t236 - t250 * t252 + t260 * t251;
t170 = t192 * t211 + t193 * t241 + t206 * t280;
t169 = -t192 * t252 + t193 * t236 + t206 * t248;
t168 = t190 * t211 + t191 * t241 - t202 * t280;
t167 = -t190 * t252 + t191 * t236 - t202 * t248;
t166 = t175 * t240 + t189 * t235;
t163 = t187 * t236 + t247 * t252 - t259 * t251;
t157 = t183 * t235 + t240 * t246;
t155 = t170 * t240 + t180 * t235;
t154 = t168 * t240 + t179 * t235;
t153 = t178 * t235 + t240 * t244;
t152 = t178 * t240 - t235 * t244;
t148 = t153 * t239 + t163 * t234;
t147 = -t153 * t234 + t163 * t239;
t1 = [(g(1) * t238 - g(2) * t243) * MDP(2) + (g(1) * t243 + g(2) * t238) * MDP(3) + (g(1) * t202 + g(2) * t206) * MDP(9) + (-g(1) * t201 - g(2) * t204) * MDP(10) + (-g(1) * t186 - g(2) * t187) * MDP(11) + (-g(1) * t185 + g(2) * t247) * MDP(12) + t286 * MDP(13) + (-g(1) * (-t238 * pkin(1) - t202 * pkin(2) + pkin(10) * t277) - g(2) * (t243 * pkin(1) - pkin(2) * t206 + pkin(10) * t278) + t286 * qJ(3)) * MDP(14) + (-g(1) * t162 - g(2) * t244) * MDP(20) + (g(1) * t161 + g(2) * t163) * MDP(21) + (-g(1) * t151 - g(2) * t153) * MDP(27) + (g(1) * t289 - g(2) * t152) * MDP(28) + (-g(1) * t292 - g(2) * t148) * MDP(34) + (g(1) * t293 - g(2) * t147) * MDP(35); (-g(1) * t193 - g(2) * t191 - g(3) * t198) * MDP(11) + (-g(1) * t192 - g(2) * t190 - g(3) * t197) * MDP(12) + (-g(1) * t170 - g(2) * t168 - g(3) * t175) * MDP(20) + (g(1) * t169 + g(2) * t167 + g(3) * t174) * MDP(21) + (-g(1) * t155 - g(2) * t154 - g(3) * t166) * MDP(27) + (-g(1) * (-t170 * t235 + t180 * t240) - g(2) * (-t168 * t235 + t179 * t240) - g(3) * (-t175 * t235 + t189 * t240)) * MDP(28) + (-g(1) * (t155 * t239 + t169 * t234) - g(2) * (t154 * t239 + t167 * t234) - g(3) * (t166 * t239 + t174 * t234)) * MDP(34) + (-g(1) * (-t155 * t234 + t169 * t239) - g(2) * (-t154 * t234 + t167 * t239) - g(3) * (-t166 * t234 + t174 * t239)) * MDP(35) + (MDP(14) * pkin(2) + MDP(9)) * (-g(1) * t204 + g(2) * t201 - g(3) * t212) + (-MDP(10) + (MDP(14) * qJ(3) + MDP(13)) * t228) * (g(1) * t206 - g(2) * t202 + g(3) * t215); (-g(1) * t259 + g(2) * t199 + g(3) * t260) * MDP(14); (g(1) * t244 - g(2) * t162 + g(3) * t246) * MDP(21) + (-g(1) * (-t163 * t275 + t234 * t244) - g(2) * (t161 * t275 - t162 * t234) - g(3) * (-t171 * t275 + t234 * t246)) * MDP(34) + (-g(1) * (t163 * t276 + t239 * t244) - g(2) * (-t161 * t276 - t162 * t239) - g(3) * (t171 * t276 + t239 * t246)) * MDP(35) + (t240 * MDP(27) - MDP(28) * t235 + MDP(20)) * (g(1) * t163 - g(2) * t161 + g(3) * t171); (g(1) * t153 - g(2) * t151 + g(3) * t157) * MDP(28) + (-MDP(34) * t239 + MDP(35) * t234 - MDP(27)) * (g(1) * t152 + g(2) * t289 + g(3) * (t183 * t240 - t235 * t246)); (-g(1) * t147 - g(2) * t293 - g(3) * (-t157 * t234 + t171 * t239)) * MDP(34) + (g(1) * t148 - g(2) * t292 - g(3) * (-t157 * t239 - t171 * t234)) * MDP(35);];
taug  = t1;
