% Calculate Gravitation load on the joints for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:41
% EndTime: 2019-03-10 03:20:46
% DurationCPUTime: 1.86s
% Computational Cost: add. (1323->191), mult. (3663->308), div. (0->0), fcn. (4724->14), ass. (0->87)
t232 = cos(qJ(1));
t271 = cos(pkin(6));
t282 = cos(qJ(2));
t255 = t271 * t282;
t279 = sin(qJ(2));
t280 = sin(qJ(1));
t217 = -t232 * t255 + t279 * t280;
t254 = t271 * t279;
t218 = t232 * t254 + t280 * t282;
t229 = sin(qJ(3));
t226 = sin(pkin(6));
t269 = sin(pkin(7));
t260 = t226 * t269;
t256 = t232 * t260;
t270 = cos(pkin(7));
t258 = t229 * t270;
t281 = cos(qJ(3));
t188 = -t217 * t258 + t218 * t281 - t229 * t256;
t261 = t226 * t270;
t207 = t217 * t269 - t232 * t261;
t228 = sin(qJ(4));
t231 = cos(qJ(4));
t167 = t188 * t231 + t207 * t228;
t253 = t270 * t281;
t187 = t217 * t253 + t218 * t229 + t256 * t281;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t150 = t167 * t227 - t187 * t230;
t151 = t167 * t230 + t187 * t227;
t285 = MDP(24) - MDP(33);
t284 = MDP(30) + MDP(32);
t283 = MDP(31) - MDP(34);
t168 = -t188 * t228 + t207 * t231;
t235 = t232 * t279 + t255 * t280;
t286 = t235 * t270 - t280 * t260;
t278 = pkin(9) * t226;
t266 = t227 * t231;
t265 = t230 * t231;
t264 = pkin(10) * t269;
t263 = t226 * t282;
t262 = t226 * t279;
t259 = t228 * t269;
t257 = t231 * t269;
t252 = t270 * t279;
t251 = t269 * t279;
t248 = t271 * t269;
t246 = t226 * t251;
t219 = t232 * t282 - t280 * t254;
t192 = t219 * t281 - t286 * t229;
t233 = -t235 * t269 - t261 * t280;
t171 = t192 * t231 - t228 * t233;
t191 = t219 * t229 + t286 * t281;
t154 = t171 * t227 - t191 * t230;
t206 = t229 * t248 + (t258 * t282 + t281 * t279) * t226;
t216 = -t260 * t282 + t271 * t270;
t186 = t206 * t231 + t216 * t228;
t205 = t229 * t262 - t248 * t281 - t253 * t263;
t164 = t186 * t227 - t205 * t230;
t243 = g(1) * t154 + g(2) * t150 + g(3) * t164;
t215 = (-t229 * t252 + t281 * t282) * t226;
t214 = (t229 * t282 + t281 * t252) * t226;
t200 = t215 * t231 + t228 * t246;
t199 = t215 * t228 - t231 * t246;
t198 = -t219 * t258 - t235 * t281;
t197 = t219 * t253 - t229 * t235;
t196 = -t217 * t281 - t218 * t258;
t195 = -t217 * t229 + t218 * t253;
t179 = t200 * t230 + t214 * t227;
t178 = t200 * t227 - t214 * t230;
t177 = t198 * t231 + t219 * t259;
t176 = t198 * t228 - t219 * t257;
t175 = t196 * t231 + t218 * t259;
t174 = t196 * t228 - t218 * t257;
t173 = -t205 * t265 + t206 * t227;
t172 = -t205 * t266 - t206 * t230;
t170 = t192 * t228 + t231 * t233;
t165 = t186 * t230 + t205 * t227;
t163 = t177 * t230 + t197 * t227;
t162 = t177 * t227 - t197 * t230;
t161 = t175 * t230 + t195 * t227;
t160 = t175 * t227 - t195 * t230;
t159 = -t191 * t265 + t192 * t227;
t158 = -t191 * t266 - t192 * t230;
t157 = -t187 * t265 + t188 * t227;
t156 = -t187 * t266 - t188 * t230;
t155 = t171 * t230 + t191 * t227;
t1 = [(g(1) * t280 - g(2) * t232) * MDP(2) + (g(1) * t232 + g(2) * t280) * MDP(3) + (g(1) * t218 - g(2) * t219) * MDP(9) + (-g(1) * t217 + g(2) * t235) * MDP(10) + (g(1) * t188 - g(2) * t192) * MDP(16) + (-g(1) * t187 + g(2) * t191) * MDP(17) + (g(1) * t167 - g(2) * t171) * MDP(23) + (-g(1) * (-pkin(1) * t280 - t218 * pkin(2) - pkin(3) * t188 - pkin(4) * t167 - pkin(5) * t151 - pkin(11) * t187 + t168 * pkin(12) - qJ(6) * t150 + t232 * t278) - g(2) * (t232 * pkin(1) + t219 * pkin(2) + t192 * pkin(3) + t171 * pkin(4) + t155 * pkin(5) + t191 * pkin(11) + t170 * pkin(12) + t154 * qJ(6) + t278 * t280) + (g(1) * t207 + g(2) * t233) * pkin(10)) * MDP(35) + t283 * (-g(1) * t150 + g(2) * t154) + t285 * (g(1) * t168 + g(2) * t170) + t284 * (g(1) * t151 - g(2) * t155); (g(1) * t235 + g(2) * t217 - g(3) * t263) * MDP(9) + (g(1) * t219 + g(2) * t218 + g(3) * t262) * MDP(10) + (-g(1) * t198 - g(2) * t196 - g(3) * t215) * MDP(16) + (g(1) * t197 + g(2) * t195 + g(3) * t214) * MDP(17) + (-g(1) * t177 - g(2) * t175 - g(3) * t200) * MDP(23) + (-g(1) * (-pkin(2) * t235 + t198 * pkin(3) + t177 * pkin(4) + t163 * pkin(5) + t197 * pkin(11) + t176 * pkin(12) + t162 * qJ(6) + t219 * t264) - g(2) * (-t217 * pkin(2) + t196 * pkin(3) + t175 * pkin(4) + t161 * pkin(5) + t195 * pkin(11) + t174 * pkin(12) + t160 * qJ(6) + t218 * t264) - g(3) * (t215 * pkin(3) + t200 * pkin(4) + t179 * pkin(5) + t214 * pkin(11) + t199 * pkin(12) + t178 * qJ(6) + (pkin(2) * t282 + pkin(10) * t251) * t226)) * MDP(35) + t283 * (g(1) * t162 + g(2) * t160 + g(3) * t178) + t285 * (g(1) * t176 + g(2) * t174 + g(3) * t199) + t284 * (-g(1) * t163 - g(2) * t161 - g(3) * t179); (g(1) * t192 + g(2) * t188 + g(3) * t206) * MDP(17) + (-g(1) * (pkin(5) * t159 + pkin(11) * t192 + qJ(6) * t158) - g(2) * (pkin(5) * t157 + pkin(11) * t188 + qJ(6) * t156) - g(3) * (pkin(5) * t173 + pkin(11) * t206 + qJ(6) * t172)) * MDP(35) + t283 * (g(1) * t158 + g(2) * t156 + g(3) * t172) + t284 * (-g(1) * t159 - g(2) * t157 - g(3) * t173) + (MDP(16) + t231 * MDP(23) + (pkin(4) * t231 + pkin(12) * t228 + pkin(3)) * MDP(35) - t285 * t228) * (g(1) * t191 + g(2) * t187 + g(3) * t205); (-MDP(35) * pkin(12) + t285) * (g(1) * t171 + g(2) * t167 + g(3) * t186) + (MDP(23) + MDP(35) * (pkin(5) * t230 + qJ(6) * t227 + pkin(4)) + t284 * t230 - t283 * t227) * (-g(3) * (-t206 * t228 + t216 * t231) - g(2) * t168 + g(1) * t170); (-g(1) * (-pkin(5) * t154 + qJ(6) * t155) - g(2) * (-pkin(5) * t150 + qJ(6) * t151) - g(3) * (-pkin(5) * t164 + qJ(6) * t165)) * MDP(35) + t284 * t243 + t283 * (g(1) * t155 + g(2) * t151 + g(3) * t165); -t243 * MDP(35);];
taug  = t1;
