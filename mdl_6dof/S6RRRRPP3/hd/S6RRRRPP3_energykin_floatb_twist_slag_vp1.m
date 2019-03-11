% Calculate kinetic energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:39
% EndTime: 2019-03-09 20:53:41
% DurationCPUTime: 2.69s
% Computational Cost: add. (1877->285), mult. (2296->408), div. (0->0), fcn. (2216->8), ass. (0->145)
t317 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t316 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t315 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t314 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t313 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t312 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t311 = rSges(7,1) + pkin(5);
t310 = rSges(7,3) + qJ(6);
t242 = qJ(2) + qJ(3);
t239 = cos(t242);
t246 = cos(qJ(4));
t248 = cos(qJ(1));
t283 = t246 * t248;
t243 = sin(qJ(4));
t245 = sin(qJ(1));
t285 = t243 * t245;
t196 = t239 * t285 + t283;
t282 = t248 * t243;
t284 = t245 * t246;
t197 = t239 * t284 - t282;
t238 = sin(t242);
t288 = t238 * t245;
t309 = t316 * t196 + t197 * t317 - t315 * t288;
t198 = t239 * t282 - t284;
t199 = t239 * t283 + t285;
t286 = t238 * t248;
t308 = t316 * t198 + t199 * t317 - t315 * t286;
t307 = t196 * t314 + t197 * t316 - t288 * t313;
t306 = t198 * t314 + t199 * t316 - t286 * t313;
t305 = -t196 * t313 - t197 * t315 - t288 * t312;
t304 = -t198 * t313 - t199 * t315 - t286 * t312;
t303 = t312 * t239 + (-t243 * t313 - t246 * t315) * t238;
t302 = t313 * t239 + (t243 * t314 + t246 * t316) * t238;
t301 = t315 * t239 + (t316 * t243 + t246 * t317) * t238;
t244 = sin(qJ(2));
t296 = pkin(2) * t244;
t247 = cos(qJ(2));
t295 = pkin(2) * t247;
t293 = Icges(2,4) * t245;
t292 = Icges(3,4) * t244;
t291 = Icges(3,4) * t247;
t290 = Icges(4,4) * t238;
t289 = Icges(4,4) * t239;
t287 = t238 * t246;
t281 = rSges(7,2) * t196 + t310 * t197 + t311 * t288;
t280 = rSges(7,2) * t198 + t310 * t199 + t311 * t286;
t279 = (rSges(7,2) * t243 + rSges(7,3) * t246) * t238 + qJ(6) * t287 - t311 * t239;
t170 = -pkin(8) * t248 + t245 * t295;
t228 = t245 * pkin(1) - t248 * pkin(7);
t278 = -t170 - t228;
t277 = qJD(4) * t238;
t276 = V_base(5) * pkin(6) + V_base(1);
t231 = qJD(2) * t245 + V_base(4);
t235 = V_base(6) + qJD(1);
t230 = -qJD(2) * t248 + V_base(5);
t273 = t230 * t296 + t276;
t208 = qJD(3) * t245 + t231;
t272 = pkin(3) * t239 + pkin(9) * t238;
t271 = rSges(3,1) * t247 - rSges(3,2) * t244;
t270 = rSges(4,1) * t239 - rSges(4,2) * t238;
t269 = Icges(3,1) * t247 - t292;
t268 = Icges(4,1) * t239 - t290;
t267 = -Icges(3,2) * t244 + t291;
t266 = -Icges(4,2) * t238 + t289;
t265 = Icges(3,5) * t247 - Icges(3,6) * t244;
t264 = Icges(4,5) * t239 - Icges(4,6) * t238;
t229 = t248 * pkin(1) + t245 * pkin(7);
t263 = -V_base(4) * pkin(6) + t235 * t229 + V_base(2);
t262 = V_base(4) * t228 - t229 * V_base(5) + V_base(3);
t207 = V_base(5) + (-qJD(2) - qJD(3)) * t248;
t261 = (-Icges(4,3) * t248 + t245 * t264) * t207 + (Icges(4,3) * t245 + t248 * t264) * t208 + (Icges(4,5) * t238 + Icges(4,6) * t239) * t235;
t260 = (-Icges(3,3) * t248 + t245 * t265) * t230 + (Icges(3,3) * t245 + t248 * t265) * t231 + (Icges(3,5) * t244 + Icges(3,6) * t247) * t235;
t194 = t272 * t245;
t206 = pkin(3) * t238 - pkin(9) * t239;
t259 = t207 * t206 + (-t194 + t278) * t235 + t273;
t171 = pkin(8) * t245 + t248 * t295;
t258 = t231 * t170 - t171 * t230 + t262;
t257 = t235 * t171 - t231 * t296 + t263;
t181 = t245 * t277 + t207;
t193 = (pkin(4) * t246 + qJ(5) * t243) * t238;
t256 = qJD(5) * t198 + t181 * t193 + t259;
t195 = t272 * t248;
t255 = t208 * t194 - t195 * t207 + t258;
t254 = t235 * t195 - t206 * t208 + t257;
t147 = pkin(4) * t197 + qJ(5) * t196;
t182 = t248 * t277 + t208;
t253 = qJD(5) * t238 * t243 + t182 * t147 + t255;
t148 = pkin(4) * t199 + qJ(5) * t198;
t214 = -qJD(4) * t239 + t235;
t252 = qJD(5) * t196 + t214 * t148 + t254;
t175 = -Icges(4,6) * t248 + t245 * t266;
t176 = Icges(4,6) * t245 + t248 * t266;
t177 = -Icges(4,5) * t248 + t245 * t268;
t178 = Icges(4,5) * t245 + t248 * t268;
t203 = Icges(4,2) * t239 + t290;
t204 = Icges(4,1) * t238 + t289;
t251 = (-t176 * t238 + t178 * t239) * t208 + (-t175 * t238 + t177 * t239) * t207 + (-t203 * t238 + t204 * t239) * t235;
t185 = -Icges(3,6) * t248 + t245 * t267;
t186 = Icges(3,6) * t245 + t248 * t267;
t187 = -Icges(3,5) * t248 + t245 * t269;
t188 = Icges(3,5) * t245 + t248 * t269;
t218 = Icges(3,2) * t247 + t292;
t221 = Icges(3,1) * t244 + t291;
t250 = (-t186 * t244 + t188 * t247) * t231 + (-t185 * t244 + t187 * t247) * t230 + (-t218 * t244 + t221 * t247) * t235;
t240 = Icges(2,4) * t248;
t226 = rSges(2,1) * t248 - rSges(2,2) * t245;
t225 = rSges(2,1) * t245 + rSges(2,2) * t248;
t224 = rSges(3,1) * t244 + rSges(3,2) * t247;
t223 = Icges(2,1) * t248 - t293;
t222 = Icges(2,1) * t245 + t240;
t220 = -Icges(2,2) * t245 + t240;
t219 = Icges(2,2) * t248 + t293;
t213 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t212 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t211 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t205 = rSges(4,1) * t238 + rSges(4,2) * t239;
t192 = rSges(3,3) * t245 + t248 * t271;
t191 = -rSges(3,3) * t248 + t245 * t271;
t180 = rSges(4,3) * t245 + t248 * t270;
t179 = -rSges(4,3) * t248 + t245 * t270;
t169 = -rSges(6,1) * t239 + (-rSges(6,2) * t246 + rSges(6,3) * t243) * t238;
t167 = -rSges(5,3) * t239 + (rSges(5,1) * t246 - rSges(5,2) * t243) * t238;
t157 = V_base(5) * rSges(2,3) - t225 * t235 + t276;
t156 = t226 * t235 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t145 = rSges(5,1) * t199 - rSges(5,2) * t198 + rSges(5,3) * t286;
t144 = rSges(5,1) * t197 - rSges(5,2) * t196 + rSges(5,3) * t288;
t143 = rSges(6,1) * t286 - rSges(6,2) * t199 + rSges(6,3) * t198;
t141 = rSges(6,1) * t288 - rSges(6,2) * t197 + rSges(6,3) * t196;
t119 = t224 * t230 + (-t191 - t228) * t235 + t276;
t118 = t192 * t235 - t224 * t231 + t263;
t117 = t191 * t231 - t192 * t230 + t262;
t116 = t205 * t207 + (-t179 + t278) * t235 + t273;
t115 = t180 * t235 - t205 * t208 + t257;
t114 = t179 * t208 - t180 * t207 + t258;
t113 = -t144 * t214 + t167 * t181 + t259;
t112 = t145 * t214 - t167 * t182 + t254;
t111 = t144 * t182 - t145 * t181 + t255;
t110 = t169 * t181 + (-t141 - t147) * t214 + t256;
t109 = t143 * t214 + (-t169 - t193) * t182 + t252;
t108 = t141 * t182 + (-t143 - t148) * t181 + t253;
t107 = qJD(6) * t199 + t279 * t181 + (-t147 - t281) * t214 + t256;
t106 = qJD(6) * t197 + t280 * t214 + (-t193 - t279) * t182 + t252;
t105 = qJD(6) * t287 + t281 * t182 + (-t148 - t280) * t181 + t253;
t1 = t208 * (t245 * t261 + t248 * t251) / 0.2e1 + t207 * (t245 * t251 - t248 * t261) / 0.2e1 + t231 * (t245 * t260 + t248 * t250) / 0.2e1 + t230 * (t245 * t250 - t260 * t248) / 0.2e1 + m(1) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(2) * (t154 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + ((-t219 * t245 + t222 * t248 + Icges(1,4)) * V_base(5) + (-t220 * t245 + t223 * t248 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t219 * t248 + t222 * t245 + Icges(1,2)) * V_base(5) + (t220 * t248 + t223 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t196 * t302 + t197 * t301 + t288 * t303) * t214 + (t196 * t306 + t197 * t308 + t288 * t304) * t182 + (t307 * t196 + t309 * t197 + t305 * t288) * t181) * t181 / 0.2e1 + ((t198 * t302 + t199 * t301 + t286 * t303) * t214 + (t306 * t198 + t308 * t199 + t304 * t286) * t182 + (t307 * t198 + t199 * t309 + t305 * t286) * t181) * t182 / 0.2e1 + ((-t181 * t305 - t182 * t304 - t214 * t303) * t239 + ((t243 * t302 + t246 * t301) * t214 + (t243 * t306 + t246 * t308) * t182 + (t307 * t243 + t246 * t309) * t181) * t238) * t214 / 0.2e1 + ((t186 * t247 + t188 * t244) * t231 + (t185 * t247 + t187 * t244) * t230 + (t176 * t239 + t178 * t238) * t208 + (t175 * t239 + t177 * t238) * t207 + (t203 * t239 + t204 * t238 + t218 * t247 + t221 * t244 + Icges(2,3)) * t235) * t235 / 0.2e1 + t235 * V_base(4) * (Icges(2,5) * t248 - Icges(2,6) * t245) + t235 * V_base(5) * (Icges(2,5) * t245 + Icges(2,6) * t248) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
