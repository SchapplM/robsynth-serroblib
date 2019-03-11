% Calculate kinetic energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:27
% EndTime: 2019-03-09 15:56:31
% DurationCPUTime: 3.49s
% Computational Cost: add. (1777->364), mult. (3340->520), div. (0->0), fcn. (3620->10), ass. (0->167)
t331 = Icges(4,1) + Icges(5,1);
t330 = -Icges(4,4) + Icges(5,5);
t329 = Icges(5,4) + Icges(4,5);
t328 = Icges(4,2) + Icges(5,3);
t327 = -Icges(5,6) + Icges(4,6);
t326 = Icges(6,3) + Icges(4,3) + Icges(5,2);
t267 = sin(qJ(3));
t270 = cos(qJ(3));
t272 = cos(qJ(1));
t269 = sin(qJ(1));
t271 = cos(qJ(2));
t303 = t269 * t271;
t219 = t267 * t303 + t270 * t272;
t220 = -t267 * t272 + t270 * t303;
t268 = sin(qJ(2));
t305 = t268 * t269;
t325 = t328 * t219 + t330 * t220 - t327 * t305;
t302 = t271 * t272;
t221 = t267 * t302 - t269 * t270;
t222 = t269 * t267 + t270 * t302;
t304 = t268 * t272;
t324 = t328 * t221 + t330 * t222 - t327 * t304;
t323 = t330 * t219 + t331 * t220 + t329 * t305;
t322 = t330 * t221 + t331 * t222 + t329 * t304;
t321 = t327 * t271 + (t328 * t267 + t330 * t270) * t268;
t320 = -t329 * t271 + (t330 * t267 + t331 * t270) * t268;
t264 = sin(pkin(10));
t265 = cos(pkin(10));
t176 = t219 * t265 - t220 * t264;
t308 = t219 * t264;
t177 = t220 * t265 + t308;
t319 = Icges(6,5) * t177 + Icges(6,6) * t176 + t327 * t219 - t329 * t220 - t326 * t305;
t178 = t221 * t265 - t222 * t264;
t307 = t221 * t264;
t179 = t222 * t265 + t307;
t318 = Icges(6,5) * t179 + Icges(6,6) * t178 + t327 * t221 - t329 * t222 - t326 * t304;
t212 = (-t264 * t270 + t265 * t267) * t268;
t306 = t264 * t267;
t213 = (t265 * t270 + t306) * t268;
t317 = Icges(6,5) * t213 + Icges(6,6) * t212 + (t327 * t267 - t329 * t270) * t268 + t326 * t271;
t312 = pkin(5) * t265;
t311 = Icges(2,4) * t269;
t310 = Icges(3,4) * t268;
t309 = Icges(3,4) * t271;
t181 = pkin(3) * t220 + qJ(4) * t219;
t188 = pkin(4) * t220 - qJ(5) * t305;
t300 = -t181 - t188;
t182 = pkin(3) * t222 + qJ(4) * t221;
t189 = t222 * pkin(4) - qJ(5) * t304;
t299 = -t182 - t189;
t223 = (pkin(3) * t270 + qJ(4) * t267) * t268;
t229 = pkin(4) * t268 * t270 + qJ(5) * t271;
t298 = -t223 - t229;
t297 = qJD(3) * t268;
t296 = qJD(5) * t268;
t295 = qJD(6) * t268;
t294 = V_base(5) * pkin(6) + V_base(1);
t250 = qJD(2) * t269 + V_base(4);
t258 = V_base(6) + qJD(1);
t291 = t268 * pkin(9);
t218 = t272 * t297 + t250;
t290 = pkin(2) * t271 + pkin(8) * t268;
t249 = -qJD(2) * t272 + V_base(5);
t289 = rSges(3,1) * t271 - rSges(3,2) * t268;
t288 = Icges(3,1) * t271 - t310;
t287 = -Icges(3,2) * t268 + t309;
t286 = Icges(3,5) * t271 - Icges(3,6) * t268;
t217 = t269 * t297 + t249;
t248 = pkin(1) * t272 + t269 * pkin(7);
t285 = -V_base(4) * pkin(6) + t258 * t248 + V_base(2);
t247 = t269 * pkin(1) - pkin(7) * t272;
t284 = V_base(4) * t247 - t248 * V_base(5) + V_base(3);
t225 = t290 * t269;
t246 = pkin(2) * t268 - pkin(8) * t271;
t283 = t249 * t246 + (-t225 - t247) * t258 + t294;
t282 = (-Icges(3,3) * t272 + t269 * t286) * t249 + (Icges(3,3) * t269 + t272 * t286) * t250 + (Icges(3,5) * t268 + Icges(3,6) * t271) * t258;
t226 = t290 * t272;
t281 = t258 * t226 - t246 * t250 + t285;
t280 = qJD(4) * t221 + t217 * t223 + t283;
t279 = t250 * t225 - t226 * t249 + t284;
t242 = -qJD(3) * t271 + t258;
t278 = qJD(4) * t219 + t242 * t182 + t281;
t277 = qJD(4) * t268 * t267 + t218 * t181 + t279;
t276 = qJD(5) * t271 + t218 * t188 + t277;
t275 = t217 * t229 - t272 * t296 + t280;
t274 = t242 * t189 - t269 * t296 + t278;
t202 = -Icges(3,6) * t272 + t269 * t287;
t203 = Icges(3,6) * t269 + t272 * t287;
t206 = -Icges(3,5) * t272 + t269 * t288;
t207 = Icges(3,5) * t269 + t272 * t288;
t236 = Icges(3,2) * t271 + t310;
t239 = Icges(3,1) * t268 + t309;
t273 = (-t203 * t268 + t207 * t271) * t250 + (-t202 * t268 + t206 * t271) * t249 + (-t236 * t268 + t239 * t271) * t258;
t263 = pkin(10) + qJ(6);
t261 = Icges(2,4) * t272;
t257 = cos(t263);
t256 = sin(t263);
t245 = rSges(2,1) * t272 - t269 * rSges(2,2);
t244 = t269 * rSges(2,1) + rSges(2,2) * t272;
t243 = rSges(3,1) * t268 + rSges(3,2) * t271;
t241 = Icges(2,1) * t272 - t311;
t240 = Icges(2,1) * t269 + t261;
t238 = -Icges(2,2) * t269 + t261;
t237 = Icges(2,2) * t272 + t311;
t232 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t231 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t230 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t227 = (-qJD(3) + qJD(6)) * t271 + t258;
t211 = t269 * rSges(3,3) + t272 * t289;
t210 = -rSges(3,3) * t272 + t269 * t289;
t209 = -rSges(4,3) * t271 + (rSges(4,1) * t270 - rSges(4,2) * t267) * t268;
t208 = -rSges(5,2) * t271 + (rSges(5,1) * t270 + rSges(5,3) * t267) * t268;
t194 = (t256 * t267 + t257 * t270) * t268;
t193 = (-t256 * t270 + t257 * t267) * t268;
t192 = -t272 * t295 + t218;
t191 = -t269 * t295 + t217;
t187 = V_base(5) * rSges(2,3) - t244 * t258 + t294;
t186 = t245 * t258 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t184 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t175 = pkin(9) * t271 + (pkin(5) * t306 + t270 * t312) * t268;
t173 = t222 * rSges(4,1) - t221 * rSges(4,2) + rSges(4,3) * t304;
t172 = t222 * rSges(5,1) + rSges(5,2) * t304 + t221 * rSges(5,3);
t171 = rSges(4,1) * t220 - rSges(4,2) * t219 + rSges(4,3) * t305;
t170 = rSges(5,1) * t220 + rSges(5,2) * t305 + rSges(5,3) * t219;
t157 = t221 * t256 + t222 * t257;
t156 = t221 * t257 - t222 * t256;
t155 = t219 * t256 + t220 * t257;
t154 = t219 * t257 - t220 * t256;
t152 = rSges(6,1) * t213 + rSges(6,2) * t212 + rSges(6,3) * t271;
t151 = Icges(6,1) * t213 + Icges(6,4) * t212 + Icges(6,5) * t271;
t150 = Icges(6,4) * t213 + Icges(6,2) * t212 + Icges(6,6) * t271;
t147 = rSges(7,1) * t194 + rSges(7,2) * t193 + rSges(7,3) * t271;
t146 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t271;
t145 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t271;
t144 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t271;
t143 = t243 * t249 + (-t210 - t247) * t258 + t294;
t142 = t211 * t258 - t243 * t250 + t285;
t141 = pkin(5) * t307 + t222 * t312 - t272 * t291;
t140 = pkin(5) * t308 + t220 * t312 - t269 * t291;
t139 = t210 * t250 - t211 * t249 + t284;
t138 = t179 * rSges(6,1) + t178 * rSges(6,2) - rSges(6,3) * t304;
t137 = rSges(6,1) * t177 + rSges(6,2) * t176 - rSges(6,3) * t305;
t136 = Icges(6,1) * t179 + Icges(6,4) * t178 - Icges(6,5) * t304;
t135 = Icges(6,1) * t177 + Icges(6,4) * t176 - Icges(6,5) * t305;
t134 = Icges(6,4) * t179 + Icges(6,2) * t178 - Icges(6,6) * t304;
t133 = Icges(6,4) * t177 + Icges(6,2) * t176 - Icges(6,6) * t305;
t130 = t157 * rSges(7,1) + t156 * rSges(7,2) - rSges(7,3) * t304;
t129 = rSges(7,1) * t155 + rSges(7,2) * t154 - rSges(7,3) * t305;
t128 = Icges(7,1) * t157 + Icges(7,4) * t156 - Icges(7,5) * t304;
t127 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t305;
t126 = Icges(7,4) * t157 + Icges(7,2) * t156 - Icges(7,6) * t304;
t125 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t305;
t124 = Icges(7,5) * t157 + Icges(7,6) * t156 - Icges(7,3) * t304;
t123 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t305;
t122 = -t171 * t242 + t209 * t217 + t283;
t121 = t173 * t242 - t209 * t218 + t281;
t120 = t171 * t218 - t173 * t217 + t279;
t119 = t208 * t217 + (-t170 - t181) * t242 + t280;
t118 = t172 * t242 + (-t208 - t223) * t218 + t278;
t117 = t170 * t218 + (-t172 - t182) * t217 + t277;
t116 = t217 * t152 + (-t137 + t300) * t242 + t275;
t115 = t138 * t242 + (-t152 + t298) * t218 + t274;
t114 = t137 * t218 + (-t138 + t299) * t217 + t276;
t113 = -t227 * t129 + t191 * t147 + t217 * t175 + (-t140 + t300) * t242 + t275;
t112 = t130 * t227 + t141 * t242 - t147 * t192 + (-t175 + t298) * t218 + t274;
t111 = t129 * t192 - t130 * t191 + t140 * t218 + (-t141 + t299) * t217 + t276;
t1 = t191 * ((-t124 * t305 + t126 * t154 + t128 * t155) * t192 + (-t123 * t305 + t125 * t154 + t155 * t127) * t191 + (-t144 * t305 + t145 * t154 + t146 * t155) * t227) / 0.2e1 + t192 * ((-t124 * t304 + t156 * t126 + t157 * t128) * t192 + (-t123 * t304 + t156 * t125 + t157 * t127) * t191 + (-t144 * t304 + t156 * t145 + t157 * t146) * t227) / 0.2e1 + t250 * (t282 * t269 + t273 * t272) / 0.2e1 + t249 * (t273 * t269 - t282 * t272) / 0.2e1 + t227 * ((t124 * t271 + t126 * t193 + t128 * t194) * t192 + (t123 * t271 + t125 * t193 + t127 * t194) * t191 + (t271 * t144 + t193 * t145 + t194 * t146) * t227) / 0.2e1 + m(1) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(2) * (t184 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(3) * (t139 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + ((t203 * t271 + t207 * t268) * t250 + (t202 * t271 + t206 * t268) * t249 + (t271 * t236 + t268 * t239 + Icges(2,3)) * t258) * t258 / 0.2e1 + ((-t269 * t237 + t240 * t272 + Icges(1,4)) * V_base(5) + (-t269 * t238 + t272 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t272 * t237 + t269 * t240 + Icges(1,2)) * V_base(5) + (t238 * t272 + t269 * t241 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t150 * t176 + t151 * t177 + t219 * t321 + t220 * t320 - t305 * t317) * t242 + (t134 * t176 + t136 * t177 + t219 * t324 + t220 * t322 - t305 * t318) * t218 + (t176 * t133 + t177 * t135 + t219 * t325 + t323 * t220 - t319 * t305) * t217) * t217 / 0.2e1 + ((t178 * t150 + t179 * t151 + t221 * t321 + t222 * t320 - t304 * t317) * t242 + (t178 * t134 + t179 * t136 + t221 * t324 + t222 * t322 - t304 * t318) * t218 + (t178 * t133 + t179 * t135 + t221 * t325 + t323 * t222 - t319 * t304) * t217) * t218 / 0.2e1 + ((t134 * t212 + t136 * t213) * t218 + (t133 * t212 + t135 * t213) * t217 + (t212 * t150 + t213 * t151) * t242 + (t217 * t319 + t218 * t318 + t242 * t317) * t271 + ((t267 * t321 + t270 * t320) * t242 + (t267 * t324 + t270 * t322) * t218 + (t267 * t325 + t323 * t270) * t217) * t268) * t242 / 0.2e1 + V_base(4) * t258 * (Icges(2,5) * t272 - Icges(2,6) * t269) + V_base(5) * t258 * (Icges(2,5) * t269 + Icges(2,6) * t272) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
