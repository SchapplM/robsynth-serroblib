% Calculate kinetic energy for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:14
% EndTime: 2019-03-09 09:28:18
% DurationCPUTime: 3.78s
% Computational Cost: add. (1874->344), mult. (4158->485), div. (0->0), fcn. (4704->10), ass. (0->155)
t353 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t352 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t351 = Icges(3,5) + Icges(5,5) - Icges(4,4);
t350 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t349 = Icges(3,6) - Icges(4,5) - Icges(5,4);
t348 = Icges(3,3) + Icges(5,1) + Icges(4,1);
t297 = cos(pkin(6));
t301 = sin(qJ(1));
t303 = cos(qJ(2));
t326 = t301 * t303;
t300 = sin(qJ(2));
t304 = cos(qJ(1));
t327 = t300 * t304;
t263 = t297 * t326 + t327;
t325 = t303 * t304;
t328 = t300 * t301;
t264 = -t297 * t328 + t325;
t296 = sin(pkin(6));
t331 = t296 * t301;
t347 = -t352 * t263 + t264 * t353 + t351 * t331;
t261 = -t297 * t325 + t328;
t262 = t297 * t327 + t326;
t329 = t296 * t304;
t346 = -t352 * t261 + t262 * t353 - t351 * t329;
t345 = t263 * t350 - t264 * t352 - t331 * t349;
t344 = t261 * t350 - t262 * t352 + t329 * t349;
t343 = -t263 * t349 + t264 * t351 + t331 * t348;
t342 = -t261 * t349 + t262 * t351 - t329 * t348;
t341 = t348 * t297 + (t300 * t351 + t303 * t349) * t296;
t340 = t349 * t297 + (t300 * t352 + t303 * t350) * t296;
t339 = t351 * t297 + (t300 * t353 + t352 * t303) * t296;
t335 = cos(qJ(5));
t334 = pkin(8) * t297;
t333 = Icges(2,4) * t301;
t332 = t296 * t300;
t330 = t296 * t303;
t218 = pkin(2) * t262 + qJ(3) * t261;
t233 = -pkin(3) * t329 + t262 * qJ(4);
t324 = -t218 - t233;
t219 = pkin(2) * t264 + qJ(3) * t263;
t232 = pkin(3) * t331 + qJ(4) * t264;
t323 = -t219 - t232;
t265 = (pkin(2) * t300 - qJ(3) * t303) * t296;
t267 = pkin(3) * t297 + qJ(4) * t332;
t322 = -t265 - t267;
t321 = qJD(2) * t296;
t320 = V_base(5) * pkin(7) + V_base(1);
t317 = t296 * t335;
t275 = t301 * t321 + V_base(4);
t293 = V_base(6) + qJD(1);
t223 = -qJD(5) * t263 + t275;
t276 = qJD(2) * t297 + t293;
t257 = qJD(5) * t330 + t276;
t274 = -t304 * t321 + V_base(5);
t269 = t301 * pkin(1) - pkin(8) * t329;
t316 = -t269 * t293 + V_base(5) * t334 + t320;
t270 = pkin(1) * t304 + pkin(8) * t331;
t315 = V_base(4) * t269 - t270 * V_base(5) + V_base(3);
t222 = -qJD(5) * t261 + t274;
t314 = qJD(3) * t263 + t274 * t265 + t316;
t313 = t293 * t270 + V_base(2) + (-pkin(7) - t334) * V_base(4);
t312 = qJD(4) * t264 + t274 * t267 + t314;
t311 = qJD(3) * t261 + t276 * t219 + t313;
t310 = -qJD(3) * t330 + t275 * t218 + t315;
t309 = qJD(4) * t262 + t276 * t232 + t311;
t308 = qJD(4) * t332 + t275 * t233 + t310;
t235 = -pkin(4) * t329 - t261 * pkin(9);
t268 = pkin(4) * t297 + pkin(9) * t330;
t307 = t274 * t268 + (-t235 + t324) * t276 + t312;
t234 = pkin(4) * t331 - pkin(9) * t263;
t306 = t276 * t234 + (-t268 + t322) * t275 + t309;
t305 = t275 * t235 + (-t234 + t323) * t274 + t308;
t302 = cos(qJ(6));
t299 = sin(qJ(5));
t298 = sin(qJ(6));
t294 = Icges(2,4) * t304;
t284 = rSges(2,1) * t304 - t301 * rSges(2,2);
t283 = t301 * rSges(2,1) + rSges(2,2) * t304;
t282 = Icges(2,1) * t304 - t333;
t281 = Icges(2,1) * t301 + t294;
t280 = -Icges(2,2) * t301 + t294;
t279 = Icges(2,2) * t304 + t333;
t273 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t272 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t271 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t260 = t297 * t335 + t299 * t332;
t259 = t297 * t299 - t300 * t317;
t249 = rSges(4,1) * t297 + (-rSges(4,2) * t300 - rSges(4,3) * t303) * t296;
t248 = rSges(5,1) * t297 + (-rSges(5,2) * t303 + rSges(5,3) * t300) * t296;
t247 = rSges(3,3) * t297 + (rSges(3,1) * t300 + rSges(3,2) * t303) * t296;
t231 = V_base(5) * rSges(2,3) - t283 * t293 + t320;
t230 = t284 * t293 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t228 = t283 * V_base(4) - t284 * V_base(5) + V_base(3);
t227 = t262 * t299 - t304 * t317;
t226 = t262 * t335 + t299 * t329;
t225 = t264 * t299 + t301 * t317;
t224 = -t264 * t335 + t299 * t331;
t221 = t260 * t302 + t298 * t330;
t220 = -t260 * t298 + t302 * t330;
t216 = qJD(6) * t259 + t257;
t215 = pkin(5) * t260 + pkin(10) * t259;
t211 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t331;
t210 = t262 * rSges(3,1) - t261 * rSges(3,2) - rSges(3,3) * t329;
t209 = -rSges(4,1) * t329 - t262 * rSges(4,2) + t261 * rSges(4,3);
t208 = -rSges(5,1) * t329 + t261 * rSges(5,2) + t262 * rSges(5,3);
t207 = rSges(4,1) * t331 - rSges(4,2) * t264 + rSges(4,3) * t263;
t206 = rSges(5,1) * t331 + rSges(5,2) * t263 + rSges(5,3) * t264;
t187 = rSges(6,1) * t260 - rSges(6,2) * t259 + rSges(6,3) * t330;
t186 = Icges(6,1) * t260 - Icges(6,4) * t259 + Icges(6,5) * t330;
t185 = Icges(6,4) * t260 - Icges(6,2) * t259 + Icges(6,6) * t330;
t184 = Icges(6,5) * t260 - Icges(6,6) * t259 + Icges(6,3) * t330;
t181 = t227 * t302 - t261 * t298;
t180 = -t227 * t298 - t261 * t302;
t179 = t225 * t302 - t263 * t298;
t178 = -t225 * t298 - t263 * t302;
t177 = qJD(6) * t224 + t223;
t176 = -qJD(6) * t226 + t222;
t175 = pkin(5) * t227 - pkin(10) * t226;
t174 = pkin(5) * t225 + pkin(10) * t224;
t173 = rSges(6,1) * t227 + rSges(6,2) * t226 - rSges(6,3) * t261;
t172 = rSges(6,1) * t225 - rSges(6,2) * t224 - rSges(6,3) * t263;
t171 = Icges(6,1) * t227 + Icges(6,4) * t226 - Icges(6,5) * t261;
t170 = Icges(6,1) * t225 - Icges(6,4) * t224 - Icges(6,5) * t263;
t169 = Icges(6,4) * t227 + Icges(6,2) * t226 - Icges(6,6) * t261;
t168 = Icges(6,4) * t225 - Icges(6,2) * t224 - Icges(6,6) * t263;
t167 = Icges(6,5) * t227 + Icges(6,6) * t226 - Icges(6,3) * t261;
t166 = Icges(6,5) * t225 - Icges(6,6) * t224 - Icges(6,3) * t263;
t165 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t259;
t164 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t259;
t163 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t259;
t162 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t259;
t161 = -t210 * t276 + t247 * t274 + t316;
t160 = t211 * t276 - t247 * t275 + t313;
t159 = rSges(7,1) * t181 + rSges(7,2) * t180 - rSges(7,3) * t226;
t158 = rSges(7,1) * t179 + rSges(7,2) * t178 + rSges(7,3) * t224;
t157 = Icges(7,1) * t181 + Icges(7,4) * t180 - Icges(7,5) * t226;
t156 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t224;
t155 = Icges(7,4) * t181 + Icges(7,2) * t180 - Icges(7,6) * t226;
t154 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t224;
t153 = Icges(7,5) * t181 + Icges(7,6) * t180 - Icges(7,3) * t226;
t152 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t224;
t151 = t210 * t275 - t211 * t274 + t315;
t150 = t249 * t274 + (-t209 - t218) * t276 + t314;
t149 = t207 * t276 + (-t249 - t265) * t275 + t311;
t148 = t209 * t275 + (-t207 - t219) * t274 + t310;
t147 = t248 * t274 + (-t208 + t324) * t276 + t312;
t146 = t206 * t276 + (-t248 + t322) * t275 + t309;
t145 = t208 * t275 + (-t206 + t323) * t274 + t308;
t144 = -t173 * t257 + t187 * t222 + t307;
t143 = t172 * t257 - t187 * t223 + t306;
t142 = -t172 * t222 + t173 * t223 + t305;
t141 = -t159 * t216 + t165 * t176 - t175 * t257 + t215 * t222 + t307;
t140 = t158 * t216 - t165 * t177 + t174 * t257 - t215 * t223 + t306;
t139 = -t158 * t176 + t159 * t177 - t174 * t222 + t175 * t223 + t305;
t1 = t257 * ((t166 * t330 - t168 * t259 + t170 * t260) * t223 + (t167 * t330 - t169 * t259 + t171 * t260) * t222 + (t184 * t330 - t185 * t259 + t186 * t260) * t257) / 0.2e1 + m(3) * (t151 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(4) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(7) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(1) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + t223 * ((-t263 * t166 - t224 * t168 + t225 * t170) * t223 + (-t167 * t263 - t169 * t224 + t171 * t225) * t222 + (-t184 * t263 - t185 * t224 + t186 * t225) * t257) / 0.2e1 + t222 * ((-t166 * t261 + t168 * t226 + t170 * t227) * t223 + (-t261 * t167 + t226 * t169 + t227 * t171) * t222 + (-t184 * t261 + t185 * t226 + t186 * t227) * t257) / 0.2e1 + t176 * ((-t152 * t226 + t154 * t180 + t156 * t181) * t177 + (-t226 * t153 + t180 * t155 + t181 * t157) * t176 + (-t162 * t226 + t163 * t180 + t164 * t181) * t216) / 0.2e1 + m(2) * (t228 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t216 * ((t152 * t259 + t154 * t220 + t156 * t221) * t177 + (t153 * t259 + t155 * t220 + t157 * t221) * t176 + (t259 * t162 + t220 * t163 + t221 * t164) * t216) / 0.2e1 + t177 * ((t224 * t152 + t178 * t154 + t179 * t156) * t177 + (t153 * t224 + t155 * t178 + t157 * t179) * t176 + (t162 * t224 + t163 * t178 + t164 * t179) * t216) / 0.2e1 + ((-t301 * t279 + t281 * t304 + Icges(1,4)) * V_base(5) + (-t301 * t280 + t282 * t304 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t279 * t304 + t301 * t281 + Icges(1,2)) * V_base(5) + (t280 * t304 + t301 * t282 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t340 * t261 + t339 * t262 - t341 * t329) * t276 + (t345 * t261 + t347 * t262 - t343 * t329) * t275 + (t344 * t261 + t346 * t262 - t342 * t329) * t274) * t274 / 0.2e1 + ((-t340 * t263 + t339 * t264 + t341 * t331) * t276 + (t345 * t263 + t347 * t264 + t343 * t331) * t275 + (t344 * t263 + t346 * t264 + t342 * t331) * t274) * t275 / 0.2e1 + ((t342 * t274 + t343 * t275 + t341 * t276) * t297 + ((t339 * t300 + t340 * t303) * t276 + (t347 * t300 - t345 * t303) * t275 + (t346 * t300 - t344 * t303) * t274) * t296) * t276 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t301 + Icges(2,6) * t304) * V_base(5) + (Icges(2,5) * t304 - Icges(2,6) * t301) * V_base(4) + Icges(2,3) * t293 / 0.2e1) * t293;
T  = t1;
