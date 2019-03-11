% Calculate kinetic energy for
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:11:55
% EndTime: 2019-03-10 01:11:58
% DurationCPUTime: 3.41s
% Computational Cost: add. (2318->333), mult. (2384->495), div. (0->0), fcn. (2288->10), ass. (0->169)
t343 = Icges(6,1) + Icges(7,1);
t342 = -Icges(6,4) + Icges(7,5);
t341 = Icges(7,4) + Icges(6,5);
t340 = Icges(6,2) + Icges(7,3);
t339 = -Icges(7,6) + Icges(6,6);
t338 = -Icges(6,3) - Icges(7,2);
t337 = rSges(7,1) + pkin(5);
t336 = rSges(7,3) + qJ(6);
t260 = qJ(4) + qJ(5);
t256 = cos(t260);
t261 = qJ(2) + qJ(3);
t257 = cos(t261);
t267 = cos(qJ(1));
t254 = sin(t260);
t264 = sin(qJ(1));
t312 = t254 * t264;
t197 = t256 * t267 + t257 * t312;
t306 = t264 * t256;
t198 = -t254 * t267 + t257 * t306;
t255 = sin(t261);
t311 = t255 * t264;
t335 = t340 * t197 + t342 * t198 - t339 * t311;
t309 = t257 * t267;
t199 = t254 * t309 - t306;
t200 = t256 * t309 + t312;
t310 = t255 * t267;
t334 = t340 * t199 + t342 * t200 - t339 * t310;
t333 = -t339 * t197 + t341 * t198 - t338 * t311;
t332 = -t339 * t199 + t341 * t200 - t338 * t310;
t331 = t342 * t197 + t343 * t198 + t341 * t311;
t330 = t342 * t199 + t343 * t200 + t341 * t310;
t329 = t339 * t257 + (t340 * t254 + t342 * t256) * t255;
t328 = t338 * t257 + (-t339 * t254 + t341 * t256) * t255;
t327 = -t341 * t257 + (t342 * t254 + t343 * t256) * t255;
t263 = sin(qJ(2));
t322 = pkin(2) * t263;
t266 = cos(qJ(2));
t321 = pkin(2) * t266;
t265 = cos(qJ(4));
t320 = pkin(4) * t265;
t317 = Icges(2,4) * t264;
t316 = Icges(3,4) * t263;
t315 = Icges(3,4) * t266;
t314 = Icges(4,4) * t255;
t313 = Icges(4,4) * t257;
t262 = sin(qJ(4));
t308 = t262 * t264;
t307 = t262 * t267;
t305 = t264 * t265;
t304 = t265 * t267;
t303 = rSges(7,2) * t311 + t336 * t197 + t337 * t198;
t302 = rSges(7,2) * t310 + t336 * t199 + t337 * t200;
t301 = -rSges(7,2) * t257 + (t336 * t254 + t337 * t256) * t255;
t183 = -pkin(8) * t267 + t264 * t321;
t244 = t264 * pkin(1) - t267 * pkin(7);
t300 = -t183 - t244;
t299 = qJD(4) * t255;
t298 = qJD(5) * t255;
t297 = V_base(5) * pkin(6) + V_base(1);
t247 = qJD(2) * t264 + V_base(4);
t251 = V_base(6) + qJD(1);
t246 = -qJD(2) * t267 + V_base(5);
t294 = t246 * t322 + t297;
t223 = qJD(3) * t264 + t247;
t293 = pkin(3) * t257 + pkin(9) * t255;
t292 = rSges(3,1) * t266 - rSges(3,2) * t263;
t291 = rSges(4,1) * t257 - rSges(4,2) * t255;
t196 = t267 * t299 + t223;
t290 = Icges(3,1) * t266 - t316;
t289 = Icges(4,1) * t257 - t314;
t288 = -Icges(3,2) * t263 + t315;
t287 = -Icges(4,2) * t255 + t313;
t286 = Icges(3,5) * t266 - Icges(3,6) * t263;
t285 = Icges(4,5) * t257 - Icges(4,6) * t255;
t245 = t267 * pkin(1) + t264 * pkin(7);
t284 = -V_base(4) * pkin(6) + t251 * t245 + V_base(2);
t283 = V_base(4) * t244 - t245 * V_base(5) + V_base(3);
t222 = V_base(5) + (-qJD(2) - qJD(3)) * t267;
t195 = t264 * t299 + t222;
t282 = pkin(10) * t255 + t257 * t320;
t281 = (-Icges(4,3) * t267 + t264 * t285) * t222 + (Icges(4,3) * t264 + t267 * t285) * t223 + (Icges(4,5) * t255 + Icges(4,6) * t257) * t251;
t280 = (-Icges(3,3) * t267 + t264 * t286) * t246 + (Icges(3,3) * t264 + t267 * t286) * t247 + (Icges(3,5) * t263 + Icges(3,6) * t266) * t251;
t209 = t293 * t264;
t221 = pkin(3) * t255 - pkin(9) * t257;
t279 = t222 * t221 + (-t209 + t300) * t251 + t294;
t184 = pkin(8) * t264 + t267 * t321;
t278 = t247 * t183 - t184 * t246 + t283;
t277 = t251 * t184 - t247 * t322 + t284;
t149 = -pkin(4) * t307 + t264 * t282;
t162 = -pkin(10) * t257 + t255 * t320;
t229 = -qJD(4) * t257 + t251;
t276 = -t149 * t229 + t195 * t162 + t279;
t210 = t293 * t267;
t275 = t223 * t209 - t210 * t222 + t278;
t274 = t251 * t210 - t221 * t223 + t277;
t150 = pkin(4) * t308 + t267 * t282;
t273 = t196 * t149 - t150 * t195 + t275;
t272 = t229 * t150 - t162 * t196 + t274;
t188 = -Icges(4,6) * t267 + t264 * t287;
t189 = Icges(4,6) * t264 + t267 * t287;
t190 = -Icges(4,5) * t267 + t264 * t289;
t191 = Icges(4,5) * t264 + t267 * t289;
t218 = Icges(4,2) * t257 + t314;
t219 = Icges(4,1) * t255 + t313;
t271 = (-t189 * t255 + t191 * t257) * t223 + (-t188 * t255 + t190 * t257) * t222 + (-t218 * t255 + t219 * t257) * t251;
t203 = -Icges(3,6) * t267 + t264 * t288;
t204 = Icges(3,6) * t264 + t267 * t288;
t205 = -Icges(3,5) * t267 + t264 * t290;
t206 = Icges(3,5) * t264 + t267 * t290;
t233 = Icges(3,2) * t266 + t316;
t236 = Icges(3,1) * t263 + t315;
t270 = (-t204 * t263 + t206 * t266) * t247 + (-t203 * t263 + t205 * t266) * t246 + (-t233 * t263 + t236 * t266) * t251;
t258 = Icges(2,4) * t267;
t241 = rSges(2,1) * t267 - rSges(2,2) * t264;
t240 = rSges(2,1) * t264 + rSges(2,2) * t267;
t239 = rSges(3,1) * t263 + rSges(3,2) * t266;
t238 = Icges(2,1) * t267 - t317;
t237 = Icges(2,1) * t264 + t258;
t235 = -Icges(2,2) * t264 + t258;
t234 = Icges(2,2) * t267 + t317;
t228 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t227 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t226 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t220 = rSges(4,1) * t255 + rSges(4,2) * t257;
t215 = t257 * t304 + t308;
t214 = -t257 * t307 + t305;
t213 = t257 * t305 - t307;
t212 = -t257 * t308 - t304;
t211 = (-qJD(4) - qJD(5)) * t257 + t251;
t208 = rSges(3,3) * t264 + t267 * t292;
t207 = -rSges(3,3) * t267 + t264 * t292;
t193 = rSges(4,3) * t264 + t267 * t291;
t192 = -rSges(4,3) * t267 + t264 * t291;
t182 = -rSges(5,3) * t257 + (rSges(5,1) * t265 - rSges(5,2) * t262) * t255;
t181 = -Icges(5,5) * t257 + (Icges(5,1) * t265 - Icges(5,4) * t262) * t255;
t180 = -Icges(5,6) * t257 + (Icges(5,4) * t265 - Icges(5,2) * t262) * t255;
t179 = -Icges(5,3) * t257 + (Icges(5,5) * t265 - Icges(5,6) * t262) * t255;
t178 = V_base(5) * rSges(2,3) - t240 * t251 + t297;
t177 = t241 * t251 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t175 = t240 * V_base(4) - t241 * V_base(5) + V_base(3);
t174 = -rSges(6,3) * t257 + (rSges(6,1) * t256 - rSges(6,2) * t254) * t255;
t165 = t267 * t298 + t196;
t164 = t264 * t298 + t195;
t158 = rSges(5,1) * t215 + rSges(5,2) * t214 + rSges(5,3) * t310;
t157 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t311;
t156 = Icges(5,1) * t215 + Icges(5,4) * t214 + Icges(5,5) * t310;
t155 = Icges(5,1) * t213 + Icges(5,4) * t212 + Icges(5,5) * t311;
t154 = Icges(5,4) * t215 + Icges(5,2) * t214 + Icges(5,6) * t310;
t153 = Icges(5,4) * t213 + Icges(5,2) * t212 + Icges(5,6) * t311;
t152 = Icges(5,5) * t215 + Icges(5,6) * t214 + Icges(5,3) * t310;
t151 = Icges(5,5) * t213 + Icges(5,6) * t212 + Icges(5,3) * t311;
t147 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t310;
t145 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t311;
t130 = t239 * t246 + (-t207 - t244) * t251 + t297;
t129 = t208 * t251 - t239 * t247 + t284;
t127 = t207 * t247 - t208 * t246 + t283;
t126 = t220 * t222 + (-t192 + t300) * t251 + t294;
t125 = t193 * t251 - t220 * t223 + t277;
t124 = t192 * t223 - t193 * t222 + t278;
t123 = -t157 * t229 + t182 * t195 + t279;
t122 = t158 * t229 - t182 * t196 + t274;
t121 = t157 * t196 - t158 * t195 + t275;
t120 = -t145 * t211 + t164 * t174 + t276;
t119 = t147 * t211 - t165 * t174 + t272;
t118 = t145 * t165 - t147 * t164 + t273;
t117 = qJD(6) * t199 + t164 * t301 - t211 * t303 + t276;
t116 = qJD(6) * t197 - t165 * t301 + t211 * t302 + t272;
t115 = qJD(6) * t254 * t255 - t164 * t302 + t165 * t303 + t273;
t1 = t196 * ((t152 * t310 + t214 * t154 + t215 * t156) * t196 + (t151 * t310 + t153 * t214 + t155 * t215) * t195 + (t179 * t310 + t180 * t214 + t181 * t215) * t229) / 0.2e1 + t195 * ((t152 * t311 + t154 * t212 + t156 * t213) * t196 + (t151 * t311 + t212 * t153 + t213 * t155) * t195 + (t179 * t311 + t180 * t212 + t181 * t213) * t229) / 0.2e1 + t247 * (t280 * t264 + t270 * t267) / 0.2e1 + t246 * (t270 * t264 - t280 * t267) / 0.2e1 + t223 * (t281 * t264 + t271 * t267) / 0.2e1 + t222 * (t271 * t264 - t281 * t267) / 0.2e1 + t229 * ((-t151 * t195 - t152 * t196 - t179 * t229) * t257 + ((-t154 * t262 + t156 * t265) * t196 + (-t153 * t262 + t155 * t265) * t195 + (-t180 * t262 + t181 * t265) * t229) * t255) / 0.2e1 + m(1) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + m(2) * (t175 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(5) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(4) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(7) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + ((t197 * t329 + t198 * t327 + t311 * t328) * t211 + (t197 * t334 + t198 * t330 + t311 * t332) * t165 + (t335 * t197 + t331 * t198 + t333 * t311) * t164) * t164 / 0.2e1 + ((t199 * t329 + t200 * t327 + t310 * t328) * t211 + (t334 * t199 + t330 * t200 + t332 * t310) * t165 + (t199 * t335 + t331 * t200 + t333 * t310) * t164) * t165 / 0.2e1 + ((-t164 * t333 - t165 * t332 - t211 * t328) * t257 + ((t254 * t329 + t256 * t327) * t211 + (t254 * t334 + t256 * t330) * t165 + (t254 * t335 + t331 * t256) * t164) * t255) * t211 / 0.2e1 + ((-t234 * t264 + t237 * t267 + Icges(1,4)) * V_base(5) + (-t264 * t235 + t267 * t238 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t234 + t264 * t237 + Icges(1,2)) * V_base(5) + (t235 * t267 + t238 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t204 * t266 + t206 * t263) * t247 + (t203 * t266 + t205 * t263) * t246 + (t189 * t257 + t191 * t255) * t223 + (t188 * t257 + t190 * t255) * t222 + (t257 * t218 + t255 * t219 + t266 * t233 + t263 * t236 + Icges(2,3)) * t251) * t251 / 0.2e1 + t251 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t264) + t251 * V_base(5) * (Icges(2,5) * t264 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
