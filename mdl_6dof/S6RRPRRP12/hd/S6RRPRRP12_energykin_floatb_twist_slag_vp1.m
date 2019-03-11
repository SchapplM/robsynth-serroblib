% Calculate kinetic energy for
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:17
% EndTime: 2019-03-09 12:50:21
% DurationCPUTime: 4.47s
% Computational Cost: add. (1519->313), mult. (2294->446), div. (0->0), fcn. (2198->8), ass. (0->156)
t351 = Icges(3,4) + Icges(4,6);
t350 = Icges(3,1) + Icges(4,2);
t349 = -Icges(3,2) - Icges(4,3);
t256 = cos(qJ(2));
t348 = t351 * t256;
t253 = sin(qJ(2));
t347 = t351 * t253;
t346 = Icges(4,4) - Icges(3,5);
t345 = Icges(4,5) - Icges(3,6);
t344 = t349 * t253 + t348;
t343 = t350 * t256 - t347;
t342 = Icges(4,1) + Icges(3,3);
t341 = Icges(6,1) + Icges(7,1);
t340 = Icges(6,4) - Icges(7,5);
t339 = Icges(7,4) + Icges(6,5);
t338 = Icges(6,2) + Icges(7,3);
t337 = Icges(7,6) - Icges(6,6);
t336 = Icges(6,3) + Icges(7,2);
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t335 = t254 * t344 + t257 * t345;
t334 = -t254 * t345 + t257 * t344;
t333 = t343 * t254 + t257 * t346;
t332 = -t254 * t346 + t343 * t257;
t331 = t349 * t256 - t347;
t330 = t350 * t253 + t348;
t329 = t345 * t253 - t256 * t346;
t328 = rSges(7,1) + pkin(5);
t327 = rSges(7,3) + qJ(6);
t251 = qJ(4) + qJ(5);
t247 = sin(t251);
t248 = cos(t251);
t299 = t253 * t257;
t191 = t247 * t254 - t248 * t299;
t192 = t247 * t299 + t248 * t254;
t295 = t256 * t257;
t326 = t338 * t191 - t340 * t192 + t337 * t295;
t300 = t253 * t254;
t193 = t247 * t257 + t248 * t300;
t194 = t247 * t300 - t248 * t257;
t297 = t254 * t256;
t325 = -t338 * t193 - t340 * t194 + t337 * t297;
t324 = t337 * t191 + t339 * t192 + t336 * t295;
t323 = -t337 * t193 + t339 * t194 + t336 * t297;
t322 = -t340 * t191 + t341 * t192 + t339 * t295;
t321 = t340 * t193 + t341 * t194 + t339 * t297;
t320 = (t340 * t247 + t338 * t248) * t256 + t337 * t253;
t319 = (-t339 * t247 + t337 * t248) * t256 + t336 * t253;
t318 = (-t341 * t247 - t340 * t248) * t256 + t339 * t253;
t237 = -qJD(2) * t257 + V_base(5);
t238 = qJD(2) * t254 + V_base(4);
t244 = V_base(6) + qJD(1);
t317 = (t331 * t253 + t330 * t256) * t244 + (-t334 * t253 + t332 * t256) * t238 + (-t335 * t253 + t333 * t256) * t237;
t316 = (-t253 * t346 - t345 * t256) * t244 + (t342 * t254 + t329 * t257) * t238 + (t329 * t254 - t342 * t257) * t237;
t252 = sin(qJ(4));
t309 = pkin(4) * t252;
t308 = t253 * pkin(8);
t255 = cos(qJ(4));
t307 = pkin(4) * t255;
t305 = Icges(2,4) * t254;
t298 = t254 * t255;
t296 = t255 * t257;
t294 = rSges(7,2) * t295 + t327 * t191 + t328 * t192;
t293 = rSges(7,2) * t297 - t327 * t193 + t328 * t194;
t292 = rSges(7,2) * t253 + (-t328 * t247 + t327 * t248) * t256;
t280 = pkin(2) * t256 + qJ(3) * t253;
t205 = t280 * t254;
t235 = pkin(1) * t254 - pkin(7) * t257;
t291 = -t205 - t235;
t290 = qJD(3) * t253;
t289 = qJD(3) * t256;
t288 = qJD(4) * t256;
t287 = qJD(5) * t256;
t286 = V_base(5) * pkin(6) + V_base(1);
t200 = t257 * t288 + t238;
t229 = qJD(4) * t253 + t244;
t230 = pkin(2) * t253 - qJ(3) * t256;
t283 = t237 * t230 + t257 * t290 + t286;
t282 = rSges(3,1) * t256 - rSges(3,2) * t253;
t281 = -rSges(4,2) * t256 + rSges(4,3) * t253;
t199 = t254 * t288 + t237;
t236 = pkin(1) * t257 + pkin(7) * t254;
t273 = -V_base(4) * pkin(6) + t244 * t236 + V_base(2);
t272 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t271 = pkin(9) * t256 + t253 * t309;
t270 = t238 * t205 + t272;
t206 = t280 * t257;
t267 = t244 * t206 + t254 * t290 + t273;
t212 = -t257 * pkin(3) + pkin(8) * t297;
t266 = t237 * t308 + (-t212 + t291) * t244 + t283;
t211 = t254 * pkin(3) + pkin(8) * t295;
t265 = t238 * t212 + (-t206 - t211) * t237 + t270;
t153 = t254 * t271 - t257 * t307;
t197 = pkin(9) * t253 - t256 * t309;
t264 = -t153 * t229 + t199 * t197 + t266;
t263 = t244 * t211 + (-t230 - t308) * t238 + t267;
t152 = t254 * t307 + t257 * t271;
t262 = -t152 * t199 + t200 * t153 + t265;
t261 = t229 * t152 - t197 * t200 + t263;
t249 = Icges(2,4) * t257;
t234 = rSges(2,1) * t257 - rSges(2,2) * t254;
t233 = rSges(2,1) * t254 + rSges(2,2) * t257;
t232 = rSges(3,1) * t253 + rSges(3,2) * t256;
t231 = -rSges(4,2) * t253 - rSges(4,3) * t256;
t228 = Icges(2,1) * t257 - t305;
t227 = Icges(2,1) * t254 + t249;
t225 = -Icges(2,2) * t254 + t249;
t224 = Icges(2,2) * t257 + t305;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t208 = qJD(5) * t253 + t229;
t204 = t252 * t300 - t296;
t203 = t252 * t257 + t253 * t298;
t202 = t252 * t299 + t298;
t201 = -t252 * t254 + t253 * t296;
t190 = -rSges(4,1) * t257 + t254 * t281;
t189 = rSges(4,1) * t254 + t257 * t281;
t188 = rSges(3,3) * t254 + t257 * t282;
t187 = rSges(5,3) * t253 + (-rSges(5,1) * t252 - rSges(5,2) * t255) * t256;
t186 = -rSges(3,3) * t257 + t254 * t282;
t177 = Icges(5,5) * t253 + (-Icges(5,1) * t252 - Icges(5,4) * t255) * t256;
t174 = Icges(5,6) * t253 + (-Icges(5,4) * t252 - Icges(5,2) * t255) * t256;
t171 = Icges(5,3) * t253 + (-Icges(5,5) * t252 - Icges(5,6) * t255) * t256;
t168 = t257 * t287 + t200;
t167 = t254 * t287 + t199;
t166 = rSges(6,3) * t253 + (-rSges(6,1) * t247 - rSges(6,2) * t248) * t256;
t157 = V_base(5) * rSges(2,3) - t233 * t244 + t286;
t156 = t234 * t244 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t151 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t297;
t150 = rSges(5,1) * t202 + rSges(5,2) * t201 + rSges(5,3) * t295;
t149 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t297;
t148 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t295;
t147 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t297;
t146 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t295;
t145 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t297;
t144 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t295;
t140 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t297;
t138 = rSges(6,1) * t192 - rSges(6,2) * t191 + rSges(6,3) * t295;
t123 = t232 * t237 + (-t186 - t235) * t244 + t286;
t122 = t188 * t244 - t232 * t238 + t273;
t121 = t186 * t238 - t188 * t237 + t272;
t120 = t231 * t237 + (-t190 + t291) * t244 + t283;
t119 = t189 * t244 + (-t230 - t231) * t238 + t267;
t118 = -t289 + t190 * t238 + (-t189 - t206) * t237 + t270;
t117 = -t151 * t229 + t187 * t199 + t266;
t116 = t150 * t229 - t187 * t200 + t263;
t115 = -t150 * t199 + t151 * t200 + t265 - t289;
t114 = -t140 * t208 + t166 * t167 + t264;
t113 = t138 * t208 - t166 * t168 + t261;
t112 = -t138 * t167 + t140 * t168 + t262 - t289;
t111 = qJD(6) * t191 + t167 * t292 - t208 * t293 + t264;
t110 = -qJD(6) * t193 - t168 * t292 + t208 * t294 + t261;
t109 = (qJD(6) * t248 - qJD(3)) * t256 + t293 * t168 - t294 * t167 + t262;
t1 = t229 * ((t144 * t200 + t145 * t199 + t171 * t229) * t253 + ((-t146 * t255 - t148 * t252) * t200 + (-t147 * t255 - t149 * t252) * t199 + (-t174 * t255 - t177 * t252) * t229) * t256) / 0.2e1 + m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t200 * ((t144 * t295 + t201 * t146 + t202 * t148) * t200 + (t145 * t295 + t147 * t201 + t149 * t202) * t199 + (t171 * t295 + t174 * t201 + t177 * t202) * t229) / 0.2e1 + t199 * ((t144 * t297 + t146 * t203 + t148 * t204) * t200 + (t145 * t297 + t203 * t147 + t204 * t149) * t199 + (t171 * t297 + t174 * t203 + t177 * t204) * t229) / 0.2e1 + ((-t193 * t320 + t194 * t318 + t297 * t319) * t208 + (-t193 * t326 + t322 * t194 + t324 * t297) * t168 + (-t325 * t193 + t321 * t194 + t323 * t297) * t167) * t167 / 0.2e1 + ((t191 * t320 + t192 * t318 + t295 * t319) * t208 + (t326 * t191 + t322 * t192 + t324 * t295) * t168 + (t191 * t325 + t192 * t321 + t295 * t323) * t167) * t168 / 0.2e1 + (((-t247 * t318 + t248 * t320) * t208 + (-t322 * t247 + t248 * t326) * t168 + (-t247 * t321 + t248 * t325) * t167) * t256 + (t167 * t323 + t168 * t324 + t208 * t319) * t253) * t208 / 0.2e1 + (t317 * t254 - t316 * t257) * t237 / 0.2e1 + (t316 * t254 + t317 * t257) * t238 / 0.2e1 + ((-t224 * t254 + t227 * t257 + Icges(1,4)) * V_base(5) + (-t254 * t225 + t257 * t228 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t257 * t224 + t254 * t227 + Icges(1,2)) * V_base(5) + (t225 * t257 + t228 * t254 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t332 * t253 + t334 * t256) * t238 + (t333 * t253 + t335 * t256) * t237 + (t330 * t253 - t331 * t256 + Icges(2,3)) * t244) * t244 / 0.2e1 + t244 * V_base(4) * (Icges(2,5) * t257 - Icges(2,6) * t254) + t244 * V_base(5) * (Icges(2,5) * t254 + Icges(2,6) * t257) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
