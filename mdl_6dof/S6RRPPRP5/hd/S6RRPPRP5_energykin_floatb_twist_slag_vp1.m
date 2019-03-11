% Calculate kinetic energy for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:20
% EndTime: 2019-03-09 08:41:24
% DurationCPUTime: 4.17s
% Computational Cost: add. (1456->316), mult. (2189->428), div. (0->0), fcn. (2093->8), ass. (0->158)
t355 = Icges(3,4) + Icges(4,6);
t354 = Icges(3,1) + Icges(4,2);
t353 = Icges(3,2) + Icges(4,3);
t256 = cos(qJ(2));
t352 = t355 * t256;
t254 = sin(qJ(2));
t351 = t355 * t254;
t350 = Icges(4,4) - Icges(3,5);
t349 = Icges(4,5) - Icges(3,6);
t348 = t353 * t254 - t352;
t347 = t354 * t256 - t351;
t346 = Icges(4,1) + Icges(3,3);
t345 = Icges(6,1) + Icges(7,1);
t344 = Icges(6,4) - Icges(7,5);
t343 = Icges(7,4) + Icges(6,5);
t342 = Icges(6,2) + Icges(7,3);
t341 = Icges(7,6) - Icges(6,6);
t340 = Icges(6,3) + Icges(7,2);
t255 = sin(qJ(1));
t257 = cos(qJ(1));
t339 = t348 * t255 - t349 * t257;
t338 = t349 * t255 + t348 * t257;
t337 = -t350 * t255 + t347 * t257;
t336 = t347 * t255 + t350 * t257;
t335 = -t353 * t256 - t351;
t334 = t354 * t254 + t352;
t333 = t349 * t254 - t350 * t256;
t332 = rSges(7,1) + pkin(5);
t331 = rSges(7,3) + qJ(6);
t250 = pkin(9) + qJ(5);
t244 = cos(t250);
t303 = t254 * t257;
t243 = sin(t250);
t304 = t243 * t255;
t175 = -t244 * t303 + t304;
t302 = t255 * t244;
t176 = t243 * t303 + t302;
t298 = t256 * t257;
t330 = t342 * t175 - t344 * t176 + t341 * t298;
t177 = t243 * t257 + t254 * t302;
t178 = -t244 * t257 + t254 * t304;
t299 = t255 * t256;
t329 = -t342 * t177 - t344 * t178 + t341 * t299;
t328 = t341 * t175 + t343 * t176 + t340 * t298;
t327 = -t341 * t177 + t343 * t178 + t340 * t299;
t326 = -t344 * t175 + t345 * t176 + t343 * t298;
t325 = t344 * t177 + t345 * t178 + t343 * t299;
t324 = (t344 * t243 + t342 * t244) * t256 + t341 * t254;
t323 = (-t343 * t243 + t341 * t244) * t256 + t340 * t254;
t322 = (-t345 * t243 - t344 * t244) * t256 + t343 * t254;
t211 = -pkin(3) * t257 + qJ(4) * t299;
t237 = qJD(2) * t255 + V_base(4);
t321 = qJD(4) * t254 + t237 * t211;
t236 = -qJD(2) * t257 + V_base(5);
t245 = V_base(6) + qJD(1);
t320 = (t335 * t254 + t334 * t256) * t245 + (t338 * t254 + t337 * t256) * t237 + (t339 * t254 + t336 * t256) * t236;
t319 = (-t350 * t254 - t349 * t256) * t245 + (t346 * t255 + t333 * t257) * t237 + (t333 * t255 - t346 * t257) * t236;
t251 = sin(pkin(9));
t312 = pkin(4) * t251;
t252 = cos(pkin(9));
t311 = pkin(4) * t252;
t310 = Icges(2,4) * t255;
t305 = qJ(4) * t254;
t301 = t255 * t251;
t300 = t255 * t252;
t296 = rSges(7,2) * t298 + t331 * t175 + t332 * t176;
t295 = rSges(7,2) * t299 - t331 * t177 + t332 * t178;
t294 = rSges(7,2) * t254 + (-t332 * t243 + t331 * t244) * t256;
t279 = pkin(2) * t256 + qJ(3) * t254;
t205 = t279 * t255;
t234 = t255 * pkin(1) - pkin(7) * t257;
t293 = -t205 - t234;
t206 = t279 * t257;
t210 = t255 * pkin(3) + qJ(4) * t298;
t292 = -t206 - t210;
t291 = qJD(3) * t254;
t290 = qJD(3) * t256;
t289 = qJD(4) * t256;
t288 = qJD(5) * t256;
t287 = V_base(5) * pkin(6) + V_base(1);
t284 = -t211 + t293;
t229 = pkin(2) * t254 - qJ(3) * t256;
t283 = -t229 - t305;
t282 = t236 * t229 + t257 * t291 + t287;
t281 = rSges(3,1) * t256 - rSges(3,2) * t254;
t280 = -rSges(4,2) * t256 + rSges(4,3) * t254;
t235 = pkin(1) * t257 + t255 * pkin(7);
t272 = -V_base(4) * pkin(6) + t245 * t235 + V_base(2);
t271 = V_base(4) * t234 - t235 * V_base(5) + V_base(3);
t270 = t236 * t305 + t257 * t289 + t282;
t269 = t237 * t205 + t271;
t268 = pkin(8) * t256 + t254 * t312;
t265 = t245 * t206 + t255 * t291 + t272;
t264 = t269 - t290;
t263 = t245 * t210 + t255 * t289 + t265;
t155 = t255 * t268 - t257 * t311;
t196 = pkin(8) * t254 - t256 * t312;
t262 = t236 * t196 + (-t155 + t284) * t245 + t270;
t154 = t255 * t311 + t257 * t268;
t261 = t237 * t155 + (-t154 + t292) * t236 + t269 + t321;
t260 = t245 * t154 + (-t196 + t283) * t237 + t263;
t248 = Icges(2,4) * t257;
t233 = rSges(2,1) * t257 - t255 * rSges(2,2);
t232 = t255 * rSges(2,1) + rSges(2,2) * t257;
t231 = rSges(3,1) * t254 + rSges(3,2) * t256;
t230 = -rSges(4,2) * t254 - rSges(4,3) * t256;
t228 = qJD(5) * t254 + t245;
t227 = Icges(2,1) * t257 - t310;
t226 = Icges(2,1) * t255 + t248;
t224 = -Icges(2,2) * t255 + t248;
t223 = Icges(2,2) * t257 + t310;
t214 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t213 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t212 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t204 = t257 * t288 + t237;
t203 = t255 * t288 + t236;
t202 = -t252 * t257 + t254 * t301;
t201 = t251 * t257 + t254 * t300;
t200 = t251 * t303 + t300;
t199 = t252 * t303 - t301;
t195 = -rSges(4,1) * t257 + t255 * t280;
t194 = t255 * rSges(4,1) + t257 * t280;
t193 = t255 * rSges(3,3) + t257 * t281;
t192 = -rSges(3,3) * t257 + t255 * t281;
t174 = rSges(5,3) * t254 + (-rSges(5,1) * t251 - rSges(5,2) * t252) * t256;
t172 = Icges(5,5) * t254 + (-Icges(5,1) * t251 - Icges(5,4) * t252) * t256;
t171 = Icges(5,6) * t254 + (-Icges(5,4) * t251 - Icges(5,2) * t252) * t256;
t170 = Icges(5,3) * t254 + (-Icges(5,5) * t251 - Icges(5,6) * t252) * t256;
t167 = rSges(6,3) * t254 + (-rSges(6,1) * t243 - rSges(6,2) * t244) * t256;
t159 = V_base(5) * rSges(2,3) - t232 * t245 + t287;
t158 = t233 * t245 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t156 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t153 = rSges(5,1) * t202 + rSges(5,2) * t201 + rSges(5,3) * t299;
t152 = t200 * rSges(5,1) + t199 * rSges(5,2) + rSges(5,3) * t298;
t151 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t299;
t150 = Icges(5,1) * t200 + Icges(5,4) * t199 + Icges(5,5) * t298;
t149 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t299;
t148 = Icges(5,4) * t200 + Icges(5,2) * t199 + Icges(5,6) * t298;
t147 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t299;
t146 = Icges(5,5) * t200 + Icges(5,6) * t199 + Icges(5,3) * t298;
t141 = rSges(6,1) * t178 + rSges(6,2) * t177 + rSges(6,3) * t299;
t139 = t176 * rSges(6,1) - t175 * rSges(6,2) + rSges(6,3) * t298;
t125 = t231 * t236 + (-t192 - t234) * t245 + t287;
t124 = t193 * t245 - t231 * t237 + t272;
t123 = t192 * t237 - t193 * t236 + t271;
t122 = t230 * t236 + (-t195 + t293) * t245 + t282;
t121 = t194 * t245 + (-t229 - t230) * t237 + t265;
t120 = t195 * t237 + (-t194 - t206) * t236 + t264;
t119 = t174 * t236 + (-t153 + t284) * t245 + t270;
t118 = t152 * t245 + (-t174 + t283) * t237 + t263;
t117 = t153 * t237 + (-t152 + t292) * t236 + t264 + t321;
t116 = -t141 * t228 + t167 * t203 + t262;
t115 = t139 * t228 - t167 * t204 + t260;
t114 = -t139 * t203 + t141 * t204 + t261 - t290;
t113 = qJD(6) * t175 + t203 * t294 - t228 * t295 + t262;
t112 = -qJD(6) * t177 - t204 * t294 + t228 * t296 + t260;
t111 = (qJD(6) * t244 - qJD(3)) * t256 + t295 * t204 - t296 * t203 + t261;
t1 = m(1) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(2) * (t156 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((-t177 * t324 + t178 * t322 + t299 * t323) * t228 + (-t177 * t330 + t326 * t178 + t328 * t299) * t204 + (-t329 * t177 + t325 * t178 + t327 * t299) * t203) * t203 / 0.2e1 + ((t175 * t324 + t176 * t322 + t298 * t323) * t228 + (t330 * t175 + t326 * t176 + t328 * t298) * t204 + (t175 * t329 + t176 * t325 + t298 * t327) * t203) * t204 / 0.2e1 + (((-t243 * t322 + t244 * t324) * t228 + (-t326 * t243 + t244 * t330) * t204 + (-t243 * t325 + t244 * t329) * t203) * t256 + (t203 * t327 + t204 * t328 + t228 * t323) * t254) * t228 / 0.2e1 + ((-t255 * t223 + t226 * t257 + Icges(1,4)) * V_base(5) + (-t255 * t224 + t227 * t257 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t223 * t257 + t255 * t226 + Icges(1,2)) * V_base(5) + (t224 * t257 + t255 * t227 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t146 * t299 + t148 * t201 + t150 * t202) * t237 + (t147 * t299 + t149 * t201 + t151 * t202) * t236 + (t170 * t299 + t171 * t201 + t172 * t202) * t245 - t319 * t257 + t320 * t255) * t236 / 0.2e1 + ((t146 * t298 + t199 * t148 + t200 * t150) * t237 + (t147 * t298 + t199 * t149 + t200 * t151) * t236 + (t170 * t298 + t199 * t171 + t200 * t172) * t245 + t320 * t257 + t319 * t255) * t237 / 0.2e1 + (((-t148 * t252 - t150 * t251 - t338) * t256 + (t146 + t337) * t254) * t237 + ((-t149 * t252 - t151 * t251 - t339) * t256 + (t147 + t336) * t254) * t236 + (Icges(2,3) + (-t171 * t252 - t172 * t251 - t335) * t256 + (t170 + t334) * t254) * t245) * t245 / 0.2e1 + t245 * V_base(4) * (Icges(2,5) * t257 - Icges(2,6) * t255) + V_base(5) * t245 * (Icges(2,5) * t255 + Icges(2,6) * t257) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
