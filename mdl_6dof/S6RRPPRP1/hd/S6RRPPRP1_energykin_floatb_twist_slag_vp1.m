% Calculate kinetic energy for
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:20
% EndTime: 2019-03-09 08:25:23
% DurationCPUTime: 3.13s
% Computational Cost: add. (2128->332), mult. (2215->456), div. (0->0), fcn. (2119->10), ass. (0->165)
t343 = Icges(6,1) + Icges(7,1);
t342 = -Icges(6,4) + Icges(7,5);
t341 = Icges(7,4) + Icges(6,5);
t340 = Icges(6,2) + Icges(7,3);
t339 = -Icges(7,6) + Icges(6,6);
t338 = Icges(3,3) + Icges(4,3);
t337 = -Icges(6,3) - Icges(7,2);
t252 = qJ(2) + pkin(9);
t243 = sin(t252);
t245 = cos(t252);
t257 = sin(qJ(2));
t259 = cos(qJ(2));
t336 = Icges(3,5) * t259 + Icges(4,5) * t245 - Icges(3,6) * t257 - Icges(4,6) * t243;
t335 = rSges(7,1) + pkin(5);
t334 = rSges(7,3) + qJ(6);
t251 = pkin(10) + qJ(5);
t244 = cos(t251);
t260 = cos(qJ(1));
t242 = sin(t251);
t258 = sin(qJ(1));
t302 = t258 * t242;
t189 = t244 * t260 + t245 * t302;
t301 = t258 * t244;
t190 = -t242 * t260 + t245 * t301;
t307 = t243 * t258;
t333 = t340 * t189 + t342 * t190 - t339 * t307;
t305 = t245 * t260;
t191 = t242 * t305 - t301;
t192 = t244 * t305 + t302;
t306 = t243 * t260;
t332 = t340 * t191 + t342 * t192 - t339 * t306;
t331 = -t339 * t189 + t341 * t190 - t337 * t307;
t330 = -t339 * t191 + t341 * t192 - t337 * t306;
t329 = t342 * t189 + t343 * t190 + t341 * t307;
t328 = t342 * t191 + t343 * t192 + t341 * t306;
t327 = t339 * t245 + (t340 * t242 + t342 * t244) * t243;
t326 = t337 * t245 + (-t339 * t242 + t341 * t244) * t243;
t325 = -t341 * t245 + (t342 * t242 + t343 * t244) * t243;
t308 = Icges(4,4) * t245;
t278 = -Icges(4,2) * t243 + t308;
t181 = -Icges(4,6) * t260 + t258 * t278;
t182 = Icges(4,6) * t258 + t260 * t278;
t309 = Icges(4,4) * t243;
t280 = Icges(4,1) * t245 - t309;
t183 = -Icges(4,5) * t260 + t258 * t280;
t184 = Icges(4,5) * t258 + t260 * t280;
t310 = Icges(3,4) * t259;
t279 = -Icges(3,2) * t257 + t310;
t195 = -Icges(3,6) * t260 + t258 * t279;
t196 = Icges(3,6) * t258 + t260 * t279;
t311 = Icges(3,4) * t257;
t281 = Icges(3,1) * t259 - t311;
t197 = -Icges(3,5) * t260 + t258 * t281;
t198 = Icges(3,5) * t258 + t260 * t281;
t211 = Icges(4,2) * t245 + t309;
t212 = Icges(4,1) * t243 + t308;
t224 = Icges(3,2) * t259 + t311;
t227 = Icges(3,1) * t257 + t310;
t238 = -qJD(2) * t260 + V_base(5);
t239 = qJD(2) * t258 + V_base(4);
t246 = V_base(6) + qJD(1);
t324 = (-t211 * t243 + t212 * t245 - t224 * t257 + t227 * t259) * t246 + (-t182 * t243 + t184 * t245 - t196 * t257 + t198 * t259) * t239 + (-t181 * t243 + t183 * t245 - t195 * t257 + t197 * t259) * t238;
t323 = (Icges(3,5) * t257 + Icges(4,5) * t243 + Icges(3,6) * t259 + Icges(4,6) * t245) * t246 + (t338 * t258 + t336 * t260) * t239 + (t336 * t258 - t338 * t260) * t238;
t316 = pkin(2) * t257;
t315 = pkin(2) * t259;
t254 = cos(pkin(10));
t314 = pkin(4) * t254;
t312 = Icges(2,4) * t258;
t253 = sin(pkin(10));
t304 = t253 * t260;
t303 = t254 * t260;
t300 = t258 * t253;
t299 = t258 * t254;
t297 = rSges(7,2) * t307 + t334 * t189 + t335 * t190;
t296 = rSges(7,2) * t306 + t334 * t191 + t335 * t192;
t295 = -rSges(7,2) * t245 + (t334 * t242 + t335 * t244) * t243;
t177 = -qJ(3) * t260 + t258 * t315;
t235 = t258 * pkin(1) - pkin(7) * t260;
t294 = -t177 - t235;
t178 = qJ(3) * t258 + t260 * t315;
t282 = pkin(3) * t245 + qJ(4) * t243;
t200 = t282 * t260;
t293 = -t178 - t200;
t292 = qJD(4) * t243;
t291 = qJD(5) * t243;
t290 = V_base(5) * pkin(6) + V_base(1);
t199 = t282 * t258;
t287 = -t199 + t294;
t213 = pkin(3) * t243 - qJ(4) * t245;
t286 = -t213 - t316;
t285 = qJD(3) * t258 + t238 * t316 + t290;
t284 = rSges(3,1) * t259 - rSges(3,2) * t257;
t283 = rSges(4,1) * t245 - rSges(4,2) * t243;
t236 = pkin(1) * t260 + t258 * pkin(7);
t275 = -V_base(4) * pkin(6) + t246 * t236 + V_base(2);
t274 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t273 = t238 * t213 + t260 * t292 + t285;
t272 = t239 * t177 + t274;
t269 = pkin(8) * t243 + t245 * t314;
t268 = -qJD(3) * t260 + t246 * t178 + t275;
t267 = -qJD(4) * t245 + t239 * t199 + t272;
t266 = t246 * t200 + t258 * t292 + t268;
t144 = -pkin(4) * t304 + t258 * t269;
t157 = -pkin(8) * t245 + t243 * t314;
t265 = t238 * t157 + (-t144 + t287) * t246 + t273;
t145 = pkin(4) * t300 + t260 * t269;
t264 = t239 * t144 + (-t145 + t293) * t238 + t267;
t263 = t246 * t145 + (-t157 + t286) * t239 + t266;
t249 = Icges(2,4) * t260;
t234 = rSges(2,1) * t260 - t258 * rSges(2,2);
t233 = t258 * rSges(2,1) + rSges(2,2) * t260;
t232 = rSges(3,1) * t257 + rSges(3,2) * t259;
t229 = Icges(2,1) * t260 - t312;
t228 = Icges(2,1) * t258 + t249;
t226 = -Icges(2,2) * t258 + t249;
t225 = Icges(2,2) * t260 + t312;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t217 = -qJD(5) * t245 + t246;
t214 = rSges(4,1) * t243 + rSges(4,2) * t245;
t208 = t260 * t291 + t239;
t207 = t258 * t291 + t238;
t206 = t245 * t303 + t300;
t205 = -t245 * t304 + t299;
t204 = t245 * t299 - t304;
t203 = -t245 * t300 - t303;
t202 = t258 * rSges(3,3) + t260 * t284;
t201 = -rSges(3,3) * t260 + t258 * t284;
t187 = t258 * rSges(4,3) + t260 * t283;
t186 = -rSges(4,3) * t260 + t258 * t283;
t175 = V_base(5) * rSges(2,3) - t233 * t246 + t290;
t174 = t234 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -rSges(5,3) * t245 + (rSges(5,1) * t254 - rSges(5,2) * t253) * t243;
t172 = -Icges(5,5) * t245 + (Icges(5,1) * t254 - Icges(5,4) * t253) * t243;
t171 = -Icges(5,6) * t245 + (Icges(5,4) * t254 - Icges(5,2) * t253) * t243;
t170 = -Icges(5,3) * t245 + (Icges(5,5) * t254 - Icges(5,6) * t253) * t243;
t168 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t166 = -rSges(6,3) * t245 + (rSges(6,1) * t244 - rSges(6,2) * t242) * t243;
t153 = t206 * rSges(5,1) + t205 * rSges(5,2) + rSges(5,3) * t306;
t152 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t307;
t151 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t306;
t150 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t307;
t149 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t306;
t148 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t307;
t147 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t306;
t146 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t307;
t142 = t192 * rSges(6,1) - t191 * rSges(6,2) + rSges(6,3) * t306;
t140 = rSges(6,1) * t190 - rSges(6,2) * t189 + rSges(6,3) * t307;
t125 = t232 * t238 + (-t201 - t235) * t246 + t290;
t124 = t202 * t246 - t232 * t239 + t275;
t123 = t201 * t239 - t202 * t238 + t274;
t122 = t214 * t238 + (-t186 + t294) * t246 + t285;
t121 = t246 * t187 + (-t214 - t316) * t239 + t268;
t120 = t186 * t239 + (-t178 - t187) * t238 + t272;
t119 = t173 * t238 + (-t152 + t287) * t246 + t273;
t118 = t246 * t153 + (-t173 + t286) * t239 + t266;
t117 = t152 * t239 + (-t153 + t293) * t238 + t267;
t116 = -t140 * t217 + t166 * t207 + t265;
t115 = t217 * t142 - t208 * t166 + t263;
t114 = t140 * t208 - t142 * t207 + t264;
t113 = qJD(6) * t191 + t207 * t295 - t217 * t297 + t265;
t112 = qJD(6) * t189 - t208 * t295 + t217 * t296 + t263;
t111 = qJD(6) * t242 * t243 - t207 * t296 + t208 * t297 + t264;
t1 = m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(2) * (t168 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + ((t189 * t327 + t190 * t325 + t307 * t326) * t217 + (t189 * t332 + t190 * t328 + t307 * t330) * t208 + (t333 * t189 + t329 * t190 + t331 * t307) * t207) * t207 / 0.2e1 + ((t191 * t327 + t192 * t325 + t306 * t326) * t217 + (t332 * t191 + t328 * t192 + t330 * t306) * t208 + (t191 * t333 + t329 * t192 + t331 * t306) * t207) * t208 / 0.2e1 + ((-t207 * t331 - t208 * t330 - t217 * t326) * t245 + ((t242 * t327 + t244 * t325) * t217 + (t242 * t332 + t244 * t328) * t208 + (t242 * t333 + t329 * t244) * t207) * t243) * t217 / 0.2e1 + ((-t258 * t225 + t228 * t260 + Icges(1,4)) * V_base(5) + (-t258 * t226 + t260 * t229 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t260 * t225 + t258 * t228 + Icges(1,2)) * V_base(5) + (t226 * t260 + t258 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t147 * t307 + t203 * t149 + t204 * t151) * t239 + (t146 * t307 + t203 * t148 + t204 * t150) * t238 + (t170 * t307 + t203 * t171 + t204 * t172) * t246 - t323 * t260 + t324 * t258) * t238 / 0.2e1 + ((t147 * t306 + t205 * t149 + t206 * t151) * t239 + (t146 * t306 + t205 * t148 + t206 * t150) * t238 + (t170 * t306 + t205 * t171 + t206 * t172) * t246 + t324 * t260 + t323 * t258) * t239 / 0.2e1 + ((t259 * t196 + t257 * t198 + (-t147 + t182) * t245 + (-t149 * t253 + t151 * t254 + t184) * t243) * t239 + (t259 * t195 + t257 * t197 + (-t146 + t181) * t245 + (-t148 * t253 + t150 * t254 + t183) * t243) * t238 + (t259 * t224 + t257 * t227 + Icges(2,3) + (-t170 + t211) * t245 + (-t171 * t253 + t172 * t254 + t212) * t243) * t246) * t246 / 0.2e1 + t246 * V_base(4) * (Icges(2,5) * t260 - Icges(2,6) * t258) + V_base(5) * t246 * (Icges(2,5) * t258 + Icges(2,6) * t260) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
