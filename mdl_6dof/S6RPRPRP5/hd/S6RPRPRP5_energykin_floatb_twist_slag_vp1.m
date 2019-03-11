% Calculate kinetic energy for
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:22
% EndTime: 2019-03-09 03:14:25
% DurationCPUTime: 3.41s
% Computational Cost: add. (2102->338), mult. (2189->476), div. (0->0), fcn. (2093->10), ass. (0->168)
t338 = Icges(6,1) + Icges(7,1);
t337 = -Icges(6,4) + Icges(7,5);
t336 = Icges(7,4) + Icges(6,5);
t335 = Icges(6,2) + Icges(7,3);
t334 = -Icges(7,6) + Icges(6,6);
t333 = -Icges(6,3) - Icges(7,2);
t332 = rSges(7,1) + pkin(5);
t331 = rSges(7,3) + qJ(6);
t251 = pkin(10) + qJ(5);
t244 = cos(t251);
t252 = pkin(9) + qJ(3);
t245 = cos(t252);
t260 = cos(qJ(1));
t242 = sin(t251);
t259 = sin(qJ(1));
t302 = t259 * t242;
t189 = t244 * t260 + t245 * t302;
t301 = t259 * t244;
t190 = -t242 * t260 + t245 * t301;
t243 = sin(t252);
t307 = t243 * t259;
t330 = t335 * t189 + t337 * t190 - t334 * t307;
t305 = t245 * t260;
t191 = t242 * t305 - t301;
t192 = t244 * t305 + t302;
t306 = t243 * t260;
t329 = t335 * t191 + t337 * t192 - t334 * t306;
t328 = -t334 * t189 + t336 * t190 - t333 * t307;
t327 = -t334 * t191 + t336 * t192 - t333 * t306;
t326 = t337 * t189 + t338 * t190 + t336 * t307;
t325 = t337 * t191 + t338 * t192 + t336 * t306;
t324 = t334 * t245 + (t335 * t242 + t337 * t244) * t243;
t323 = t333 * t245 + (-t334 * t242 + t336 * t244) * t243;
t322 = -t336 * t245 + (t337 * t242 + t338 * t244) * t243;
t254 = sin(pkin(9));
t315 = pkin(2) * t254;
t256 = cos(pkin(9));
t314 = pkin(2) * t256;
t255 = cos(pkin(10));
t313 = pkin(4) * t255;
t312 = Icges(2,4) * t259;
t311 = Icges(3,4) * t254;
t310 = Icges(3,4) * t256;
t309 = Icges(4,4) * t243;
t308 = Icges(4,4) * t245;
t253 = sin(pkin(10));
t304 = t253 * t260;
t303 = t255 * t260;
t300 = t259 * t253;
t299 = t259 * t255;
t296 = rSges(7,2) * t307 + t331 * t189 + t190 * t332;
t295 = rSges(7,2) * t306 + t331 * t191 + t192 * t332;
t294 = -rSges(7,2) * t245 + (t331 * t242 + t244 * t332) * t243;
t176 = -pkin(7) * t260 + t259 * t314;
t232 = t259 * pkin(1) - qJ(2) * t260;
t293 = -t176 - t232;
t292 = qJD(4) * t243;
t291 = qJD(5) * t243;
t290 = V_base(4) * t232 + V_base(3);
t289 = V_base(5) * pkin(6) + V_base(1);
t281 = pkin(3) * t245 + qJ(4) * t243;
t201 = t281 * t259;
t286 = -t201 + t293;
t238 = qJD(3) * t259 + V_base(4);
t246 = V_base(6) + qJD(1);
t285 = qJD(2) * t259 + t289;
t284 = V_base(5) * t315 + t285;
t237 = -qJD(3) * t260 + V_base(5);
t283 = rSges(3,1) * t256 - rSges(3,2) * t254;
t282 = rSges(4,1) * t245 - rSges(4,2) * t243;
t280 = Icges(3,1) * t256 - t311;
t279 = Icges(4,1) * t245 - t309;
t278 = -Icges(3,2) * t254 + t310;
t277 = -Icges(4,2) * t243 + t308;
t276 = Icges(3,5) * t256 - Icges(3,6) * t254;
t275 = Icges(4,5) * t245 - Icges(4,6) * t243;
t234 = pkin(1) * t260 + t259 * qJ(2);
t274 = -qJD(2) * t260 + t246 * t234 + V_base(2);
t213 = pkin(3) * t243 - qJ(4) * t245;
t273 = t237 * t213 + t260 * t292 + t284;
t272 = (-Icges(4,3) * t260 + t259 * t275) * t237 + (Icges(4,3) * t259 + t260 * t275) * t238 + (Icges(4,5) * t243 + Icges(4,6) * t245) * t246;
t271 = pkin(8) * t243 + t245 * t313;
t177 = pkin(7) * t259 + t260 * t314;
t270 = V_base(4) * t176 + (-t177 - t234) * V_base(5) + t290;
t269 = (-Icges(3,3) * t260 + t259 * t276) * V_base(5) + (Icges(3,3) * t259 + t260 * t276) * V_base(4) + (Icges(3,5) * t254 + Icges(3,6) * t256) * t246;
t268 = t246 * t177 + (-pkin(6) - t315) * V_base(4) + t274;
t267 = -qJD(4) * t245 + t238 * t201 + t270;
t144 = -pkin(4) * t304 + t259 * t271;
t157 = -pkin(8) * t245 + t243 * t313;
t266 = t237 * t157 + (-t144 + t286) * t246 + t273;
t202 = t281 * t260;
t265 = t246 * t202 + t259 * t292 + t268;
t145 = pkin(4) * t300 + t260 * t271;
t264 = t238 * t144 + (-t145 - t202) * t237 + t267;
t263 = t246 * t145 + (-t157 - t213) * t238 + t265;
t181 = -Icges(4,6) * t260 + t259 * t277;
t182 = Icges(4,6) * t259 + t260 * t277;
t183 = -Icges(4,5) * t260 + t259 * t279;
t184 = Icges(4,5) * t259 + t260 * t279;
t211 = Icges(4,2) * t245 + t309;
t212 = Icges(4,1) * t243 + t308;
t262 = (-t182 * t243 + t184 * t245) * t238 + (-t181 * t243 + t183 * t245) * t237 + (-t211 * t243 + t212 * t245) * t246;
t195 = -Icges(3,6) * t260 + t259 * t278;
t196 = Icges(3,6) * t259 + t260 * t278;
t197 = -Icges(3,5) * t260 + t259 * t280;
t198 = Icges(3,5) * t259 + t260 * t280;
t221 = Icges(3,2) * t256 + t311;
t222 = Icges(3,1) * t254 + t310;
t261 = (-t196 * t254 + t198 * t256) * V_base(4) + (-t195 * t254 + t197 * t256) * V_base(5) + (-t221 * t254 + t222 * t256) * t246;
t249 = Icges(2,4) * t260;
t235 = rSges(2,1) * t260 - t259 * rSges(2,2);
t233 = t259 * rSges(2,1) + rSges(2,2) * t260;
t229 = Icges(2,1) * t260 - t312;
t228 = Icges(2,1) * t259 + t249;
t227 = -Icges(2,2) * t259 + t249;
t226 = Icges(2,2) * t260 + t312;
t225 = Icges(2,5) * t260 - Icges(2,6) * t259;
t224 = Icges(2,5) * t259 + Icges(2,6) * t260;
t223 = rSges(3,1) * t254 + rSges(3,2) * t256;
t219 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t218 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t217 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t216 = -qJD(5) * t245 + t246;
t214 = rSges(4,1) * t243 + rSges(4,2) * t245;
t208 = t260 * t291 + t238;
t207 = t259 * t291 + t237;
t206 = t245 * t303 + t300;
t205 = -t245 * t304 + t299;
t204 = t245 * t299 - t304;
t203 = -t245 * t300 - t303;
t200 = t259 * rSges(3,3) + t260 * t283;
t199 = -rSges(3,3) * t260 + t259 * t283;
t187 = t259 * rSges(4,3) + t260 * t282;
t186 = -rSges(4,3) * t260 + t259 * t282;
t175 = V_base(5) * rSges(2,3) - t233 * t246 + t289;
t174 = t235 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -rSges(5,3) * t245 + (rSges(5,1) * t255 - rSges(5,2) * t253) * t243;
t172 = -Icges(5,5) * t245 + (Icges(5,1) * t255 - Icges(5,4) * t253) * t243;
t171 = -Icges(5,6) * t245 + (Icges(5,4) * t255 - Icges(5,2) * t253) * t243;
t170 = -Icges(5,3) * t245 + (Icges(5,5) * t255 - Icges(5,6) * t253) * t243;
t168 = t233 * V_base(4) - t235 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t245 + (rSges(6,1) * t244 - rSges(6,2) * t242) * t243;
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
t125 = t223 * V_base(5) + (-t199 - t232) * t246 + t285;
t124 = t246 * t200 + (-pkin(6) - t223) * V_base(4) + t274;
t123 = t199 * V_base(4) + (-t200 - t234) * V_base(5) + t290;
t122 = t214 * t237 + (-t186 + t293) * t246 + t284;
t121 = t246 * t187 - t238 * t214 + t268;
t120 = t186 * t238 - t187 * t237 + t270;
t119 = t173 * t237 + (-t152 + t286) * t246 + t273;
t118 = t246 * t153 + (-t173 - t213) * t238 + t265;
t117 = t152 * t238 + (-t153 - t202) * t237 + t267;
t116 = -t140 * t216 + t165 * t207 + t266;
t115 = t216 * t142 - t208 * t165 + t263;
t114 = t140 * t208 - t142 * t207 + t264;
t113 = qJD(6) * t191 + t207 * t294 - t216 * t296 + t266;
t112 = qJD(6) * t189 - t208 * t294 + t216 * t295 + t263;
t111 = qJD(6) * t242 * t243 - t207 * t295 + t208 * t296 + t264;
t1 = m(2) * (t168 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(1) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + ((t189 * t324 + t190 * t322 + t307 * t323) * t216 + (t189 * t329 + t190 * t325 + t307 * t327) * t208 + (t189 * t330 + t326 * t190 + t328 * t307) * t207) * t207 / 0.2e1 + ((t191 * t324 + t192 * t322 + t306 * t323) * t216 + (t191 * t329 + t192 * t325 + t306 * t327) * t208 + (t191 * t330 + t326 * t192 + t328 * t306) * t207) * t208 / 0.2e1 + ((-t207 * t328 - t208 * t327 - t216 * t323) * t245 + ((t242 * t324 + t244 * t322) * t216 + (t242 * t329 + t244 * t325) * t208 + (t242 * t330 + t326 * t244) * t207) * t243) * t216 / 0.2e1 + ((t147 * t307 + t149 * t203 + t151 * t204) * t238 + (t146 * t307 + t203 * t148 + t204 * t150) * t237 + (t170 * t307 + t171 * t203 + t172 * t204) * t246 + t259 * t262 - t260 * t272) * t237 / 0.2e1 + ((t147 * t306 + t205 * t149 + t206 * t151) * t238 + (t146 * t306 + t205 * t148 + t206 * t150) * t237 + (t170 * t306 + t205 * t171 + t206 * t172) * t246 + t259 * t272 + t260 * t262) * t238 / 0.2e1 + (t225 * t246 + t259 * t269 + t260 * t261 + (-t259 * t226 + t228 * t260 + Icges(1,4)) * V_base(5) + (-t259 * t227 + t260 * t229 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t224 * t246 + t259 * t261 - t269 * t260 + (t260 * t226 + t259 * t228 + Icges(1,2)) * V_base(5) + (t227 * t260 + t259 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t146 * t237 - t147 * t238) * t245 + ((-t149 * t253 + t151 * t255) * t238 + (-t148 * t253 + t150 * t255) * t237) * t243 + (t182 * t245 + t184 * t243) * t238 + (t181 * t245 + t183 * t243) * t237 + (t195 * t256 + t197 * t254 + t224) * V_base(5) + (t196 * t256 + t198 * t254 + t225) * V_base(4) + (t256 * t221 + t254 * t222 + Icges(2,3) + (-t170 + t211) * t245 + (-t171 * t253 + t172 * t255 + t212) * t243) * t246) * t246 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
