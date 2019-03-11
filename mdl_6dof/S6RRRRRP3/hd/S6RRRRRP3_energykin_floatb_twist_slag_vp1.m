% Calculate kinetic energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:46
% EndTime: 2019-03-10 01:05:49
% DurationCPUTime: 3.47s
% Computational Cost: add. (2365->338), mult. (2417->504), div. (0->0), fcn. (2311->10), ass. (0->170)
t335 = Icges(6,1) + Icges(7,1);
t334 = Icges(6,4) + Icges(7,4);
t333 = -Icges(7,5) - Icges(6,5);
t332 = Icges(6,2) + Icges(7,2);
t331 = -Icges(7,6) - Icges(6,6);
t330 = -Icges(7,3) - Icges(6,3);
t249 = qJ(4) + qJ(5);
t241 = sin(t249);
t243 = cos(t249);
t256 = cos(qJ(1));
t250 = qJ(2) + qJ(3);
t244 = cos(t250);
t253 = sin(qJ(1));
t303 = t244 * t253;
t183 = -t241 * t303 - t243 * t256;
t184 = -t241 * t256 + t243 * t303;
t242 = sin(t250);
t305 = t242 * t253;
t329 = -t331 * t183 - t333 * t184 - t330 * t305;
t302 = t244 * t256;
t185 = -t241 * t302 + t243 * t253;
t186 = t241 * t253 + t243 * t302;
t304 = t242 * t256;
t328 = -t331 * t185 - t333 * t186 - t330 * t304;
t327 = t332 * t183 + t334 * t184 - t331 * t305;
t326 = t332 * t185 + t334 * t186 - t331 * t304;
t325 = t334 * t183 + t335 * t184 - t333 * t305;
t324 = t334 * t185 + t335 * t186 - t333 * t304;
t323 = t330 * t244 + (t331 * t241 - t333 * t243) * t242;
t322 = t331 * t244 + (-t332 * t241 + t334 * t243) * t242;
t321 = t333 * t244 + (-t334 * t241 + t335 * t243) * t242;
t252 = sin(qJ(2));
t316 = pkin(2) * t252;
t255 = cos(qJ(2));
t314 = pkin(2) * t255;
t254 = cos(qJ(4));
t313 = t254 * pkin(4);
t310 = Icges(2,4) * t253;
t309 = Icges(3,4) * t252;
t308 = Icges(3,4) * t255;
t307 = Icges(4,4) * t242;
t306 = Icges(4,4) * t244;
t251 = sin(qJ(4));
t301 = t251 * t253;
t300 = t251 * t256;
t299 = t253 * t254;
t298 = t254 * t256;
t293 = pkin(5) * t243;
t269 = qJ(6) * t242 + t244 * t293;
t285 = pkin(5) * t241;
t297 = rSges(7,1) * t184 + rSges(7,2) * t183 + rSges(7,3) * t305 + t253 * t269 - t256 * t285;
t296 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t304 + t253 * t285 + t256 * t269;
t295 = (-qJ(6) - rSges(7,3)) * t244 + (rSges(7,1) * t243 - rSges(7,2) * t241 + t293) * t242;
t170 = -pkin(8) * t256 + t253 * t314;
t232 = t253 * pkin(1) - t256 * pkin(7);
t294 = -t170 - t232;
t291 = qJD(4) * t242;
t290 = qJD(5) * t242;
t289 = qJD(6) * t242;
t288 = V_base(5) * pkin(6) + V_base(1);
t235 = qJD(2) * t253 + V_base(4);
t238 = V_base(6) + qJD(1);
t234 = -qJD(2) * t256 + V_base(5);
t284 = t234 * t316 + t288;
t209 = qJD(3) * t253 + t235;
t283 = pkin(3) * t244 + pkin(9) * t242;
t282 = rSges(3,1) * t255 - rSges(3,2) * t252;
t281 = rSges(4,1) * t244 - rSges(4,2) * t242;
t182 = t256 * t291 + t209;
t280 = Icges(3,1) * t255 - t309;
t279 = Icges(4,1) * t244 - t307;
t278 = -Icges(3,2) * t252 + t308;
t277 = -Icges(4,2) * t242 + t306;
t276 = Icges(3,5) * t255 - Icges(3,6) * t252;
t275 = Icges(4,5) * t244 - Icges(4,6) * t242;
t233 = t256 * pkin(1) + t253 * pkin(7);
t274 = -V_base(4) * pkin(6) + t238 * t233 + V_base(2);
t273 = V_base(4) * t232 - t233 * V_base(5) + V_base(3);
t208 = V_base(5) + (-qJD(2) - qJD(3)) * t256;
t181 = t253 * t291 + t208;
t272 = pkin(10) * t242 + t244 * t313;
t271 = (-Icges(4,3) * t256 + t253 * t275) * t208 + (Icges(4,3) * t253 + t256 * t275) * t209 + (Icges(4,5) * t242 + Icges(4,6) * t244) * t238;
t270 = (-Icges(3,3) * t256 + t253 * t276) * t234 + (Icges(3,3) * t253 + t256 * t276) * t235 + (Icges(3,5) * t252 + Icges(3,6) * t255) * t238;
t195 = t283 * t253;
t207 = pkin(3) * t242 - pkin(9) * t244;
t268 = t208 * t207 + (-t195 + t294) * t238 + t284;
t171 = pkin(8) * t253 + t256 * t314;
t267 = t235 * t170 - t171 * t234 + t273;
t266 = t238 * t171 - t235 * t316 + t274;
t137 = -pkin(4) * t300 + t253 * t272;
t149 = -pkin(10) * t244 + t242 * t313;
t216 = -qJD(4) * t244 + t238;
t265 = -t137 * t216 + t181 * t149 + t268;
t196 = t283 * t256;
t264 = t209 * t195 - t196 * t208 + t267;
t263 = t238 * t196 - t207 * t209 + t266;
t138 = pkin(4) * t301 + t256 * t272;
t262 = t182 * t137 - t138 * t181 + t264;
t261 = t216 * t138 - t149 * t182 + t263;
t175 = -Icges(4,6) * t256 + t253 * t277;
t176 = Icges(4,6) * t253 + t256 * t277;
t177 = -Icges(4,5) * t256 + t253 * t279;
t178 = Icges(4,5) * t253 + t256 * t279;
t204 = Icges(4,2) * t244 + t307;
t205 = Icges(4,1) * t242 + t306;
t260 = (-t176 * t242 + t178 * t244) * t209 + (-t175 * t242 + t177 * t244) * t208 + (-t204 * t242 + t205 * t244) * t238;
t189 = -Icges(3,6) * t256 + t253 * t278;
t190 = Icges(3,6) * t253 + t256 * t278;
t191 = -Icges(3,5) * t256 + t253 * t280;
t192 = Icges(3,5) * t253 + t256 * t280;
t221 = Icges(3,2) * t255 + t309;
t224 = Icges(3,1) * t252 + t308;
t259 = (-t190 * t252 + t192 * t255) * t235 + (-t189 * t252 + t191 * t255) * t234 + (-t221 * t252 + t224 * t255) * t238;
t245 = Icges(2,4) * t256;
t229 = rSges(2,1) * t256 - rSges(2,2) * t253;
t228 = rSges(2,1) * t253 + rSges(2,2) * t256;
t227 = rSges(3,1) * t252 + rSges(3,2) * t255;
t226 = Icges(2,1) * t256 - t310;
t225 = Icges(2,1) * t253 + t245;
t223 = -Icges(2,2) * t253 + t245;
t222 = Icges(2,2) * t256 + t310;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t206 = rSges(4,1) * t242 + rSges(4,2) * t244;
t201 = t244 * t298 + t301;
t200 = -t244 * t300 + t299;
t199 = t244 * t299 - t300;
t198 = -t244 * t301 - t298;
t197 = (-qJD(4) - qJD(5)) * t244 + t238;
t194 = rSges(3,3) * t253 + t256 * t282;
t193 = -rSges(3,3) * t256 + t253 * t282;
t180 = rSges(4,3) * t253 + t256 * t281;
t179 = -rSges(4,3) * t256 + t253 * t281;
t169 = -rSges(5,3) * t244 + (rSges(5,1) * t254 - rSges(5,2) * t251) * t242;
t168 = -Icges(5,5) * t244 + (Icges(5,1) * t254 - Icges(5,4) * t251) * t242;
t167 = -Icges(5,6) * t244 + (Icges(5,4) * t254 - Icges(5,2) * t251) * t242;
t166 = -Icges(5,3) * t244 + (Icges(5,5) * t254 - Icges(5,6) * t251) * t242;
t165 = V_base(5) * rSges(2,3) - t228 * t238 + t288;
t164 = t229 * t238 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t228 * V_base(4) - t229 * V_base(5) + V_base(3);
t161 = -rSges(6,3) * t244 + (rSges(6,1) * t243 - rSges(6,2) * t241) * t242;
t152 = t256 * t290 + t182;
t151 = t253 * t290 + t181;
t146 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t304;
t145 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t305;
t144 = Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t304;
t143 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t305;
t142 = Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t304;
t141 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t305;
t140 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t304;
t139 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t305;
t135 = rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t304;
t133 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t305;
t118 = t227 * t234 + (-t193 - t232) * t238 + t288;
t117 = t194 * t238 - t227 * t235 + t274;
t113 = t193 * t235 - t194 * t234 + t273;
t112 = t206 * t208 + (-t179 + t294) * t238 + t284;
t111 = t180 * t238 - t206 * t209 + t266;
t110 = t179 * t209 - t180 * t208 + t267;
t109 = -t145 * t216 + t169 * t181 + t268;
t108 = t146 * t216 - t169 * t182 + t263;
t107 = t145 * t182 - t146 * t181 + t264;
t106 = -t133 * t197 + t151 * t161 + t265;
t105 = t135 * t197 - t152 * t161 + t261;
t104 = t133 * t152 - t135 * t151 + t262;
t103 = t151 * t295 - t197 * t297 + t256 * t289 + t265;
t102 = -t152 * t295 + t197 * t296 + t253 * t289 + t261;
t101 = -qJD(6) * t244 - t151 * t296 + t152 * t297 + t262;
t1 = t235 * (t270 * t253 + t259 * t256) / 0.2e1 + t234 * (t259 * t253 - t270 * t256) / 0.2e1 + t209 * (t271 * t253 + t260 * t256) / 0.2e1 + t208 * (t260 * t253 - t271 * t256) / 0.2e1 + m(2) * (t162 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(3) * (t113 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t182 * ((t140 * t304 + t200 * t142 + t201 * t144) * t182 + (t139 * t304 + t141 * t200 + t143 * t201) * t181 + (t166 * t304 + t167 * t200 + t168 * t201) * t216) / 0.2e1 + t181 * ((t140 * t305 + t142 * t198 + t144 * t199) * t182 + (t139 * t305 + t198 * t141 + t199 * t143) * t181 + (t166 * t305 + t167 * t198 + t168 * t199) * t216) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t216 * ((-t139 * t181 - t140 * t182 - t166 * t216) * t244 + ((-t142 * t251 + t144 * t254) * t182 + (-t141 * t251 + t143 * t254) * t181 + (-t167 * t251 + t168 * t254) * t216) * t242) / 0.2e1 + ((t183 * t322 + t184 * t321 + t305 * t323) * t197 + (t183 * t326 + t184 * t324 + t305 * t328) * t152 + (t327 * t183 + t325 * t184 + t329 * t305) * t151) * t151 / 0.2e1 + ((t185 * t322 + t186 * t321 + t304 * t323) * t197 + (t326 * t185 + t324 * t186 + t328 * t304) * t152 + (t327 * t185 + t325 * t186 + t304 * t329) * t151) * t152 / 0.2e1 + ((-t151 * t329 - t328 * t152 - t323 * t197) * t244 + ((-t241 * t322 + t243 * t321) * t197 + (-t241 * t326 + t243 * t324) * t152 + (-t241 * t327 + t243 * t325) * t151) * t242) * t197 / 0.2e1 + ((-t222 * t253 + t225 * t256 + Icges(1,4)) * V_base(5) + (-t253 * t223 + t256 * t226 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t256 * t222 + t253 * t225 + Icges(1,2)) * V_base(5) + (t223 * t256 + t226 * t253 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t190 * t255 + t192 * t252) * t235 + (t189 * t255 + t191 * t252) * t234 + (t176 * t244 + t178 * t242) * t209 + (t175 * t244 + t177 * t242) * t208 + (t244 * t204 + t242 * t205 + t255 * t221 + t252 * t224 + Icges(2,3)) * t238) * t238 / 0.2e1 + t238 * V_base(4) * (Icges(2,5) * t256 - Icges(2,6) * t253) + V_base(5) * t238 * (Icges(2,5) * t253 + Icges(2,6) * t256) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
