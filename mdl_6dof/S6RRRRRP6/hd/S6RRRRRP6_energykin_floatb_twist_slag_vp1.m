% Calculate kinetic energy for
% S6RRRRRP6
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:28
% EndTime: 2019-03-10 01:26:31
% DurationCPUTime: 4.13s
% Computational Cost: add. (2317->357), mult. (2828->537), div. (0->0), fcn. (2788->10), ass. (0->169)
t336 = Icges(6,1) + Icges(7,1);
t335 = -Icges(6,4) + Icges(7,5);
t334 = Icges(7,4) + Icges(6,5);
t333 = Icges(6,2) + Icges(7,3);
t332 = -Icges(7,6) + Icges(6,6);
t331 = -Icges(6,3) - Icges(7,2);
t330 = rSges(7,1) + pkin(5);
t329 = rSges(7,3) + qJ(6);
t263 = qJ(3) + qJ(4);
t259 = qJ(5) + t263;
t251 = sin(t259);
t252 = cos(t259);
t269 = cos(qJ(1));
t266 = sin(qJ(1));
t268 = cos(qJ(2));
t305 = t266 * t268;
t188 = t251 * t305 + t252 * t269;
t189 = -t251 * t269 + t252 * t305;
t265 = sin(qJ(2));
t307 = t265 * t266;
t328 = t333 * t188 + t335 * t189 - t332 * t307;
t304 = t268 * t269;
t190 = t251 * t304 - t266 * t252;
t191 = t251 * t266 + t252 * t304;
t306 = t265 * t269;
t327 = t333 * t190 + t335 * t191 - t332 * t306;
t326 = -t332 * t188 + t334 * t189 - t331 * t307;
t325 = -t332 * t190 + t334 * t191 - t331 * t306;
t324 = t335 * t188 + t336 * t189 + t334 * t307;
t323 = t335 * t190 + t336 * t191 + t334 * t306;
t322 = t332 * t268 + (t333 * t251 + t335 * t252) * t265;
t321 = t331 * t268 + (-t332 * t251 + t334 * t252) * t265;
t320 = -t334 * t268 + (t335 * t251 + t336 * t252) * t265;
t267 = cos(qJ(3));
t314 = t267 * pkin(3);
t312 = Icges(2,4) * t266;
t311 = Icges(3,4) * t265;
t310 = Icges(3,4) * t268;
t264 = sin(qJ(3));
t309 = t264 * t266;
t308 = t264 * t269;
t303 = rSges(7,2) * t307 + t329 * t188 + t189 * t330;
t302 = rSges(7,2) * t306 + t329 * t190 + t191 * t330;
t301 = -rSges(7,2) * t268 + (t329 * t251 + t252 * t330) * t265;
t257 = cos(t263);
t300 = pkin(4) * t257;
t298 = qJD(3) * t265;
t297 = qJD(4) * t265;
t296 = qJD(5) * t265;
t295 = -qJD(3) - qJD(4);
t294 = V_base(5) * pkin(6) + V_base(1);
t246 = qJD(2) * t266 + V_base(4);
t254 = V_base(6) + qJD(1);
t256 = sin(t263);
t291 = pkin(4) * t256;
t213 = t269 * t298 + t246;
t290 = pkin(2) * t268 + pkin(8) * t265;
t245 = -qJD(2) * t269 + V_base(5);
t289 = rSges(3,1) * t268 - rSges(3,2) * t265;
t187 = t269 * t297 + t213;
t288 = Icges(3,1) * t268 - t311;
t287 = -Icges(3,2) * t265 + t310;
t286 = Icges(3,5) * t268 - Icges(3,6) * t265;
t212 = t266 * t298 + t245;
t243 = pkin(1) * t269 + pkin(7) * t266;
t285 = -V_base(4) * pkin(6) + t254 * t243 + V_base(2);
t242 = pkin(1) * t266 - pkin(7) * t269;
t284 = V_base(4) * t242 - t243 * V_base(5) + V_base(3);
t186 = t266 * t297 + t212;
t283 = pkin(9) * t265 + t268 * t314;
t219 = t290 * t266;
t241 = t265 * pkin(2) - t268 * pkin(8);
t282 = t245 * t241 + (-t219 - t242) * t254 + t294;
t281 = (-Icges(3,3) * t269 + t266 * t286) * t245 + (Icges(3,3) * t266 + t269 * t286) * t246 + (Icges(3,5) * t265 + Icges(3,6) * t268) * t254;
t280 = pkin(10) * t265 + t268 * t300;
t220 = t290 * t269;
t279 = t254 * t220 - t241 * t246 + t285;
t278 = t246 * t219 - t220 * t245 + t284;
t163 = -pkin(3) * t308 + t266 * t283;
t180 = -pkin(9) * t268 + t265 * t314;
t237 = -qJD(3) * t268 + t254;
t277 = -t163 * t237 + t212 * t180 + t282;
t164 = pkin(3) * t309 + t269 * t283;
t276 = t237 * t164 - t180 * t213 + t279;
t275 = t213 * t163 - t164 * t212 + t278;
t122 = t266 * t280 - t269 * t291;
t166 = -pkin(10) * t268 + t265 * t300;
t221 = t268 * t295 + t254;
t274 = -t122 * t221 + t186 * t166 + t277;
t123 = t266 * t291 + t269 * t280;
t273 = t221 * t123 - t166 * t187 + t276;
t272 = t187 * t122 - t123 * t186 + t275;
t198 = -Icges(3,6) * t269 + t266 * t287;
t199 = Icges(3,6) * t266 + t269 * t287;
t201 = -Icges(3,5) * t269 + t266 * t288;
t202 = Icges(3,5) * t266 + t269 * t288;
t231 = Icges(3,2) * t268 + t311;
t234 = Icges(3,1) * t265 + t310;
t271 = (-t199 * t265 + t202 * t268) * t246 + (-t198 * t265 + t201 * t268) * t245 + (-t231 * t265 + t234 * t268) * t254;
t258 = Icges(2,4) * t269;
t240 = rSges(2,1) * t269 - rSges(2,2) * t266;
t239 = rSges(2,1) * t266 + rSges(2,2) * t269;
t238 = rSges(3,1) * t265 + rSges(3,2) * t268;
t236 = Icges(2,1) * t269 - t312;
t235 = Icges(2,1) * t266 + t258;
t233 = -Icges(2,2) * t266 + t258;
t232 = Icges(2,2) * t269 + t312;
t226 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t225 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t224 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t217 = t267 * t304 + t309;
t216 = -t264 * t304 + t266 * t267;
t215 = t267 * t305 - t308;
t214 = -t264 * t305 - t267 * t269;
t210 = (-qJD(5) + t295) * t268 + t254;
t209 = t256 * t266 + t257 * t304;
t208 = -t256 * t304 + t257 * t266;
t207 = -t256 * t269 + t257 * t305;
t206 = -t256 * t305 - t257 * t269;
t205 = rSges(3,3) * t266 + t269 * t289;
t204 = -rSges(3,3) * t269 + t266 * t289;
t203 = -rSges(4,3) * t268 + (rSges(4,1) * t267 - rSges(4,2) * t264) * t265;
t200 = -Icges(4,5) * t268 + (Icges(4,1) * t267 - Icges(4,4) * t264) * t265;
t197 = -Icges(4,6) * t268 + (Icges(4,4) * t267 - Icges(4,2) * t264) * t265;
t194 = -Icges(4,3) * t268 + (Icges(4,5) * t267 - Icges(4,6) * t264) * t265;
t185 = -rSges(5,3) * t268 + (rSges(5,1) * t257 - rSges(5,2) * t256) * t265;
t183 = -Icges(5,5) * t268 + (Icges(5,1) * t257 - Icges(5,4) * t256) * t265;
t182 = -Icges(5,6) * t268 + (Icges(5,4) * t257 - Icges(5,2) * t256) * t265;
t181 = -Icges(5,3) * t268 + (Icges(5,5) * t257 - Icges(5,6) * t256) * t265;
t179 = -rSges(6,3) * t268 + (rSges(6,1) * t252 - rSges(6,2) * t251) * t265;
t171 = V_base(5) * rSges(2,3) - t239 * t254 + t294;
t170 = t240 * t254 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t169 = t239 * V_base(4) - t240 * V_base(5) + V_base(3);
t168 = t269 * t296 + t187;
t167 = t266 * t296 + t186;
t162 = rSges(4,1) * t217 + rSges(4,2) * t216 + rSges(4,3) * t306;
t161 = rSges(4,1) * t215 + rSges(4,2) * t214 + rSges(4,3) * t307;
t160 = Icges(4,1) * t217 + Icges(4,4) * t216 + Icges(4,5) * t306;
t159 = Icges(4,1) * t215 + Icges(4,4) * t214 + Icges(4,5) * t307;
t158 = Icges(4,4) * t217 + Icges(4,2) * t216 + Icges(4,6) * t306;
t157 = Icges(4,4) * t215 + Icges(4,2) * t214 + Icges(4,6) * t307;
t156 = Icges(4,5) * t217 + Icges(4,6) * t216 + Icges(4,3) * t306;
t155 = Icges(4,5) * t215 + Icges(4,6) * t214 + Icges(4,3) * t307;
t152 = rSges(5,1) * t209 + rSges(5,2) * t208 + rSges(5,3) * t306;
t151 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t307;
t150 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t306;
t149 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t307;
t148 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t306;
t147 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t307;
t146 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t306;
t145 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t307;
t143 = rSges(6,1) * t191 - rSges(6,2) * t190 + rSges(6,3) * t306;
t141 = rSges(6,1) * t189 - rSges(6,2) * t188 + rSges(6,3) * t307;
t125 = t238 * t245 + (-t204 - t242) * t254 + t294;
t124 = t205 * t254 - t238 * t246 + t285;
t121 = t204 * t246 - t205 * t245 + t284;
t118 = -t161 * t237 + t203 * t212 + t282;
t117 = t162 * t237 - t203 * t213 + t279;
t116 = t161 * t213 - t162 * t212 + t278;
t115 = -t151 * t221 + t185 * t186 + t277;
t114 = t152 * t221 - t185 * t187 + t276;
t113 = t151 * t187 - t152 * t186 + t275;
t112 = -t141 * t210 + t167 * t179 + t274;
t111 = t143 * t210 - t168 * t179 + t273;
t110 = t141 * t168 - t143 * t167 + t272;
t109 = qJD(6) * t190 + t167 * t301 - t210 * t303 + t274;
t108 = qJD(6) * t188 - t168 * t301 + t210 * t302 + t273;
t107 = qJD(6) * t251 * t265 - t167 * t302 + t168 * t303 + t272;
t1 = m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t237 * ((-t155 * t212 - t156 * t213 - t194 * t237) * t268 + ((-t158 * t264 + t160 * t267) * t213 + (-t157 * t264 + t159 * t267) * t212 + (-t197 * t264 + t200 * t267) * t237) * t265) / 0.2e1 + t221 * ((-t145 * t186 - t146 * t187 - t181 * t221) * t268 + ((-t148 * t256 + t150 * t257) * t187 + (-t147 * t256 + t149 * t257) * t186 + (-t182 * t256 + t183 * t257) * t221) * t265) / 0.2e1 + m(1) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + t213 * ((t156 * t306 + t216 * t158 + t217 * t160) * t213 + (t155 * t306 + t157 * t216 + t159 * t217) * t212 + (t194 * t306 + t197 * t216 + t200 * t217) * t237) / 0.2e1 + t187 * ((t146 * t306 + t208 * t148 + t209 * t150) * t187 + (t145 * t306 + t147 * t208 + t149 * t209) * t186 + (t181 * t306 + t182 * t208 + t183 * t209) * t221) / 0.2e1 + t212 * ((t156 * t307 + t158 * t214 + t160 * t215) * t213 + (t155 * t307 + t214 * t157 + t215 * t159) * t212 + (t194 * t307 + t197 * t214 + t200 * t215) * t237) / 0.2e1 + t186 * ((t146 * t307 + t148 * t206 + t150 * t207) * t187 + (t145 * t307 + t206 * t147 + t207 * t149) * t186 + (t181 * t307 + t182 * t206 + t183 * t207) * t221) / 0.2e1 + t246 * (t266 * t281 + t269 * t271) / 0.2e1 + t245 * (t266 * t271 - t281 * t269) / 0.2e1 + m(2) * (t169 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + ((t188 * t322 + t189 * t320 + t307 * t321) * t210 + (t188 * t327 + t189 * t323 + t307 * t325) * t168 + (t328 * t188 + t324 * t189 + t326 * t307) * t167) * t167 / 0.2e1 + ((t190 * t322 + t191 * t320 + t306 * t321) * t210 + (t327 * t190 + t323 * t191 + t325 * t306) * t168 + (t190 * t328 + t324 * t191 + t326 * t306) * t167) * t168 / 0.2e1 + ((-t167 * t326 - t168 * t325 - t210 * t321) * t268 + ((t251 * t322 + t252 * t320) * t210 + (t251 * t327 + t252 * t323) * t168 + (t251 * t328 + t324 * t252) * t167) * t265) * t210 / 0.2e1 + ((t199 * t268 + t202 * t265) * t246 + (t198 * t268 + t201 * t265) * t245 + (t268 * t231 + t265 * t234 + Icges(2,3)) * t254) * t254 / 0.2e1 + ((-t232 * t266 + t235 * t269 + Icges(1,4)) * V_base(5) + (-t233 * t266 + t236 * t269 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t232 * t269 + t235 * t266 + Icges(1,2)) * V_base(5) + (t233 * t269 + t236 * t266 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t254 * (Icges(2,5) * t269 - Icges(2,6) * t266) + V_base(5) * t254 * (Icges(2,5) * t266 + Icges(2,6) * t269) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
