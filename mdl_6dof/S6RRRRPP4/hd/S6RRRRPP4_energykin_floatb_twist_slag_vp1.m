% Calculate kinetic energy for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:23
% EndTime: 2019-03-09 20:58:26
% DurationCPUTime: 3.50s
% Computational Cost: add. (2281->350), mult. (2768->504), div. (0->0), fcn. (2728->10), ass. (0->162)
t332 = Icges(6,1) + Icges(7,1);
t331 = -Icges(6,4) + Icges(7,5);
t330 = Icges(7,4) + Icges(6,5);
t329 = Icges(6,2) + Icges(7,3);
t328 = -Icges(7,6) + Icges(6,6);
t327 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t326 = rSges(7,1) + pkin(5);
t325 = rSges(7,3) + qJ(6);
t260 = qJ(3) + qJ(4);
t251 = pkin(10) + t260;
t248 = sin(t251);
t249 = cos(t251);
t266 = cos(qJ(1));
t263 = sin(qJ(1));
t265 = cos(qJ(2));
t301 = t263 * t265;
t186 = t248 * t301 + t249 * t266;
t187 = -t248 * t266 + t249 * t301;
t262 = sin(qJ(2));
t303 = t262 * t263;
t324 = t186 * t329 + t187 * t331 - t303 * t328;
t300 = t265 * t266;
t188 = t248 * t300 - t249 * t263;
t189 = t248 * t263 + t249 * t300;
t302 = t262 * t266;
t323 = t188 * t329 + t189 * t331 - t302 * t328;
t322 = t186 * t331 + t187 * t332 + t303 * t330;
t321 = t188 * t331 + t189 * t332 + t302 * t330;
t320 = t328 * t265 + (t248 * t329 + t249 * t331) * t262;
t319 = -t330 * t265 + (t248 * t331 + t249 * t332) * t262;
t254 = sin(t260);
t255 = cos(t260);
t204 = -t254 * t301 - t255 * t266;
t205 = -t254 * t266 + t255 * t301;
t318 = Icges(5,5) * t205 + Icges(5,6) * t204 - t186 * t328 + t187 * t330 - t303 * t327;
t206 = -t254 * t300 + t255 * t263;
t207 = t254 * t263 + t255 * t300;
t317 = Icges(5,5) * t207 + Icges(5,6) * t206 - t188 * t328 + t189 * t330 - t302 * t327;
t316 = t327 * t265 + (Icges(5,5) * t255 - Icges(5,6) * t254 - t248 * t328 + t249 * t330) * t262;
t264 = cos(qJ(3));
t310 = t264 * pkin(3);
t308 = Icges(2,4) * t263;
t307 = Icges(3,4) * t262;
t306 = Icges(3,4) * t265;
t261 = sin(qJ(3));
t305 = t261 * t263;
t304 = t261 * t266;
t299 = rSges(7,2) * t303 + t186 * t325 + t187 * t326;
t298 = rSges(7,2) * t302 + t188 * t325 + t189 * t326;
t297 = -rSges(7,2) * t265 + (t248 * t325 + t249 * t326) * t262;
t296 = pkin(4) * t255;
t294 = qJD(3) * t262;
t293 = qJD(4) * t262;
t292 = qJD(5) * t262;
t291 = V_base(5) * pkin(6) + V_base(1);
t243 = qJD(2) * t263 + V_base(4);
t252 = V_base(6) + qJD(1);
t288 = pkin(4) * t254;
t210 = t266 * t294 + t243;
t287 = pkin(2) * t265 + pkin(8) * t262;
t242 = -qJD(2) * t266 + V_base(5);
t286 = rSges(3,1) * t265 - rSges(3,2) * t262;
t285 = Icges(3,1) * t265 - t307;
t284 = -Icges(3,2) * t262 + t306;
t283 = Icges(3,5) * t265 - Icges(3,6) * t262;
t209 = t263 * t294 + t242;
t241 = pkin(1) * t266 + pkin(7) * t263;
t282 = -V_base(4) * pkin(6) + t241 * t252 + V_base(2);
t240 = pkin(1) * t263 - pkin(7) * t266;
t281 = t240 * V_base(4) - t241 * V_base(5) + V_base(3);
t280 = pkin(9) * t262 + t265 * t310;
t216 = t287 * t263;
t239 = t262 * pkin(2) - pkin(8) * t265;
t279 = t242 * t239 + (-t216 - t240) * t252 + t291;
t278 = (-Icges(3,3) * t266 + t263 * t283) * t242 + (Icges(3,3) * t263 + t266 * t283) * t243 + (Icges(3,5) * t262 + Icges(3,6) * t265) * t252;
t277 = qJ(5) * t262 + t265 * t296;
t217 = t287 * t266;
t276 = t217 * t252 - t239 * t243 + t282;
t275 = t216 * t243 - t217 * t242 + t281;
t163 = -pkin(3) * t304 + t263 * t280;
t178 = -pkin(9) * t265 + t262 * t310;
t234 = -qJD(3) * t265 + t252;
t274 = -t163 * t234 + t178 * t209 + t279;
t164 = pkin(3) * t305 + t266 * t280;
t273 = t164 * t234 - t178 * t210 + t276;
t166 = -qJ(5) * t265 + t262 * t296;
t184 = t263 * t293 + t209;
t272 = t166 * t184 + t266 * t292 + t274;
t271 = t163 * t210 - t164 * t209 + t275;
t123 = t263 * t288 + t266 * t277;
t218 = (-qJD(3) - qJD(4)) * t265 + t252;
t270 = t123 * t218 + t263 * t292 + t273;
t122 = t263 * t277 - t266 * t288;
t185 = t266 * t293 + t210;
t269 = -qJD(5) * t265 + t122 * t185 + t271;
t196 = -Icges(3,6) * t266 + t263 * t284;
t197 = Icges(3,6) * t263 + t266 * t284;
t199 = -Icges(3,5) * t266 + t263 * t285;
t200 = Icges(3,5) * t263 + t266 * t285;
t228 = Icges(3,2) * t265 + t307;
t231 = Icges(3,1) * t262 + t306;
t268 = (-t197 * t262 + t200 * t265) * t243 + (-t196 * t262 + t199 * t265) * t242 + (-t228 * t262 + t231 * t265) * t252;
t256 = Icges(2,4) * t266;
t237 = rSges(2,1) * t266 - rSges(2,2) * t263;
t236 = rSges(2,1) * t263 + rSges(2,2) * t266;
t235 = rSges(3,1) * t262 + rSges(3,2) * t265;
t233 = Icges(2,1) * t266 - t308;
t232 = Icges(2,1) * t263 + t256;
t230 = -Icges(2,2) * t263 + t256;
t229 = Icges(2,2) * t266 + t308;
t223 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t222 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t221 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t214 = t264 * t300 + t305;
t213 = -t261 * t300 + t263 * t264;
t212 = t264 * t301 - t304;
t211 = -t261 * t301 - t264 * t266;
t203 = rSges(3,3) * t263 + t266 * t286;
t202 = -rSges(3,3) * t266 + t263 * t286;
t201 = -rSges(4,3) * t265 + (rSges(4,1) * t264 - rSges(4,2) * t261) * t262;
t198 = -Icges(4,5) * t265 + (Icges(4,1) * t264 - Icges(4,4) * t261) * t262;
t195 = -Icges(4,6) * t265 + (Icges(4,4) * t264 - Icges(4,2) * t261) * t262;
t192 = -Icges(4,3) * t265 + (Icges(4,5) * t264 - Icges(4,6) * t261) * t262;
t183 = -rSges(5,3) * t265 + (rSges(5,1) * t255 - rSges(5,2) * t254) * t262;
t181 = -Icges(5,5) * t265 + (Icges(5,1) * t255 - Icges(5,4) * t254) * t262;
t180 = -Icges(5,6) * t265 + (Icges(5,4) * t255 - Icges(5,2) * t254) * t262;
t177 = -rSges(6,3) * t265 + (rSges(6,1) * t249 - rSges(6,2) * t248) * t262;
t175 = V_base(5) * rSges(2,3) - t236 * t252 + t291;
t174 = t237 * t252 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t162 = rSges(4,1) * t214 + rSges(4,2) * t213 + rSges(4,3) * t302;
t161 = rSges(4,1) * t212 + rSges(4,2) * t211 + rSges(4,3) * t303;
t160 = Icges(4,1) * t214 + Icges(4,4) * t213 + Icges(4,5) * t302;
t159 = Icges(4,1) * t212 + Icges(4,4) * t211 + Icges(4,5) * t303;
t158 = Icges(4,4) * t214 + Icges(4,2) * t213 + Icges(4,6) * t302;
t157 = Icges(4,4) * t212 + Icges(4,2) * t211 + Icges(4,6) * t303;
t156 = Icges(4,5) * t214 + Icges(4,6) * t213 + Icges(4,3) * t302;
t155 = Icges(4,5) * t212 + Icges(4,6) * t211 + Icges(4,3) * t303;
t152 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t302;
t151 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t303;
t150 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t302;
t149 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t303;
t148 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t302;
t147 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t303;
t143 = rSges(6,1) * t189 - rSges(6,2) * t188 + rSges(6,3) * t302;
t141 = rSges(6,1) * t187 - rSges(6,2) * t186 + rSges(6,3) * t303;
t125 = t235 * t242 + (-t202 - t240) * t252 + t291;
t124 = t203 * t252 - t235 * t243 + t282;
t121 = t202 * t243 - t203 * t242 + t281;
t118 = -t161 * t234 + t201 * t209 + t279;
t117 = t162 * t234 - t201 * t210 + t276;
t116 = t161 * t210 - t162 * t209 + t275;
t115 = -t151 * t218 + t183 * t184 + t274;
t114 = t152 * t218 - t183 * t185 + t273;
t113 = t151 * t185 - t152 * t184 + t271;
t112 = t177 * t184 + (-t122 - t141) * t218 + t272;
t111 = t143 * t218 + (-t166 - t177) * t185 + t270;
t110 = t141 * t185 + (-t123 - t143) * t184 + t269;
t109 = qJD(6) * t188 + t297 * t184 + (-t122 - t299) * t218 + t272;
t108 = qJD(6) * t186 + t298 * t218 + (-t166 - t297) * t185 + t270;
t107 = qJD(6) * t248 * t262 + t299 * t185 + (-t123 - t298) * t184 + t269;
t1 = t243 * (t263 * t278 + t266 * t268) / 0.2e1 + t242 * (t263 * t268 - t278 * t266) / 0.2e1 + t210 * ((t156 * t302 + t213 * t158 + t214 * t160) * t210 + (t155 * t302 + t157 * t213 + t159 * t214) * t209 + (t192 * t302 + t195 * t213 + t198 * t214) * t234) / 0.2e1 + t209 * ((t156 * t303 + t158 * t211 + t160 * t212) * t210 + (t155 * t303 + t211 * t157 + t212 * t159) * t209 + (t192 * t303 + t195 * t211 + t198 * t212) * t234) / 0.2e1 + m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t234 * ((-t155 * t209 - t156 * t210 - t192 * t234) * t265 + ((-t158 * t261 + t160 * t264) * t210 + (-t157 * t261 + t159 * t264) * t209 + (-t195 * t261 + t198 * t264) * t234) * t262) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + ((t197 * t265 + t200 * t262) * t243 + (t196 * t265 + t199 * t262) * t242 + (t265 * t228 + t262 * t231 + Icges(2,3)) * t252) * t252 / 0.2e1 + ((-t229 * t263 + t232 * t266 + Icges(1,4)) * V_base(5) + (-t230 * t263 + t233 * t266 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t229 * t266 + t232 * t263 + Icges(1,2)) * V_base(5) + (t230 * t266 + t233 * t263 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t180 * t204 + t181 * t205 + t186 * t320 + t187 * t319 + t303 * t316) * t218 + (t148 * t204 + t150 * t205 + t186 * t323 + t187 * t321 + t303 * t317) * t185 + (t204 * t147 + t205 * t149 + t324 * t186 + t322 * t187 + t318 * t303) * t184) * t184 / 0.2e1 + ((t180 * t206 + t181 * t207 + t188 * t320 + t189 * t319 + t302 * t316) * t218 + (t206 * t148 + t207 * t150 + t323 * t188 + t321 * t189 + t317 * t302) * t185 + (t147 * t206 + t149 * t207 + t188 * t324 + t189 * t322 + t302 * t318) * t184) * t185 / 0.2e1 + ((-t184 * t318 - t185 * t317 - t218 * t316) * t265 + ((-t180 * t254 + t181 * t255 + t248 * t320 + t249 * t319) * t218 + (-t148 * t254 + t150 * t255 + t248 * t323 + t249 * t321) * t185 + (-t147 * t254 + t149 * t255 + t248 * t324 + t249 * t322) * t184) * t262) * t218 / 0.2e1 + V_base(4) * t252 * (Icges(2,5) * t266 - Icges(2,6) * t263) + V_base(5) * t252 * (Icges(2,5) * t263 + Icges(2,6) * t266) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
