% Calculate kinetic energy for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:28:55
% EndTime: 2019-03-10 04:28:59
% DurationCPUTime: 3.96s
% Computational Cost: add. (3639->442), mult. (7098->660), div. (0->0), fcn. (8608->14), ass. (0->193)
t361 = cos(qJ(3));
t317 = cos(pkin(6));
t359 = pkin(8) * t317;
t322 = cos(qJ(4));
t358 = t322 * pkin(4);
t321 = sin(qJ(1));
t356 = Icges(2,4) * t321;
t323 = cos(qJ(2));
t324 = cos(qJ(1));
t347 = t323 * t324;
t320 = sin(qJ(2));
t350 = t320 * t321;
t274 = -t317 * t347 + t350;
t318 = sin(qJ(4));
t355 = t274 * t318;
t348 = t321 * t323;
t349 = t320 * t324;
t276 = t317 * t348 + t349;
t354 = t276 * t318;
t316 = sin(pkin(6));
t353 = t316 * t321;
t352 = t316 * t323;
t351 = t316 * t324;
t315 = qJ(4) + qJ(5);
t309 = cos(t315);
t346 = pkin(5) * t309;
t344 = qJD(2) * t316;
t343 = V_base(5) * pkin(7) + V_base(1);
t340 = t318 * t352;
t339 = t316 * t361;
t288 = t321 * t344 + V_base(4);
t307 = V_base(6) + qJD(1);
t308 = sin(t315);
t338 = pkin(5) * t308;
t253 = qJD(3) * t276 + t288;
t289 = qJD(2) * t317 + t307;
t277 = -t317 * t350 + t347;
t319 = sin(qJ(3));
t256 = t277 * t319 - t321 * t339;
t210 = qJD(4) * t256 + t253;
t287 = -t324 * t344 + V_base(5);
t186 = qJD(5) * t256 + t210;
t280 = pkin(1) * t321 - pkin(8) * t351;
t337 = -t280 * t307 + t359 * V_base(5) + t343;
t281 = pkin(1) * t324 + pkin(8) * t353;
t336 = t280 * V_base(4) - t281 * V_base(5) + V_base(3);
t252 = qJD(3) * t274 + t287;
t275 = t317 * t349 + t348;
t254 = t275 * t319 + t324 * t339;
t209 = qJD(4) * t254 + t252;
t270 = -qJD(3) * t352 + t289;
t335 = t307 * t281 + V_base(2) + (-pkin(7) - t359) * V_base(4);
t185 = qJD(5) * t254 + t209;
t272 = t316 * t319 * t320 - t317 * t361;
t239 = qJD(4) * t272 + t270;
t242 = pkin(2) * t275 + pkin(9) * t274;
t279 = (pkin(2) * t320 - pkin(9) * t323) * t316;
t334 = -t242 * t289 + t279 * t287 + t337;
t243 = pkin(2) * t277 + pkin(9) * t276;
t333 = t242 * t288 - t243 * t287 + t336;
t219 = qJD(5) * t272 + t239;
t332 = t243 * t289 - t279 * t288 + t335;
t255 = t275 * t361 - t319 * t351;
t207 = pkin(3) * t255 + pkin(10) * t254;
t273 = t317 * t319 + t320 * t339;
t238 = pkin(3) * t273 + pkin(10) * t272;
t331 = -t207 * t270 + t238 * t252 + t334;
t257 = t277 * t361 + t319 * t353;
t208 = pkin(3) * t257 + pkin(10) * t256;
t330 = t207 * t253 - t208 * t252 + t333;
t329 = t208 * t270 - t238 * t253 + t332;
t155 = pkin(4) * t355 + pkin(11) * t254 + t255 * t358;
t191 = -pkin(4) * t340 + pkin(11) * t272 + t273 * t358;
t328 = -t155 * t239 + t191 * t209 + t331;
t156 = pkin(4) * t354 + pkin(11) * t256 + t257 * t358;
t327 = t155 * t210 - t156 * t209 + t330;
t326 = t156 * t239 - t191 * t210 + t329;
t311 = qJ(6) + t315;
t310 = Icges(2,4) * t324;
t304 = cos(t311);
t303 = sin(t311);
t297 = rSges(2,1) * t324 - rSges(2,2) * t321;
t296 = rSges(2,1) * t321 + rSges(2,2) * t324;
t295 = Icges(2,1) * t324 - t356;
t294 = Icges(2,1) * t321 + t310;
t293 = -Icges(2,2) * t321 + t310;
t292 = Icges(2,2) * t324 + t356;
t285 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t284 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t283 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t265 = rSges(3,3) * t317 + (rSges(3,1) * t320 + rSges(3,2) * t323) * t316;
t264 = Icges(3,5) * t317 + (Icges(3,1) * t320 + Icges(3,4) * t323) * t316;
t263 = Icges(3,6) * t317 + (Icges(3,4) * t320 + Icges(3,2) * t323) * t316;
t262 = Icges(3,3) * t317 + (Icges(3,5) * t320 + Icges(3,6) * t323) * t316;
t261 = V_base(5) * rSges(2,3) - t296 * t307 + t343;
t260 = t297 * t307 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t258 = t296 * V_base(4) - t297 * V_base(5) + V_base(3);
t251 = t273 * t322 - t340;
t250 = -t273 * t318 - t322 * t352;
t245 = t273 * t309 - t308 * t352;
t244 = -t273 * t308 - t309 * t352;
t241 = t273 * t304 - t303 * t352;
t240 = -t273 * t303 - t304 * t352;
t237 = rSges(3,1) * t277 - rSges(3,2) * t276 + rSges(3,3) * t353;
t236 = rSges(3,1) * t275 - rSges(3,2) * t274 - rSges(3,3) * t351;
t235 = Icges(3,1) * t277 - Icges(3,4) * t276 + Icges(3,5) * t353;
t234 = Icges(3,1) * t275 - Icges(3,4) * t274 - Icges(3,5) * t351;
t233 = Icges(3,4) * t277 - Icges(3,2) * t276 + Icges(3,6) * t353;
t232 = Icges(3,4) * t275 - Icges(3,2) * t274 - Icges(3,6) * t351;
t231 = Icges(3,5) * t277 - Icges(3,6) * t276 + Icges(3,3) * t353;
t230 = Icges(3,5) * t275 - Icges(3,6) * t274 - Icges(3,3) * t351;
t229 = rSges(4,1) * t273 - rSges(4,2) * t272 - rSges(4,3) * t352;
t228 = Icges(4,1) * t273 - Icges(4,4) * t272 - Icges(4,5) * t352;
t227 = Icges(4,4) * t273 - Icges(4,2) * t272 - Icges(4,6) * t352;
t226 = Icges(4,5) * t273 - Icges(4,6) * t272 - Icges(4,3) * t352;
t223 = t257 * t322 + t354;
t222 = -t257 * t318 + t276 * t322;
t221 = t255 * t322 + t355;
t220 = -t255 * t318 + t274 * t322;
t218 = t257 * t309 + t276 * t308;
t217 = -t257 * t308 + t276 * t309;
t216 = t255 * t309 + t274 * t308;
t215 = -t255 * t308 + t274 * t309;
t214 = t257 * t304 + t276 * t303;
t213 = -t257 * t303 + t276 * t304;
t212 = t255 * t304 + t274 * t303;
t211 = -t255 * t303 + t274 * t304;
t205 = qJD(6) * t272 + t219;
t204 = rSges(4,1) * t257 - rSges(4,2) * t256 + rSges(4,3) * t276;
t203 = rSges(4,1) * t255 - rSges(4,2) * t254 + rSges(4,3) * t274;
t201 = Icges(4,1) * t257 - Icges(4,4) * t256 + Icges(4,5) * t276;
t200 = Icges(4,1) * t255 - Icges(4,4) * t254 + Icges(4,5) * t274;
t199 = Icges(4,4) * t257 - Icges(4,2) * t256 + Icges(4,6) * t276;
t198 = Icges(4,4) * t255 - Icges(4,2) * t254 + Icges(4,6) * t274;
t197 = Icges(4,5) * t257 - Icges(4,6) * t256 + Icges(4,3) * t276;
t196 = Icges(4,5) * t255 - Icges(4,6) * t254 + Icges(4,3) * t274;
t195 = rSges(5,1) * t251 + rSges(5,2) * t250 + rSges(5,3) * t272;
t194 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t272;
t193 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t272;
t192 = Icges(5,5) * t251 + Icges(5,6) * t250 + Icges(5,3) * t272;
t190 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t272;
t189 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t272;
t188 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t272;
t187 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t272;
t183 = rSges(7,1) * t241 + rSges(7,2) * t240 + rSges(7,3) * t272;
t182 = Icges(7,1) * t241 + Icges(7,4) * t240 + Icges(7,5) * t272;
t181 = Icges(7,4) * t241 + Icges(7,2) * t240 + Icges(7,6) * t272;
t180 = Icges(7,5) * t241 + Icges(7,6) * t240 + Icges(7,3) * t272;
t179 = -t236 * t289 + t265 * t287 + t337;
t178 = t237 * t289 - t265 * t288 + t335;
t177 = pkin(12) * t272 + t273 * t346 - t338 * t352;
t176 = qJD(6) * t256 + t186;
t175 = qJD(6) * t254 + t185;
t174 = rSges(5,1) * t223 + rSges(5,2) * t222 + rSges(5,3) * t256;
t173 = rSges(5,1) * t221 + rSges(5,2) * t220 + rSges(5,3) * t254;
t172 = Icges(5,1) * t223 + Icges(5,4) * t222 + Icges(5,5) * t256;
t171 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t254;
t170 = Icges(5,4) * t223 + Icges(5,2) * t222 + Icges(5,6) * t256;
t169 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t254;
t168 = Icges(5,5) * t223 + Icges(5,6) * t222 + Icges(5,3) * t256;
t167 = Icges(5,5) * t221 + Icges(5,6) * t220 + Icges(5,3) * t254;
t166 = t236 * t288 - t237 * t287 + t336;
t165 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t256;
t164 = rSges(6,1) * t216 + rSges(6,2) * t215 + rSges(6,3) * t254;
t163 = Icges(6,1) * t218 + Icges(6,4) * t217 + Icges(6,5) * t256;
t162 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t254;
t161 = Icges(6,4) * t218 + Icges(6,2) * t217 + Icges(6,6) * t256;
t160 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t254;
t159 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t256;
t158 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t254;
t154 = rSges(7,1) * t214 + rSges(7,2) * t213 + rSges(7,3) * t256;
t153 = rSges(7,1) * t212 + rSges(7,2) * t211 + rSges(7,3) * t254;
t152 = Icges(7,1) * t214 + Icges(7,4) * t213 + Icges(7,5) * t256;
t151 = Icges(7,1) * t212 + Icges(7,4) * t211 + Icges(7,5) * t254;
t150 = Icges(7,4) * t214 + Icges(7,2) * t213 + Icges(7,6) * t256;
t149 = Icges(7,4) * t212 + Icges(7,2) * t211 + Icges(7,6) * t254;
t148 = Icges(7,5) * t214 + Icges(7,6) * t213 + Icges(7,3) * t256;
t147 = Icges(7,5) * t212 + Icges(7,6) * t211 + Icges(7,3) * t254;
t145 = pkin(12) * t256 + t257 * t346 + t276 * t338;
t144 = pkin(12) * t254 + t255 * t346 + t274 * t338;
t142 = -t203 * t270 + t229 * t252 + t334;
t141 = t204 * t270 - t229 * t253 + t332;
t140 = t203 * t253 - t204 * t252 + t333;
t139 = -t173 * t239 + t195 * t209 + t331;
t138 = t174 * t239 - t195 * t210 + t329;
t137 = t173 * t210 - t174 * t209 + t330;
t136 = -t164 * t219 + t185 * t190 + t328;
t135 = t165 * t219 - t186 * t190 + t326;
t134 = t164 * t186 - t165 * t185 + t327;
t133 = -t144 * t219 - t153 * t205 + t175 * t183 + t177 * t185 + t328;
t132 = t145 * t219 + t154 * t205 - t176 * t183 - t177 * t186 + t326;
t131 = t144 * t186 - t145 * t185 + t153 * t176 - t154 * t175 + t327;
t1 = t287 * ((-t231 * t351 - t233 * t274 + t235 * t275) * t288 + (-t230 * t351 - t232 * t274 + t234 * t275) * t287 + (-t262 * t351 - t263 * t274 + t264 * t275) * t289) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((t292 * t324 + t294 * t321 + Icges(1,2)) * V_base(5) + (t293 * t324 + t295 * t321 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t292 * t321 + t294 * t324 + Icges(1,4)) * V_base(5) + (-t293 * t321 + t295 * t324 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((Icges(2,5) * t321 + Icges(2,6) * t324) * V_base(5) + (Icges(2,5) * t324 - Icges(2,6) * t321) * V_base(4) + Icges(2,3) * t307 / 0.2e1) * t307 + t289 * ((t230 * t287 + t231 * t288 + t262 * t289) * t317 + ((t233 * t323 + t235 * t320) * t288 + (t232 * t323 + t234 * t320) * t287 + (t263 * t323 + t264 * t320) * t289) * t316) / 0.2e1 + t270 * ((-t197 * t352 - t199 * t272 + t201 * t273) * t253 + (-t196 * t352 - t198 * t272 + t200 * t273) * t252 + (-t226 * t352 - t227 * t272 + t228 * t273) * t270) / 0.2e1 + t288 * ((t231 * t353 - t233 * t276 + t235 * t277) * t288 + (t230 * t353 - t232 * t276 + t234 * t277) * t287 + (t262 * t353 - t263 * t276 + t264 * t277) * t289) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(3) * (t166 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t175 * ((t148 * t254 + t150 * t211 + t152 * t212) * t176 + (t147 * t254 + t149 * t211 + t151 * t212) * t175 + (t180 * t254 + t181 * t211 + t182 * t212) * t205) / 0.2e1 + t185 * ((t159 * t254 + t161 * t215 + t163 * t216) * t186 + (t158 * t254 + t160 * t215 + t162 * t216) * t185 + (t187 * t254 + t188 * t215 + t189 * t216) * t219) / 0.2e1 + t209 * ((t168 * t254 + t170 * t220 + t172 * t221) * t210 + (t167 * t254 + t169 * t220 + t171 * t221) * t209 + (t192 * t254 + t193 * t220 + t194 * t221) * t239) / 0.2e1 + t176 * ((t148 * t256 + t150 * t213 + t152 * t214) * t176 + (t147 * t256 + t149 * t213 + t151 * t214) * t175 + (t180 * t256 + t181 * t213 + t182 * t214) * t205) / 0.2e1 + t186 * ((t159 * t256 + t161 * t217 + t163 * t218) * t186 + (t158 * t256 + t160 * t217 + t162 * t218) * t185 + (t187 * t256 + t188 * t217 + t189 * t218) * t219) / 0.2e1 + t210 * ((t168 * t256 + t170 * t222 + t172 * t223) * t210 + (t167 * t256 + t169 * t222 + t171 * t223) * t209 + (t192 * t256 + t193 * t222 + t194 * t223) * t239) / 0.2e1 + m(2) * (t258 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t205 * ((t148 * t272 + t150 * t240 + t152 * t241) * t176 + (t147 * t272 + t149 * t240 + t151 * t241) * t175 + (t180 * t272 + t181 * t240 + t182 * t241) * t205) / 0.2e1 + t219 * ((t159 * t272 + t161 * t244 + t163 * t245) * t186 + (t158 * t272 + t160 * t244 + t162 * t245) * t185 + (t187 * t272 + t188 * t244 + t189 * t245) * t219) / 0.2e1 + t239 * ((t168 * t272 + t170 * t250 + t172 * t251) * t210 + (t167 * t272 + t169 * t250 + t171 * t251) * t209 + (t192 * t272 + t193 * t250 + t194 * t251) * t239) / 0.2e1 + t252 * ((t197 * t274 - t199 * t254 + t201 * t255) * t253 + (t196 * t274 - t198 * t254 + t200 * t255) * t252 + (t226 * t274 - t227 * t254 + t228 * t255) * t270) / 0.2e1 + t253 * ((t197 * t276 - t199 * t256 + t201 * t257) * t253 + (t196 * t276 - t198 * t256 + t200 * t257) * t252 + (t226 * t276 - t227 * t256 + t228 * t257) * t270) / 0.2e1 + m(1) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1;
T  = t1;
