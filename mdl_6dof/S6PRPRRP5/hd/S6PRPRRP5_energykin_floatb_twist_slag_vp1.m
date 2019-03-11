% Calculate kinetic energy for
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:57
% EndTime: 2019-03-08 20:14:01
% DurationCPUTime: 3.96s
% Computational Cost: add. (2301->345), mult. (5348->481), div. (0->0), fcn. (6303->10), ass. (0->159)
t371 = Icges(3,1) + Icges(4,2);
t370 = Icges(6,1) + Icges(7,1);
t369 = Icges(3,4) + Icges(4,6);
t368 = Icges(6,4) + Icges(7,4);
t367 = Icges(3,5) - Icges(4,4);
t366 = Icges(6,5) + Icges(7,5);
t365 = Icges(3,2) + Icges(4,3);
t364 = Icges(6,2) + Icges(7,2);
t363 = Icges(3,6) - Icges(4,5);
t362 = Icges(6,6) + Icges(7,6);
t361 = Icges(3,3) + Icges(4,1);
t360 = Icges(6,3) + Icges(7,3);
t359 = rSges(7,3) + qJ(6);
t286 = sin(pkin(10));
t288 = cos(pkin(10));
t293 = sin(qJ(2));
t289 = cos(pkin(6));
t295 = cos(qJ(2));
t319 = t289 * t295;
t252 = t286 * t319 + t288 * t293;
t292 = sin(qJ(4));
t287 = sin(pkin(6));
t332 = cos(qJ(4));
t308 = t287 * t332;
t226 = t252 * t292 + t286 * t308;
t320 = t289 * t293;
t253 = -t286 * t320 + t288 * t295;
t291 = sin(qJ(5));
t294 = cos(qJ(5));
t190 = -t226 * t291 + t253 * t294;
t326 = t253 * t291;
t191 = t226 * t294 + t326;
t323 = t287 * t292;
t225 = -t252 * t332 + t286 * t323;
t357 = t190 * t362 + t191 * t366 + t225 * t360;
t250 = t286 * t293 - t288 * t319;
t228 = t250 * t292 - t288 * t308;
t251 = t286 * t295 + t288 * t320;
t192 = -t228 * t291 + t251 * t294;
t327 = t251 * t291;
t193 = t228 * t294 + t327;
t227 = t250 * t332 + t288 * t323;
t356 = t192 * t362 + t193 * t366 - t227 * t360;
t355 = t190 * t364 + t191 * t368 + t225 * t362;
t354 = t192 * t364 + t193 * t368 - t227 * t362;
t353 = t190 * t368 + t191 * t370 + t225 * t366;
t352 = t192 * t368 + t193 * t370 - t227 * t366;
t321 = t287 * t295;
t258 = t289 * t332 - t292 * t321;
t322 = t287 * t293;
t229 = -t258 * t291 + t294 * t322;
t309 = t291 * t322;
t230 = t258 * t294 + t309;
t257 = t289 * t292 + t295 * t308;
t351 = t229 * t362 + t230 * t366 + t257 * t360;
t350 = t229 * t364 + t230 * t368 + t257 * t362;
t349 = t229 * t368 + t230 * t370 + t257 * t366;
t325 = t286 * t287;
t348 = t252 * t365 - t253 * t369 - t325 * t363;
t324 = t287 * t288;
t347 = t250 * t365 - t251 * t369 + t324 * t363;
t346 = -t369 * t252 + t253 * t371 + t367 * t325;
t345 = -t369 * t250 + t251 * t371 - t367 * t324;
t344 = -t252 * t363 + t253 * t367 + t325 * t361;
t343 = -t250 * t363 + t251 * t367 - t324 * t361;
t342 = t361 * t289 + (t293 * t367 + t295 * t363) * t287;
t341 = t363 * t289 + (t293 * t369 + t295 * t365) * t287;
t340 = t367 * t289 + (t293 * t371 + t369 * t295) * t287;
t331 = pkin(7) * t289;
t330 = pkin(5) * t294;
t328 = Icges(2,4) * t286;
t318 = rSges(7,1) * t191 + rSges(7,2) * t190 + pkin(5) * t326 + t225 * t359 + t226 * t330;
t317 = rSges(7,1) * t193 + rSges(7,2) * t192 + pkin(5) * t327 - t227 * t359 + t228 * t330;
t316 = rSges(7,1) * t230 + rSges(7,2) * t229 + pkin(5) * t309 + t257 * t359 + t258 * t330;
t315 = qJD(2) * t287;
t314 = V_base(5) * qJ(1) + V_base(1);
t310 = qJD(1) + V_base(3);
t267 = t286 * t315 + V_base(4);
t278 = qJD(2) * t289 + V_base(6);
t224 = qJD(4) * t253 + t267;
t254 = qJD(4) * t322 + t278;
t266 = -t288 * t315 + V_base(5);
t223 = qJD(4) * t251 + t266;
t260 = pkin(1) * t286 - pkin(7) * t324;
t307 = -t260 * V_base(6) + V_base(5) * t331 + t314;
t261 = pkin(1) * t288 + pkin(7) * t325;
t306 = V_base(4) * t260 - V_base(5) * t261 + t310;
t305 = V_base(6) * t261 + V_base(2) + (-qJ(1) - t331) * V_base(4);
t259 = (pkin(2) * t293 - qJ(3) * t295) * t287;
t304 = qJD(3) * t252 + t266 * t259 + t307;
t218 = pkin(2) * t253 + qJ(3) * t252;
t303 = qJD(3) * t250 + t278 * t218 + t305;
t217 = pkin(2) * t251 + qJ(3) * t250;
t302 = -qJD(3) * t321 + t267 * t217 + t306;
t233 = -pkin(3) * t324 + pkin(8) * t251;
t262 = pkin(3) * t289 + pkin(8) * t322;
t301 = t266 * t262 + (-t217 - t233) * t278 + t304;
t232 = pkin(3) * t325 + pkin(8) * t253;
t300 = t278 * t232 + (-t259 - t262) * t267 + t303;
t299 = t267 * t233 + (-t218 - t232) * t266 + t302;
t187 = pkin(4) * t228 - pkin(9) * t227;
t219 = pkin(4) * t258 + pkin(9) * t257;
t298 = -t187 * t254 + t223 * t219 + t301;
t186 = pkin(4) * t226 + pkin(9) * t225;
t297 = t254 * t186 - t219 * t224 + t300;
t296 = -t223 * t186 + t224 * t187 + t299;
t284 = Icges(2,4) * t288;
t275 = rSges(2,1) * t288 - rSges(2,2) * t286;
t274 = rSges(2,1) * t286 + rSges(2,2) * t288;
t273 = Icges(2,1) * t288 - t328;
t272 = Icges(2,1) * t286 + t284;
t271 = -Icges(2,2) * t286 + t284;
t270 = Icges(2,2) * t288 + t328;
t265 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t264 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t263 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t244 = t289 * rSges(4,1) + (-rSges(4,2) * t293 - rSges(4,3) * t295) * t287;
t243 = t289 * rSges(3,3) + (rSges(3,1) * t293 + rSges(3,2) * t295) * t287;
t235 = V_base(5) * rSges(2,3) - t274 * V_base(6) + t314;
t234 = t275 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t222 = t274 * V_base(4) - t275 * V_base(5) + t310;
t221 = qJD(5) * t257 + t254;
t215 = rSges(5,1) * t258 - rSges(5,2) * t257 + rSges(5,3) * t322;
t214 = Icges(5,1) * t258 - Icges(5,4) * t257 + Icges(5,5) * t322;
t213 = Icges(5,4) * t258 - Icges(5,2) * t257 + Icges(5,6) * t322;
t212 = Icges(5,5) * t258 - Icges(5,6) * t257 + Icges(5,3) * t322;
t211 = rSges(3,1) * t253 - rSges(3,2) * t252 + rSges(3,3) * t325;
t210 = rSges(3,1) * t251 - rSges(3,2) * t250 - rSges(3,3) * t324;
t209 = -rSges(4,1) * t324 - rSges(4,2) * t251 + rSges(4,3) * t250;
t208 = rSges(4,1) * t325 - rSges(4,2) * t253 + rSges(4,3) * t252;
t189 = qJD(5) * t225 + t224;
t188 = -qJD(5) * t227 + t223;
t184 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t257;
t175 = rSges(5,1) * t228 + rSges(5,2) * t227 + rSges(5,3) * t251;
t174 = rSges(5,1) * t226 - rSges(5,2) * t225 + rSges(5,3) * t253;
t173 = Icges(5,1) * t228 + Icges(5,4) * t227 + Icges(5,5) * t251;
t172 = Icges(5,1) * t226 - Icges(5,4) * t225 + Icges(5,5) * t253;
t171 = Icges(5,4) * t228 + Icges(5,2) * t227 + Icges(5,6) * t251;
t170 = Icges(5,4) * t226 - Icges(5,2) * t225 + Icges(5,6) * t253;
t169 = Icges(5,5) * t228 + Icges(5,6) * t227 + Icges(5,3) * t251;
t168 = Icges(5,5) * t226 - Icges(5,6) * t225 + Icges(5,3) * t253;
t165 = -t210 * t278 + t243 * t266 + t307;
t164 = t211 * t278 - t243 * t267 + t305;
t163 = rSges(6,1) * t193 + rSges(6,2) * t192 - rSges(6,3) * t227;
t161 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t225;
t147 = t210 * t267 - t211 * t266 + t306;
t144 = t244 * t266 + (-t209 - t217) * t278 + t304;
t143 = t208 * t278 + (-t244 - t259) * t267 + t303;
t142 = t267 * t209 + (-t208 - t218) * t266 + t302;
t141 = -t175 * t254 + t215 * t223 + t301;
t140 = t174 * t254 - t215 * t224 + t300;
t139 = -t223 * t174 + t224 * t175 + t299;
t138 = -t163 * t221 + t184 * t188 + t298;
t137 = t161 * t221 - t184 * t189 + t297;
t136 = -t188 * t161 + t189 * t163 + t296;
t135 = qJD(6) * t225 + t188 * t316 - t221 * t317 + t298;
t134 = -qJD(6) * t227 - t189 * t316 + t221 * t318 + t297;
t133 = qJD(6) * t257 - t188 * t318 + t189 * t317 + t296;
t1 = m(3) * (t147 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + t224 * ((t253 * t168 - t225 * t170 + t226 * t172) * t224 + (t169 * t253 - t171 * t225 + t173 * t226) * t223 + (t212 * t253 - t213 * t225 + t214 * t226) * t254) / 0.2e1 + t223 * ((t168 * t251 + t170 * t227 + t172 * t228) * t224 + (t251 * t169 + t227 * t171 + t228 * t173) * t223 + (t212 * t251 + t213 * t227 + t214 * t228) * t254) / 0.2e1 + m(5) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(4) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(6) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(2) * (t222 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(7) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(1) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t254 * ((t168 * t322 - t170 * t257 + t172 * t258) * t224 + (t169 * t322 - t171 * t257 + t173 * t258) * t223 + (t212 * t322 - t257 * t213 + t258 * t214) * t254) / 0.2e1 + ((t192 * t350 + t193 * t349 - t227 * t351) * t221 + (t192 * t355 + t193 * t353 - t227 * t357) * t189 + (t354 * t192 + t352 * t193 - t356 * t227) * t188) * t188 / 0.2e1 + ((t190 * t350 + t191 * t349 + t225 * t351) * t221 + (t355 * t190 + t353 * t191 + t357 * t225) * t189 + (t190 * t354 + t191 * t352 + t225 * t356) * t188) * t189 / 0.2e1 + ((t350 * t229 + t349 * t230 + t351 * t257) * t221 + (t229 * t355 + t230 * t353 + t257 * t357) * t189 + (t229 * t354 + t230 * t352 + t257 * t356) * t188) * t221 / 0.2e1 + ((-t250 * t341 + t251 * t340 - t324 * t342) * t278 + (t250 * t348 + t251 * t346 - t324 * t344) * t267 + (t347 * t250 + t345 * t251 - t343 * t324) * t266) * t266 / 0.2e1 + ((-t252 * t341 + t253 * t340 + t325 * t342) * t278 + (t348 * t252 + t346 * t253 + t344 * t325) * t267 + (t252 * t347 + t253 * t345 + t325 * t343) * t266) * t267 / 0.2e1 + ((t266 * t343 + t267 * t344 + t278 * t342) * t289 + ((t293 * t340 + t295 * t341) * t278 + (t293 * t346 - t295 * t348) * t267 + (t293 * t345 - t347 * t295) * t266) * t287) * t278 / 0.2e1 + ((-t270 * t286 + t272 * t288 + Icges(1,4)) * V_base(5) + (-t271 * t286 + t273 * t288 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t270 * t288 + t272 * t286 + Icges(1,2)) * V_base(5) + (t271 * t288 + t273 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t288 - Icges(2,6) * t286 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t286 + Icges(2,6) * t288 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
