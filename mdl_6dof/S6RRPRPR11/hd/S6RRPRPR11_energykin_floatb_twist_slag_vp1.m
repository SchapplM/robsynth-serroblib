% Calculate kinetic energy for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:36
% EndTime: 2019-03-09 11:11:40
% DurationCPUTime: 5.00s
% Computational Cost: add. (1669->357), mult. (2292->510), div. (0->0), fcn. (2186->10), ass. (0->172)
t360 = Icges(3,4) + Icges(4,6);
t359 = Icges(3,1) + Icges(4,2);
t358 = -Icges(3,2) - Icges(4,3);
t274 = cos(qJ(2));
t357 = t360 * t274;
t271 = sin(qJ(2));
t356 = t360 * t271;
t355 = Icges(4,4) - Icges(3,5);
t354 = Icges(4,5) - Icges(3,6);
t353 = t358 * t271 + t357;
t352 = t359 * t274 - t356;
t351 = Icges(4,1) + Icges(3,3);
t350 = Icges(6,3) + Icges(5,3);
t272 = sin(qJ(1));
t275 = cos(qJ(1));
t349 = t272 * t353 + t275 * t354;
t348 = -t272 * t354 + t275 * t353;
t347 = t352 * t272 + t275 * t355;
t346 = -t272 * t355 + t352 * t275;
t345 = t358 * t274 - t356;
t344 = t359 * t271 + t357;
t343 = t354 * t271 - t274 * t355;
t268 = qJ(4) + pkin(10);
t258 = cos(t268);
t257 = sin(t268);
t319 = t272 * t257;
t322 = t271 * t275;
t180 = t258 * t322 - t319;
t318 = t272 * t258;
t181 = t257 * t322 + t318;
t273 = cos(qJ(4));
t314 = t273 * t275;
t270 = sin(qJ(4));
t317 = t272 * t270;
t209 = t271 * t314 - t317;
t316 = t272 * t273;
t210 = t270 * t322 + t316;
t313 = t274 * t275;
t342 = Icges(5,5) * t210 + Icges(6,5) * t181 + Icges(5,6) * t209 + Icges(6,6) * t180 + t350 * t313;
t182 = t257 * t275 + t271 * t318;
t183 = -t258 * t275 + t271 * t319;
t211 = t270 * t275 + t271 * t316;
t212 = t271 * t317 - t314;
t315 = t272 * t274;
t341 = Icges(5,5) * t212 + Icges(6,5) * t183 + Icges(5,6) * t211 + Icges(6,6) * t182 + t350 * t315;
t340 = (-Icges(5,5) * t270 - Icges(6,5) * t257 - Icges(5,6) * t273 - Icges(6,6) * t258) * t274 + t350 * t271;
t246 = -qJD(2) * t275 + V_base(5);
t247 = qJD(2) * t272 + V_base(4);
t260 = V_base(6) + qJD(1);
t339 = (t345 * t271 + t344 * t274) * t260 + (-t348 * t271 + t346 * t274) * t247 + (-t349 * t271 + t347 * t274) * t246;
t338 = (-t271 * t355 - t354 * t274) * t260 + (t351 * t272 + t343 * t275) * t247 + (t343 * t272 - t351 * t275) * t246;
t331 = pkin(4) * t270;
t330 = pkin(8) * t271;
t329 = t273 * pkin(4);
t327 = Icges(2,4) * t272;
t259 = qJ(6) + t268;
t254 = sin(t259);
t321 = t272 * t254;
t255 = cos(t259);
t320 = t272 * t255;
t298 = pkin(2) * t274 + qJ(3) * t271;
t213 = t298 * t272;
t244 = t272 * pkin(1) - pkin(7) * t275;
t312 = -t213 - t244;
t311 = pkin(5) * t258;
t309 = qJD(3) * t271;
t308 = qJD(4) * t274;
t307 = qJD(5) * t274;
t306 = qJD(6) * t274;
t305 = V_base(5) * pkin(6) + V_base(1);
t302 = pkin(5) * t257;
t208 = t275 * t308 + t247;
t238 = qJD(4) * t271 + t260;
t239 = pkin(2) * t271 - qJ(3) * t274;
t301 = t246 * t239 + t275 * t309 + t305;
t300 = rSges(3,1) * t274 - rSges(3,2) * t271;
t299 = -rSges(4,2) * t274 + rSges(4,3) * t271;
t207 = t272 * t308 + t246;
t245 = pkin(1) * t275 + t272 * pkin(7);
t291 = -V_base(4) * pkin(6) + t260 * t245 + V_base(2);
t290 = V_base(4) * t244 - t245 * V_base(5) + V_base(3);
t289 = qJ(5) * t274 + t271 * t331;
t214 = t298 * t275;
t286 = t260 * t214 + t272 * t309 + t291;
t285 = pkin(9) * t274 + t271 * t302;
t284 = -qJD(3) * t274 + t247 * t213 + t290;
t221 = -pkin(3) * t275 + pkin(8) * t315;
t283 = t246 * t330 + (-t221 + t312) * t260 + t301;
t205 = qJ(5) * t271 - t274 * t331;
t282 = t207 * t205 + t275 * t307 + t283;
t220 = t272 * pkin(3) + pkin(8) * t313;
t281 = t260 * t220 + (-t239 - t330) * t247 + t286;
t280 = t247 * t221 + (-t214 - t220) * t246 + t284;
t156 = t272 * t329 + t275 * t289;
t279 = t238 * t156 + t272 * t307 + t281;
t157 = t272 * t289 - t275 * t329;
t278 = qJD(5) * t271 + t208 * t157 + t280;
t264 = Icges(2,4) * t275;
t243 = rSges(2,1) * t275 - t272 * rSges(2,2);
t242 = t272 * rSges(2,1) + rSges(2,2) * t275;
t241 = rSges(3,1) * t271 + rSges(3,2) * t274;
t240 = -rSges(4,2) * t271 - rSges(4,3) * t274;
t237 = Icges(2,1) * t275 - t327;
t236 = Icges(2,1) * t272 + t264;
t234 = -Icges(2,2) * t272 + t264;
t233 = Icges(2,2) * t275 + t327;
t224 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t223 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t222 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t216 = qJD(6) * t271 + t238;
t203 = -rSges(4,1) * t275 + t272 * t299;
t202 = t272 * rSges(4,1) + t275 * t299;
t201 = t272 * rSges(3,3) + t275 * t300;
t200 = rSges(5,3) * t271 + (-rSges(5,1) * t270 - rSges(5,2) * t273) * t274;
t199 = -rSges(3,3) * t275 + t272 * t300;
t190 = Icges(5,5) * t271 + (-Icges(5,1) * t270 - Icges(5,4) * t273) * t274;
t187 = Icges(5,6) * t271 + (-Icges(5,4) * t270 - Icges(5,2) * t273) * t274;
t177 = -t255 * t275 + t271 * t321;
t176 = t254 * t275 + t271 * t320;
t175 = t254 * t322 + t320;
t174 = t255 * t322 - t321;
t173 = t275 * t306 + t208;
t172 = t272 * t306 + t207;
t170 = rSges(6,3) * t271 + (-rSges(6,1) * t257 - rSges(6,2) * t258) * t274;
t169 = Icges(6,5) * t271 + (-Icges(6,1) * t257 - Icges(6,4) * t258) * t274;
t168 = Icges(6,6) * t271 + (-Icges(6,4) * t257 - Icges(6,2) * t258) * t274;
t166 = rSges(7,3) * t271 + (-rSges(7,1) * t254 - rSges(7,2) * t255) * t274;
t165 = V_base(5) * rSges(2,3) - t242 * t260 + t305;
t164 = t243 * t260 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t163 = Icges(7,5) * t271 + (-Icges(7,1) * t254 - Icges(7,4) * t255) * t274;
t162 = Icges(7,6) * t271 + (-Icges(7,4) * t254 - Icges(7,2) * t255) * t274;
t161 = Icges(7,3) * t271 + (-Icges(7,5) * t254 - Icges(7,6) * t255) * t274;
t160 = t242 * V_base(4) - t243 * V_base(5) + V_base(3);
t158 = pkin(9) * t271 - t274 * t302;
t155 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t315;
t154 = t210 * rSges(5,1) + t209 * rSges(5,2) + rSges(5,3) * t313;
t153 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t315;
t152 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t313;
t151 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t315;
t150 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t313;
t146 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t315;
t145 = t181 * rSges(6,1) + t180 * rSges(6,2) + rSges(6,3) * t313;
t144 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t315;
t143 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t313;
t142 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t315;
t141 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t313;
t137 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t315;
t136 = t175 * rSges(7,1) + t174 * rSges(7,2) + rSges(7,3) * t313;
t135 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t315;
t134 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t313;
t133 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t315;
t132 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t313;
t131 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t315;
t130 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t313;
t129 = t241 * t246 + (-t199 - t244) * t260 + t305;
t128 = t201 * t260 - t241 * t247 + t291;
t127 = t272 * t285 - t275 * t311;
t126 = t272 * t311 + t275 * t285;
t125 = t199 * t247 - t201 * t246 + t290;
t124 = t240 * t246 + (-t203 + t312) * t260 + t301;
t123 = t202 * t260 + (-t239 - t240) * t247 + t286;
t122 = t203 * t247 + (-t202 - t214) * t246 + t284;
t121 = -t155 * t238 + t200 * t207 + t283;
t120 = t154 * t238 - t200 * t208 + t281;
t119 = -t154 * t207 + t155 * t208 + t280;
t118 = t170 * t207 + (-t146 - t157) * t238 + t282;
t117 = t145 * t238 + (-t170 - t205) * t208 + t279;
t116 = t146 * t208 + (-t145 - t156) * t207 + t278;
t115 = -t137 * t216 + t158 * t207 + t166 * t172 + (-t127 - t157) * t238 + t282;
t114 = t126 * t238 + t136 * t216 - t166 * t173 + (-t158 - t205) * t208 + t279;
t113 = t127 * t208 - t136 * t172 + t137 * t173 + (-t126 - t156) * t207 + t278;
t1 = t173 * ((t130 * t313 + t174 * t132 + t175 * t134) * t173 + (t131 * t313 + t174 * t133 + t175 * t135) * t172 + (t161 * t313 + t174 * t162 + t175 * t163) * t216) / 0.2e1 + t172 * ((t130 * t315 + t132 * t176 + t134 * t177) * t173 + (t131 * t315 + t133 * t176 + t135 * t177) * t172 + (t161 * t315 + t162 * t176 + t163 * t177) * t216) / 0.2e1 + t216 * ((t130 * t173 + t131 * t172 + t161 * t216) * t271 + ((-t132 * t255 - t134 * t254) * t173 + (-t133 * t255 - t135 * t254) * t172 + (-t162 * t255 - t163 * t254) * t216) * t274) / 0.2e1 + m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t160 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(3) * (t125 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(7) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + ((t168 * t182 + t169 * t183 + t187 * t211 + t190 * t212 + t315 * t340) * t238 + (t141 * t182 + t143 * t183 + t150 * t211 + t152 * t212 + t315 * t342) * t208 + (t142 * t182 + t144 * t183 + t151 * t211 + t153 * t212 + t341 * t315) * t207) * t207 / 0.2e1 + ((t180 * t168 + t181 * t169 + t209 * t187 + t210 * t190 + t313 * t340) * t238 + (t180 * t141 + t181 * t143 + t209 * t150 + t210 * t152 + t342 * t313) * t208 + (t180 * t142 + t181 * t144 + t209 * t151 + t210 * t153 + t313 * t341) * t207) * t208 / 0.2e1 + (((-t168 * t258 - t169 * t257 - t187 * t273 - t190 * t270) * t238 + (-t141 * t258 - t143 * t257 - t150 * t273 - t152 * t270) * t208 + (-t142 * t258 - t144 * t257 - t151 * t273 - t153 * t270) * t207) * t274 + (t341 * t207 + t208 * t342 + t340 * t238) * t271) * t238 / 0.2e1 + (t339 * t272 - t338 * t275) * t246 / 0.2e1 + (t338 * t272 + t339 * t275) * t247 / 0.2e1 + ((-t272 * t233 + t236 * t275 + Icges(1,4)) * V_base(5) + (-t272 * t234 + t237 * t275 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t233 * t275 + t272 * t236 + Icges(1,2)) * V_base(5) + (t234 * t275 + t272 * t237 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t346 * t271 + t348 * t274) * t247 + (t347 * t271 + t349 * t274) * t246 + (t344 * t271 - t345 * t274 + Icges(2,3)) * t260) * t260 / 0.2e1 + t260 * V_base(4) * (Icges(2,5) * t275 - Icges(2,6) * t272) + V_base(5) * t260 * (Icges(2,5) * t272 + Icges(2,6) * t275) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
