% Calculate kinetic energy for
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:29
% EndTime: 2019-03-09 16:43:33
% DurationCPUTime: 3.62s
% Computational Cost: add. (1668->289), mult. (2005->403), div. (0->0), fcn. (1853->8), ass. (0->154)
t344 = Icges(4,4) + Icges(5,6);
t343 = Icges(4,1) + Icges(5,2);
t342 = -Icges(4,2) - Icges(5,3);
t237 = qJ(2) + qJ(3);
t234 = cos(t237);
t341 = t344 * t234;
t233 = sin(t237);
t340 = t344 * t233;
t339 = Icges(5,4) - Icges(4,5);
t338 = Icges(5,5) - Icges(4,6);
t337 = t342 * t233 + t341;
t336 = t343 * t234 - t340;
t335 = Icges(5,1) + Icges(4,3);
t334 = Icges(6,1) + Icges(7,1);
t333 = Icges(6,4) - Icges(7,5);
t332 = Icges(7,4) + Icges(6,5);
t331 = Icges(6,2) + Icges(7,3);
t330 = Icges(7,6) - Icges(6,6);
t329 = Icges(6,3) + Icges(7,2);
t240 = sin(qJ(1));
t243 = cos(qJ(1));
t328 = t240 * t337 + t243 * t338;
t327 = -t240 * t338 + t243 * t337;
t326 = t336 * t240 + t243 * t339;
t325 = -t240 * t339 + t336 * t243;
t324 = t342 * t234 - t340;
t323 = t343 * t233 + t341;
t322 = t338 * t233 - t234 * t339;
t321 = rSges(7,1) + pkin(5);
t320 = rSges(7,3) + qJ(6);
t241 = cos(qJ(5));
t286 = t241 * t243;
t238 = sin(qJ(5));
t289 = t238 * t240;
t186 = -t233 * t286 + t289;
t287 = t240 * t241;
t288 = t238 * t243;
t187 = t233 * t288 + t287;
t290 = t234 * t243;
t319 = t331 * t186 - t333 * t187 + t330 * t290;
t188 = t233 * t287 + t288;
t189 = t233 * t289 - t286;
t291 = t234 * t240;
t318 = -t331 * t188 - t333 * t189 + t330 * t291;
t317 = t330 * t186 + t332 * t187 + t329 * t290;
t316 = -t330 * t188 + t332 * t189 + t329 * t291;
t315 = -t333 * t186 + t334 * t187 + t332 * t290;
t314 = t333 * t188 + t334 * t189 + t332 * t291;
t313 = (t333 * t238 + t331 * t241) * t234 + t330 * t233;
t312 = (-t332 * t238 + t330 * t241) * t234 + t329 * t233;
t311 = (-t334 * t238 - t333 * t241) * t234 + t332 * t233;
t202 = V_base(5) + (-qJD(2) - qJD(3)) * t243;
t228 = qJD(2) * t240 + V_base(4);
t203 = qJD(3) * t240 + t228;
t230 = V_base(6) + qJD(1);
t310 = (t324 * t233 + t323 * t234) * t230 + (-t327 * t233 + t325 * t234) * t203 + (-t328 * t233 + t326 * t234) * t202;
t309 = (-t233 * t339 - t338 * t234) * t230 + (t335 * t240 + t322 * t243) * t203 + (t322 * t240 - t335 * t243) * t202;
t239 = sin(qJ(2));
t302 = pkin(2) * t239;
t301 = pkin(9) * t233;
t242 = cos(qJ(2));
t300 = pkin(2) * t242;
t298 = Icges(2,4) * t240;
t297 = Icges(3,4) * t239;
t296 = Icges(3,4) * t242;
t285 = rSges(7,2) * t290 + t320 * t186 + t321 * t187;
t284 = rSges(7,2) * t291 - t320 * t188 + t321 * t189;
t283 = rSges(7,2) * t233 + (-t321 * t238 + t320 * t241) * t234;
t152 = -pkin(8) * t243 + t240 * t300;
t225 = t240 * pkin(1) - t243 * pkin(7);
t282 = -t152 - t225;
t281 = qJD(4) * t233;
t280 = qJD(4) * t234;
t279 = qJD(5) * t234;
t278 = V_base(5) * pkin(6) + V_base(1);
t270 = pkin(3) * t234 + qJ(4) * t233;
t183 = t270 * t240;
t275 = -t183 + t282;
t227 = -qJD(2) * t243 + V_base(5);
t274 = t227 * t302 + t278;
t273 = rSges(3,1) * t242 - rSges(3,2) * t239;
t272 = rSges(4,1) * t234 - rSges(4,2) * t233;
t271 = -rSges(5,2) * t234 + rSges(5,3) * t233;
t269 = Icges(3,1) * t242 - t297;
t267 = -Icges(3,2) * t239 + t296;
t264 = Icges(3,5) * t242 - Icges(3,6) * t239;
t199 = pkin(3) * t233 - qJ(4) * t234;
t260 = t202 * t199 + t243 * t281 + t274;
t226 = t243 * pkin(1) + t240 * pkin(7);
t259 = -V_base(4) * pkin(6) + t230 * t226 + V_base(2);
t258 = V_base(4) * t225 - t226 * V_base(5) + V_base(3);
t255 = (-Icges(3,3) * t243 + t240 * t264) * t227 + (Icges(3,3) * t240 + t243 * t264) * t228 + (Icges(3,5) * t239 + Icges(3,6) * t242) * t230;
t153 = pkin(8) * t240 + t243 * t300;
t254 = t228 * t152 - t153 * t227 + t258;
t253 = t230 * t153 - t228 * t302 + t259;
t252 = t203 * t183 + t254;
t185 = t270 * t243;
t251 = t230 * t185 + t240 * t281 + t253;
t192 = -pkin(4) * t243 + pkin(9) * t291;
t250 = t202 * t301 + (-t192 + t275) * t230 + t260;
t191 = pkin(4) * t240 + pkin(9) * t290;
t249 = t203 * t192 + (-t185 - t191) * t202 + t252;
t248 = t230 * t191 + (-t199 - t301) * t203 + t251;
t175 = -Icges(3,6) * t243 + t240 * t267;
t176 = Icges(3,6) * t240 + t243 * t267;
t177 = -Icges(3,5) * t243 + t240 * t269;
t178 = Icges(3,5) * t240 + t243 * t269;
t214 = Icges(3,2) * t242 + t297;
t217 = Icges(3,1) * t239 + t296;
t245 = (-t176 * t239 + t178 * t242) * t228 + (-t175 * t239 + t177 * t242) * t227 + (-t214 * t239 + t217 * t242) * t230;
t235 = Icges(2,4) * t243;
t222 = rSges(2,1) * t243 - rSges(2,2) * t240;
t221 = rSges(2,1) * t240 + rSges(2,2) * t243;
t220 = rSges(3,1) * t239 + rSges(3,2) * t242;
t219 = Icges(2,1) * t243 - t298;
t218 = Icges(2,1) * t240 + t235;
t216 = -Icges(2,2) * t240 + t235;
t215 = Icges(2,2) * t243 + t298;
t209 = qJD(5) * t233 + t230;
t208 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t207 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t206 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t201 = rSges(4,1) * t233 + rSges(4,2) * t234;
t200 = -rSges(5,2) * t233 - rSges(5,3) * t234;
t180 = rSges(3,3) * t240 + t243 * t273;
t179 = -rSges(3,3) * t243 + t240 * t273;
t172 = t243 * t279 + t203;
t171 = t240 * t279 + t202;
t170 = -rSges(5,1) * t243 + t240 * t271;
t169 = rSges(5,1) * t240 + t243 * t271;
t168 = rSges(4,3) * t240 + t243 * t272;
t167 = -rSges(4,3) * t243 + t240 * t272;
t151 = rSges(6,3) * t233 + (-rSges(6,1) * t238 - rSges(6,2) * t241) * t234;
t143 = V_base(5) * rSges(2,3) - t221 * t230 + t278;
t142 = t222 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t139 = t221 * V_base(4) - t222 * V_base(5) + V_base(3);
t133 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t291;
t131 = rSges(6,1) * t187 - rSges(6,2) * t186 + rSges(6,3) * t290;
t117 = t220 * t227 + (-t179 - t225) * t230 + t278;
t116 = t180 * t230 - t220 * t228 + t259;
t115 = t179 * t228 - t180 * t227 + t258;
t114 = t201 * t202 + (-t167 + t282) * t230 + t274;
t113 = t168 * t230 - t201 * t203 + t253;
t112 = t167 * t203 - t168 * t202 + t254;
t111 = t200 * t202 + (-t170 + t275) * t230 + t260;
t110 = t169 * t230 + (-t199 - t200) * t203 + t251;
t109 = -t280 + t170 * t203 + (-t169 - t185) * t202 + t252;
t108 = -t133 * t209 + t151 * t171 + t250;
t107 = t131 * t209 - t151 * t172 + t248;
t106 = -t131 * t171 + t133 * t172 + t249 - t280;
t105 = qJD(6) * t186 + t171 * t283 - t209 * t284 + t250;
t104 = -qJD(6) * t188 - t172 * t283 + t209 * t285 + t248;
t103 = (qJD(6) * t241 - qJD(4)) * t234 + t284 * t172 - t285 * t171 + t249;
t1 = t228 * (t240 * t255 + t243 * t245) / 0.2e1 + t227 * (t240 * t245 - t243 * t255) / 0.2e1 + m(1) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t139 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + ((-t188 * t313 + t189 * t311 + t291 * t312) * t209 + (-t188 * t319 + t315 * t189 + t317 * t291) * t172 + (-t318 * t188 + t314 * t189 + t316 * t291) * t171) * t171 / 0.2e1 + ((t186 * t313 + t187 * t311 + t290 * t312) * t209 + (t319 * t186 + t315 * t187 + t317 * t290) * t172 + (t186 * t318 + t187 * t314 + t290 * t316) * t171) * t172 / 0.2e1 + (t310 * t240 - t309 * t243) * t202 / 0.2e1 + (t309 * t240 + t310 * t243) * t203 / 0.2e1 + (((-t238 * t311 + t241 * t313) * t209 + (-t315 * t238 + t241 * t319) * t172 + (-t238 * t314 + t241 * t318) * t171) * t234 + (t171 * t316 + t172 * t317 + t209 * t312) * t233) * t209 / 0.2e1 + ((-t215 * t240 + t218 * t243 + Icges(1,4)) * V_base(5) + (-t216 * t240 + t219 * t243 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t215 * t243 + t218 * t240 + Icges(1,2)) * V_base(5) + (t216 * t243 + t219 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t176 * t242 + t178 * t239) * t228 + (t175 * t242 + t177 * t239) * t227 + (t325 * t233 + t327 * t234) * t203 + (t326 * t233 + t328 * t234) * t202 + (t214 * t242 + t217 * t239 + t323 * t233 - t324 * t234 + Icges(2,3)) * t230) * t230 / 0.2e1 + t230 * V_base(4) * (Icges(2,5) * t243 - Icges(2,6) * t240) + t230 * V_base(5) * (Icges(2,5) * t240 + Icges(2,6) * t243) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
