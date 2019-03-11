% Calculate kinetic energy for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:40
% EndTime: 2019-03-09 12:16:43
% DurationCPUTime: 3.98s
% Computational Cost: add. (1478->300), mult. (3068->421), div. (0->0), fcn. (3327->8), ass. (0->147)
t344 = Icges(3,4) - Icges(4,5);
t343 = Icges(3,1) + Icges(4,1);
t342 = Icges(3,2) + Icges(4,3);
t258 = cos(qJ(2));
t341 = t344 * t258;
t255 = sin(qJ(2));
t340 = t344 * t255;
t339 = Icges(4,4) + Icges(3,5);
t338 = Icges(3,6) - Icges(4,6);
t337 = t342 * t255 - t341;
t336 = t343 * t258 - t340;
t335 = Icges(6,1) + Icges(7,1);
t334 = -Icges(6,4) + Icges(7,5);
t333 = Icges(7,4) + Icges(6,5);
t332 = Icges(4,2) + Icges(3,3);
t331 = Icges(6,2) + Icges(7,3);
t330 = Icges(7,6) - Icges(6,6);
t329 = Icges(6,3) + Icges(7,2);
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t328 = t337 * t256 + t338 * t259;
t327 = -t338 * t256 + t337 * t259;
t326 = t336 * t256 - t339 * t259;
t325 = t339 * t256 + t336 * t259;
t324 = -t342 * t258 - t340;
t323 = t343 * t255 + t341;
t322 = -t338 * t255 + t339 * t258;
t321 = rSges(7,1) + pkin(5);
t320 = rSges(7,3) + qJ(6);
t254 = sin(qJ(4));
t302 = cos(qJ(4));
t213 = t255 * t254 + t258 * t302;
t204 = t213 * t256;
t253 = sin(qJ(5));
t257 = cos(qJ(5));
t176 = t204 * t253 - t259 * t257;
t177 = t204 * t257 + t253 * t259;
t284 = t255 * t302;
t214 = -t258 * t254 + t284;
t203 = t214 * t256;
t319 = t176 * t331 + t177 * t334 - t203 * t330;
t206 = t213 * t259;
t178 = t206 * t253 + t256 * t257;
t179 = t206 * t257 - t253 * t256;
t294 = t258 * t259;
t205 = t254 * t294 - t259 * t284;
t318 = t178 * t331 + t179 * t334 + t205 * t330;
t317 = t176 * t330 + t177 * t333 - t203 * t329;
t316 = t178 * t330 + t179 * t333 + t205 * t329;
t315 = t334 * t176 + t177 * t335 - t333 * t203;
t314 = t334 * t178 + t179 * t335 + t333 * t205;
t313 = (t253 * t331 + t257 * t334) * t214 + t330 * t213;
t312 = (t253 * t330 + t257 * t333) * t214 + t329 * t213;
t311 = (t334 * t253 + t257 * t335) * t214 + t333 * t213;
t243 = -qJD(2) * t259 + V_base(5);
t244 = qJD(2) * t256 + V_base(4);
t248 = V_base(6) + qJD(1);
t310 = (t255 * t324 + t258 * t323) * t248 + (t255 * t327 + t258 * t325) * t244 + (t255 * t328 + t258 * t326) * t243;
t309 = (t339 * t255 + t338 * t258) * t248 + (t256 * t332 + t259 * t322) * t244 + (t256 * t322 - t259 * t332) * t243;
t301 = pkin(3) * t255;
t300 = Icges(2,4) * t256;
t293 = -rSges(7,2) * t203 + t320 * t176 + t321 * t177;
t292 = rSges(7,2) * t205 + t320 * t178 + t321 * t179;
t291 = rSges(7,2) * t213 + (t320 * t253 + t321 * t257) * t214;
t280 = pkin(2) * t258 + qJ(3) * t255;
t208 = t280 * t256;
t239 = t256 * pkin(1) - pkin(7) * t259;
t290 = -t208 - t239;
t289 = qJD(3) * t255;
t288 = V_base(5) * pkin(6) + V_base(1);
t217 = t256 * t258 * pkin(3) + pkin(8) * t259;
t285 = -t217 + t290;
t234 = pkin(2) * t255 - qJ(3) * t258;
t283 = t243 * t234 + t259 * t289 + t288;
t282 = rSges(3,1) * t258 - rSges(3,2) * t255;
t281 = rSges(4,1) * t258 + rSges(4,3) * t255;
t273 = t243 * t301 + t283;
t211 = qJD(4) * t259 + t243;
t212 = -qJD(4) * t256 + t244;
t240 = pkin(1) * t259 + t256 * pkin(7);
t272 = -V_base(4) * pkin(6) + t248 * t240 + V_base(2);
t271 = V_base(4) * t239 - t240 * V_base(5) + V_base(3);
t209 = t280 * t259;
t268 = t248 * t209 + t256 * t289 + t272;
t267 = -qJD(3) * t258 + t244 * t208 + t271;
t165 = pkin(4) * t204 - pkin(9) * t203;
t174 = pkin(4) * t214 + pkin(9) * t213;
t266 = t211 * t174 + (-t165 + t285) * t248 + t273;
t218 = pkin(3) * t294 - t256 * pkin(8);
t265 = t248 * t218 + (-t234 - t301) * t244 + t268;
t264 = t244 * t217 + (-t209 - t218) * t243 + t267;
t166 = pkin(4) * t206 + pkin(9) * t205;
t263 = t248 * t166 - t174 * t212 + t265;
t262 = t212 * t165 - t166 * t211 + t264;
t251 = Icges(2,4) * t259;
t238 = rSges(2,1) * t259 - t256 * rSges(2,2);
t237 = t256 * rSges(2,1) + rSges(2,2) * t259;
t236 = rSges(3,1) * t255 + rSges(3,2) * t258;
t235 = rSges(4,1) * t255 - rSges(4,3) * t258;
t233 = Icges(2,1) * t259 - t300;
t232 = Icges(2,1) * t256 + t251;
t229 = -Icges(2,2) * t256 + t251;
t228 = Icges(2,2) * t259 + t300;
t221 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t220 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t219 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t201 = t256 * rSges(3,3) + t259 * t282;
t200 = t256 * rSges(4,2) + t259 * t281;
t199 = -rSges(3,3) * t259 + t256 * t282;
t198 = -rSges(4,2) * t259 + t256 * t281;
t197 = qJD(5) * t213 + t248;
t181 = V_base(5) * rSges(2,3) - t237 * t248 + t288;
t180 = t238 * t248 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t175 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t173 = rSges(5,1) * t214 - rSges(5,2) * t213;
t172 = Icges(5,1) * t214 - Icges(5,4) * t213;
t171 = Icges(5,4) * t214 - Icges(5,2) * t213;
t170 = Icges(5,5) * t214 - Icges(5,6) * t213;
t169 = qJD(5) * t205 + t212;
t168 = -qJD(5) * t203 + t211;
t163 = rSges(5,1) * t206 - rSges(5,2) * t205 - rSges(5,3) * t256;
t162 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t259;
t161 = Icges(5,1) * t206 - Icges(5,4) * t205 - Icges(5,5) * t256;
t160 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t259;
t159 = Icges(5,4) * t206 - Icges(5,2) * t205 - Icges(5,6) * t256;
t158 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t259;
t157 = Icges(5,5) * t206 - Icges(5,6) * t205 - Icges(5,3) * t256;
t156 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t259;
t154 = rSges(6,3) * t213 + (rSges(6,1) * t257 - rSges(6,2) * t253) * t214;
t143 = t236 * t243 + (-t199 - t239) * t248 + t288;
t142 = t201 * t248 - t236 * t244 + t272;
t141 = rSges(6,1) * t179 - rSges(6,2) * t178 + rSges(6,3) * t205;
t139 = rSges(6,1) * t177 - rSges(6,2) * t176 - rSges(6,3) * t203;
t125 = t199 * t244 - t201 * t243 + t271;
t124 = t235 * t243 + (-t198 + t290) * t248 + t283;
t123 = t200 * t248 + (-t234 - t235) * t244 + t268;
t122 = t198 * t244 + (-t200 - t209) * t243 + t267;
t121 = t173 * t211 + (-t162 + t285) * t248 + t273;
t120 = t163 * t248 - t173 * t212 + t265;
t119 = t162 * t212 - t163 * t211 + t264;
t118 = -t139 * t197 + t154 * t168 + t266;
t117 = t141 * t197 - t154 * t169 + t263;
t116 = t139 * t169 - t141 * t168 + t262;
t115 = qJD(6) * t178 + t168 * t291 - t197 * t293 + t266;
t114 = qJD(6) * t176 - t169 * t291 + t197 * t292 + t263;
t113 = qJD(6) * t214 * t253 - t168 * t292 + t169 * t293 + t262;
t1 = m(7) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t125 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(2) * (t175 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + t212 * ((-t256 * t157 - t205 * t159 + t206 * t161) * t212 + (-t156 * t256 - t158 * t205 + t160 * t206) * t211 + (-t170 * t256 - t171 * t205 + t172 * t206) * t248) / 0.2e1 + t211 * ((t157 * t259 + t203 * t159 + t204 * t161) * t212 + (t259 * t156 + t203 * t158 + t204 * t160) * t211 + (t170 * t259 + t203 * t171 + t204 * t172) * t248) / 0.2e1 + ((t176 * t313 + t177 * t311 - t203 * t312) * t197 + (t176 * t318 + t177 * t314 - t203 * t316) * t169 + (t319 * t176 + t315 * t177 - t317 * t203) * t168) * t168 / 0.2e1 + ((t178 * t313 + t179 * t311 + t205 * t312) * t197 + (t318 * t178 + t314 * t179 + t316 * t205) * t169 + (t178 * t319 + t315 * t179 + t317 * t205) * t168) * t169 / 0.2e1 + (((t253 * t313 + t257 * t311) * t197 + (t253 * t318 + t257 * t314) * t169 + (t253 * t319 + t315 * t257) * t168) * t214 + (t168 * t317 + t169 * t316 + t197 * t312) * t213) * t197 / 0.2e1 + (t310 * t256 - t309 * t259) * t243 / 0.2e1 + (t309 * t256 + t310 * t259) * t244 / 0.2e1 + ((-t256 * t228 + t232 * t259 + Icges(1,4)) * V_base(5) + (-t256 * t229 + t259 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t259 * t228 + t256 * t232 + Icges(1,2)) * V_base(5) + (t229 * t259 + t256 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t159 * t213 + t161 * t214) * t212 + (-t158 * t213 + t160 * t214) * t211 + (t255 * t325 - t258 * t327) * t244 + (t255 * t326 - t258 * t328) * t243 + (-t213 * t171 + t214 * t172 + t255 * t323 - t258 * t324 + Icges(2,3)) * t248) * t248 / 0.2e1 + t248 * V_base(4) * (Icges(2,5) * t259 - Icges(2,6) * t256) + t248 * V_base(5) * (Icges(2,5) * t256 + Icges(2,6) * t259) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
