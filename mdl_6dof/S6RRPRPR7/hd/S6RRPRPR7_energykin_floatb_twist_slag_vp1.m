% Calculate kinetic energy for
% S6RRPRPR7
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:12
% EndTime: 2019-03-09 10:45:16
% DurationCPUTime: 4.36s
% Computational Cost: add. (1787->341), mult. (2604->477), div. (0->0), fcn. (2698->10), ass. (0->163)
t349 = Icges(3,4) - Icges(4,5);
t348 = Icges(3,1) + Icges(4,1);
t347 = Icges(3,2) + Icges(4,3);
t268 = sin(qJ(2));
t346 = t349 * t268;
t272 = cos(qJ(2));
t345 = t349 * t272;
t344 = Icges(4,4) + Icges(3,5);
t343 = Icges(3,6) - Icges(4,6);
t342 = t347 * t268 - t345;
t341 = t348 * t272 - t346;
t340 = Icges(4,2) + Icges(3,3);
t339 = Icges(5,3) + Icges(6,3);
t269 = sin(qJ(1));
t273 = cos(qJ(1));
t338 = t342 * t269 + t343 * t273;
t337 = -t343 * t269 + t342 * t273;
t336 = t341 * t269 - t344 * t273;
t335 = t344 * t269 + t341 * t273;
t334 = -t347 * t272 - t346;
t333 = t348 * t268 + t345;
t332 = -t343 * t268 + t344 * t272;
t307 = qJ(4) + pkin(10);
t258 = sin(t307);
t302 = cos(t307);
t296 = t268 * t302;
t216 = -t272 * t258 + t296;
t187 = t216 * t269;
t215 = t268 * t258 + t272 * t302;
t188 = t215 * t269;
t271 = cos(qJ(4));
t267 = sin(qJ(4));
t312 = t267 * t272;
t224 = t268 * t271 - t312;
t211 = t224 * t269;
t313 = t267 * t268;
t289 = t271 * t272 + t313;
t212 = t289 * t269;
t331 = Icges(5,5) * t212 + Icges(6,5) * t188 + Icges(5,6) * t211 + Icges(6,6) * t187 + t339 * t273;
t310 = t272 * t273;
t189 = t258 * t310 - t273 * t296;
t190 = t215 * t273;
t213 = t224 * t273;
t214 = t289 * t273;
t330 = Icges(5,5) * t214 + Icges(6,5) * t190 + Icges(5,6) * t213 - Icges(6,6) * t189 - t269 * t339;
t329 = Icges(5,5) * t224 + Icges(6,5) * t216 - Icges(5,6) * t289 - Icges(6,6) * t215;
t253 = -qJD(2) * t273 + V_base(5);
t254 = qJD(2) * t269 + V_base(4);
t259 = V_base(6) + qJD(1);
t328 = (t268 * t334 + t272 * t333) * t259 + (t268 * t337 + t272 * t335) * t254 + (t268 * t338 + t272 * t336) * t253;
t327 = (t344 * t268 + t343 * t272) * t259 + (t269 * t340 + t332 * t273) * t254 + (t332 * t269 - t340 * t273) * t253;
t321 = pkin(3) * t268;
t320 = pkin(4) * t271;
t318 = Icges(2,4) * t269;
t297 = pkin(2) * t272 + qJ(3) * t268;
t218 = t297 * t269;
t251 = t269 * pkin(1) - pkin(7) * t273;
t309 = -t218 - t251;
t308 = qJD(3) * t268;
t306 = V_base(5) * pkin(6) + V_base(1);
t227 = t269 * t272 * pkin(3) + pkin(8) * t273;
t303 = -t227 + t309;
t285 = pkin(4) * t313 + t272 * t320;
t164 = qJ(5) * t273 + t269 * t285;
t301 = -t164 + t303;
t246 = pkin(2) * t268 - qJ(3) * t272;
t300 = t253 * t246 + t273 * t308 + t306;
t299 = rSges(3,1) * t272 - rSges(3,2) * t268;
t298 = rSges(4,1) * t272 + rSges(4,3) * t268;
t288 = t253 * t321 + t300;
t221 = qJD(4) * t273 + t253;
t222 = -qJD(4) * t269 + t254;
t252 = pkin(1) * t273 + t269 * pkin(7);
t287 = -V_base(4) * pkin(6) + t259 * t252 + V_base(2);
t286 = V_base(4) * t251 - t252 * V_base(5) + V_base(3);
t219 = t297 * t273;
t282 = t259 * t219 + t269 * t308 + t287;
t191 = -pkin(4) * t312 + t268 * t320;
t281 = -qJD(5) * t269 + t221 * t191 + t288;
t280 = -qJD(3) * t272 + t254 * t218 + t286;
t228 = pkin(3) * t310 - t269 * pkin(8);
t279 = t259 * t228 + (-t246 - t321) * t254 + t282;
t278 = t254 * t227 + (-t219 - t228) * t253 + t280;
t277 = t222 * t164 + t278;
t165 = -qJ(5) * t269 + t273 * t285;
t276 = qJD(5) * t273 + t259 * t165 + t279;
t270 = cos(qJ(6));
t266 = sin(qJ(6));
t263 = Icges(2,4) * t273;
t250 = rSges(2,1) * t273 - t269 * rSges(2,2);
t249 = t269 * rSges(2,1) + rSges(2,2) * t273;
t248 = rSges(3,1) * t268 + rSges(3,2) * t272;
t247 = rSges(4,1) * t268 - rSges(4,3) * t272;
t245 = Icges(2,1) * t273 - t318;
t244 = Icges(2,1) * t269 + t263;
t241 = -Icges(2,2) * t269 + t263;
t240 = Icges(2,2) * t273 + t318;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = t269 * rSges(3,3) + t273 * t299;
t208 = t269 * rSges(4,2) + t273 * t298;
t207 = -rSges(3,3) * t273 + t269 * t299;
t206 = -rSges(4,2) * t273 + t269 * t298;
t185 = qJD(6) * t215 + t259;
t184 = V_base(5) * rSges(2,3) - t249 * t259 + t306;
t183 = t250 * t259 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t182 = t249 * V_base(4) - t250 * V_base(5) + V_base(3);
t181 = rSges(5,1) * t224 - rSges(5,2) * t289;
t180 = Icges(5,1) * t224 - Icges(5,4) * t289;
t179 = Icges(5,4) * t224 - Icges(5,2) * t289;
t177 = t190 * t270 - t266 * t269;
t176 = -t190 * t266 - t269 * t270;
t175 = t188 * t270 + t266 * t273;
t174 = -t188 * t266 + t270 * t273;
t172 = qJD(6) * t189 + t222;
t171 = -qJD(6) * t187 + t221;
t170 = pkin(5) * t216 + pkin(9) * t215;
t169 = rSges(6,1) * t216 - rSges(6,2) * t215;
t168 = Icges(6,1) * t216 - Icges(6,4) * t215;
t167 = Icges(6,4) * t216 - Icges(6,2) * t215;
t163 = rSges(5,1) * t214 + rSges(5,2) * t213 - rSges(5,3) * t269;
t162 = t212 * rSges(5,1) + t211 * rSges(5,2) + rSges(5,3) * t273;
t160 = Icges(5,1) * t214 + Icges(5,4) * t213 - Icges(5,5) * t269;
t159 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t273;
t158 = Icges(5,4) * t214 + Icges(5,2) * t213 - Icges(5,6) * t269;
t157 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t273;
t154 = pkin(5) * t190 + pkin(9) * t189;
t153 = pkin(5) * t188 - pkin(9) * t187;
t151 = rSges(6,1) * t190 - rSges(6,2) * t189 - rSges(6,3) * t269;
t150 = t188 * rSges(6,1) + t187 * rSges(6,2) + rSges(6,3) * t273;
t149 = Icges(6,1) * t190 - Icges(6,4) * t189 - Icges(6,5) * t269;
t148 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t273;
t147 = Icges(6,4) * t190 - Icges(6,2) * t189 - Icges(6,6) * t269;
t146 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t273;
t143 = rSges(7,3) * t215 + (rSges(7,1) * t270 - rSges(7,2) * t266) * t216;
t142 = Icges(7,5) * t215 + (Icges(7,1) * t270 - Icges(7,4) * t266) * t216;
t141 = Icges(7,6) * t215 + (Icges(7,4) * t270 - Icges(7,2) * t266) * t216;
t140 = Icges(7,3) * t215 + (Icges(7,5) * t270 - Icges(7,6) * t266) * t216;
t139 = t248 * t253 + (-t207 - t251) * t259 + t306;
t138 = t209 * t259 - t248 * t254 + t287;
t137 = t207 * t254 - t209 * t253 + t286;
t136 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t189;
t135 = rSges(7,1) * t175 + rSges(7,2) * t174 - rSges(7,3) * t187;
t134 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t189;
t133 = Icges(7,1) * t175 + Icges(7,4) * t174 - Icges(7,5) * t187;
t132 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t189;
t131 = Icges(7,4) * t175 + Icges(7,2) * t174 - Icges(7,6) * t187;
t130 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t189;
t129 = Icges(7,5) * t175 + Icges(7,6) * t174 - Icges(7,3) * t187;
t128 = t247 * t253 + (-t206 + t309) * t259 + t300;
t127 = t208 * t259 + (-t246 - t247) * t254 + t282;
t126 = t206 * t254 + (-t208 - t219) * t253 + t280;
t125 = t181 * t221 + (-t162 + t303) * t259 + t288;
t124 = t163 * t259 - t181 * t222 + t279;
t123 = t162 * t222 - t163 * t221 + t278;
t122 = t169 * t221 + (-t150 + t301) * t259 + t281;
t121 = t151 * t259 + (-t169 - t191) * t222 + t276;
t120 = t150 * t222 + (-t151 - t165) * t221 + t277;
t119 = -t135 * t185 + t143 * t171 + t170 * t221 + (-t153 + t301) * t259 + t281;
t118 = t136 * t185 - t143 * t172 + t154 * t259 + (-t170 - t191) * t222 + t276;
t117 = t135 * t172 - t136 * t171 + t153 * t222 + (-t154 - t165) * t221 + t277;
t1 = m(2) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t185 * ((t129 * t171 + t130 * t172 + t140 * t185) * t215 + ((-t132 * t266 + t134 * t270) * t172 + (-t131 * t266 + t133 * t270) * t171 + (-t141 * t266 + t142 * t270) * t185) * t216) / 0.2e1 + t172 * ((t189 * t130 + t176 * t132 + t177 * t134) * t172 + (t129 * t189 + t131 * t176 + t133 * t177) * t171 + (t140 * t189 + t141 * t176 + t142 * t177) * t185) / 0.2e1 + t171 * ((-t130 * t187 + t132 * t174 + t134 * t175) * t172 + (-t187 * t129 + t174 * t131 + t175 * t133) * t171 + (-t140 * t187 + t141 * t174 + t142 * t175) * t185) / 0.2e1 + m(3) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(5) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(7) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(6) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + ((t187 * t167 + t188 * t168 + t211 * t179 + t212 * t180 + t273 * t329) * t259 + (t187 * t147 + t188 * t149 + t211 * t158 + t212 * t160 + t273 * t330) * t222 + (t187 * t146 + t188 * t148 + t211 * t157 + t212 * t159 + t331 * t273) * t221) * t221 / 0.2e1 + ((-t189 * t167 + t168 * t190 + t179 * t213 + t180 * t214 - t269 * t329) * t259 + (-t189 * t147 + t190 * t149 + t213 * t158 + t214 * t160 - t330 * t269) * t222 + (-t146 * t189 + t148 * t190 + t157 * t213 + t159 * t214 - t269 * t331) * t221) * t222 / 0.2e1 + (t328 * t269 - t327 * t273) * t253 / 0.2e1 + (t327 * t269 + t328 * t273) * t254 / 0.2e1 + ((-t269 * t240 + t244 * t273 + Icges(1,4)) * V_base(5) + (-t269 * t241 + t273 * t245 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t273 * t240 + t269 * t244 + Icges(1,2)) * V_base(5) + (t241 * t273 + t269 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t268 * t335 - t272 * t337) * t254 + (t268 * t336 - t272 * t338) * t253 + (-t147 * t215 + t149 * t216 - t158 * t289 + t160 * t224) * t222 + (-t146 * t215 + t148 * t216 - t157 * t289 + t159 * t224) * t221 + (-t215 * t167 + t216 * t168 - t289 * t179 + t224 * t180 + t268 * t333 - t272 * t334 + Icges(2,3)) * t259) * t259 / 0.2e1 + t259 * V_base(4) * (Icges(2,5) * t273 - Icges(2,6) * t269) + t259 * V_base(5) * (Icges(2,5) * t269 + Icges(2,6) * t273) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
