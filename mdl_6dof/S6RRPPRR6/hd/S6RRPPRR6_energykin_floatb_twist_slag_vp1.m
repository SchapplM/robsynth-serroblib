% Calculate kinetic energy for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:12:57
% EndTime: 2019-03-09 09:13:02
% DurationCPUTime: 4.38s
% Computational Cost: add. (1769->340), mult. (2586->483), div. (0->0), fcn. (2680->10), ass. (0->164)
t345 = Icges(3,4) - Icges(4,5);
t344 = Icges(3,1) + Icges(4,1);
t343 = Icges(3,2) + Icges(4,3);
t269 = sin(qJ(2));
t342 = t345 * t269;
t272 = cos(qJ(2));
t341 = t345 * t272;
t340 = Icges(4,4) + Icges(3,5);
t339 = Icges(3,6) - Icges(4,6);
t338 = t343 * t269 - t341;
t337 = t344 * t272 - t342;
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t336 = t338 * t270 + t339 * t273;
t335 = -t339 * t270 + t338 * t273;
t334 = t337 * t270 - t340 * t273;
t333 = t340 * t270 + t337 * t273;
t332 = -t343 * t272 - t342;
t331 = t344 * t269 + t341;
t330 = t339 * t269 - t340 * t272;
t329 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t253 = -qJD(2) * t273 + V_base(5);
t254 = qJD(2) * t270 + V_base(4);
t259 = V_base(6) + qJD(1);
t328 = (t332 * t269 + t331 * t272) * t259 + (t335 * t269 + t333 * t272) * t254 + (t336 * t269 + t334 * t272) * t253;
t308 = pkin(10) + qJ(5);
t258 = sin(t308);
t302 = cos(t308);
t296 = t269 * t302;
t217 = -t272 * t258 + t296;
t266 = cos(pkin(10));
t265 = sin(pkin(10));
t315 = t265 * t272;
t222 = t266 * t269 - t315;
t210 = t222 * t270;
t316 = t265 * t269;
t289 = t266 * t272 + t316;
t211 = t289 * t270;
t212 = t222 * t273;
t213 = t289 * t273;
t327 = (Icges(5,5) * t222 - Icges(5,6) * t289 - t340 * t269 - t339 * t272) * t259 + (Icges(5,5) * t213 + Icges(5,6) * t212 - t329 * t270 + t330 * t273) * t254 + (Icges(5,5) * t211 + Icges(5,6) * t210 + t330 * t270 + t329 * t273) * t253;
t323 = pkin(3) * t269;
t322 = pkin(4) * t266;
t321 = Icges(2,4) * t270;
t313 = t272 * t273;
t297 = pkin(2) * t272 + qJ(3) * t269;
t218 = t297 * t270;
t251 = t270 * pkin(1) - pkin(7) * t273;
t311 = -t218 - t251;
t219 = t297 * t273;
t228 = pkin(3) * t313 - t270 * qJ(4);
t310 = -t219 - t228;
t309 = qJD(3) * t269;
t307 = V_base(5) * pkin(6) + V_base(1);
t227 = t270 * t272 * pkin(3) + qJ(4) * t273;
t304 = -t227 + t311;
t246 = pkin(2) * t269 - qJ(3) * t272;
t303 = -t246 - t323;
t286 = pkin(4) * t316 + t272 * t322;
t164 = pkin(8) * t273 + t270 * t286;
t301 = -t164 + t304;
t300 = t253 * t246 + t273 * t309 + t307;
t299 = rSges(3,1) * t272 - rSges(3,2) * t269;
t298 = rSges(4,1) * t272 + rSges(4,3) * t269;
t223 = qJD(5) * t273 + t253;
t224 = -qJD(5) * t270 + t254;
t252 = pkin(1) * t273 + t270 * pkin(7);
t288 = -V_base(4) * pkin(6) + t259 * t252 + V_base(2);
t287 = V_base(4) * t251 - t252 * V_base(5) + V_base(3);
t283 = t259 * t219 + t270 * t309 + t288;
t216 = t269 * t258 + t272 * t302;
t282 = -qJD(4) * t270 + t253 * t323 + t300;
t191 = -pkin(4) * t315 + t269 * t322;
t281 = t253 * t191 + t282;
t280 = -qJD(3) * t272 + t254 * t218 + t287;
t279 = qJD(4) * t273 + t259 * t228 + t283;
t278 = t254 * t227 + t280;
t165 = -pkin(8) * t270 + t273 * t286;
t277 = t254 * t164 + (-t165 + t310) * t253 + t278;
t276 = t259 * t165 + (-t191 + t303) * t254 + t279;
t271 = cos(qJ(6));
t268 = sin(qJ(6));
t263 = Icges(2,4) * t273;
t250 = rSges(2,1) * t273 - t270 * rSges(2,2);
t249 = t270 * rSges(2,1) + rSges(2,2) * t273;
t248 = rSges(3,1) * t269 + rSges(3,2) * t272;
t247 = rSges(4,1) * t269 - rSges(4,3) * t272;
t245 = Icges(2,1) * t273 - t321;
t244 = Icges(2,1) * t270 + t263;
t241 = -Icges(2,2) * t270 + t263;
t240 = Icges(2,2) * t273 + t321;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = t270 * rSges(3,3) + t273 * t299;
t208 = t270 * rSges(4,2) + t273 * t298;
t207 = -rSges(3,3) * t273 + t270 * t299;
t206 = -rSges(4,2) * t273 + t270 * t298;
t190 = t216 * t273;
t189 = t258 * t313 - t273 * t296;
t188 = t216 * t270;
t187 = t217 * t270;
t185 = qJD(6) * t216 + t259;
t184 = V_base(5) * rSges(2,3) - t249 * t259 + t307;
t183 = t250 * t259 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t182 = t249 * V_base(4) - t250 * V_base(5) + V_base(3);
t180 = rSges(5,1) * t222 - rSges(5,2) * t289;
t179 = Icges(5,1) * t222 - Icges(5,4) * t289;
t178 = Icges(5,4) * t222 - Icges(5,2) * t289;
t176 = t190 * t271 - t268 * t270;
t175 = -t190 * t268 - t270 * t271;
t174 = t188 * t271 + t268 * t273;
t173 = -t188 * t268 + t271 * t273;
t172 = qJD(6) * t189 + t224;
t171 = -qJD(6) * t187 + t223;
t170 = pkin(5) * t217 + pkin(9) * t216;
t169 = rSges(6,1) * t217 - rSges(6,2) * t216;
t168 = Icges(6,1) * t217 - Icges(6,4) * t216;
t167 = Icges(6,4) * t217 - Icges(6,2) * t216;
t166 = Icges(6,5) * t217 - Icges(6,6) * t216;
t162 = rSges(5,1) * t213 + rSges(5,2) * t212 - rSges(5,3) * t270;
t161 = t211 * rSges(5,1) + t210 * rSges(5,2) + rSges(5,3) * t273;
t160 = Icges(5,1) * t213 + Icges(5,4) * t212 - Icges(5,5) * t270;
t159 = Icges(5,1) * t211 + Icges(5,4) * t210 + Icges(5,5) * t273;
t158 = Icges(5,4) * t213 + Icges(5,2) * t212 - Icges(5,6) * t270;
t157 = Icges(5,4) * t211 + Icges(5,2) * t210 + Icges(5,6) * t273;
t154 = pkin(5) * t190 + pkin(9) * t189;
t153 = pkin(5) * t188 - pkin(9) * t187;
t151 = rSges(6,1) * t190 - rSges(6,2) * t189 - rSges(6,3) * t270;
t150 = t188 * rSges(6,1) + t187 * rSges(6,2) + rSges(6,3) * t273;
t149 = Icges(6,1) * t190 - Icges(6,4) * t189 - Icges(6,5) * t270;
t148 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t273;
t147 = Icges(6,4) * t190 - Icges(6,2) * t189 - Icges(6,6) * t270;
t146 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t273;
t145 = Icges(6,5) * t190 - Icges(6,6) * t189 - Icges(6,3) * t270;
t144 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t273;
t143 = rSges(7,3) * t216 + (rSges(7,1) * t271 - rSges(7,2) * t268) * t217;
t142 = Icges(7,5) * t216 + (Icges(7,1) * t271 - Icges(7,4) * t268) * t217;
t141 = Icges(7,6) * t216 + (Icges(7,4) * t271 - Icges(7,2) * t268) * t217;
t140 = Icges(7,3) * t216 + (Icges(7,5) * t271 - Icges(7,6) * t268) * t217;
t139 = t248 * t253 + (-t207 - t251) * t259 + t307;
t138 = t209 * t259 - t248 * t254 + t288;
t137 = t207 * t254 - t209 * t253 + t287;
t136 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t189;
t135 = rSges(7,1) * t174 + rSges(7,2) * t173 - rSges(7,3) * t187;
t134 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t189;
t133 = Icges(7,1) * t174 + Icges(7,4) * t173 - Icges(7,5) * t187;
t132 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t189;
t131 = Icges(7,4) * t174 + Icges(7,2) * t173 - Icges(7,6) * t187;
t130 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t189;
t129 = Icges(7,5) * t174 + Icges(7,6) * t173 - Icges(7,3) * t187;
t128 = t247 * t253 + (-t206 + t311) * t259 + t300;
t127 = t208 * t259 + (-t246 - t247) * t254 + t283;
t126 = t206 * t254 + (-t208 - t219) * t253 + t280;
t125 = t180 * t253 + (-t161 + t304) * t259 + t282;
t124 = t162 * t259 + (-t180 + t303) * t254 + t279;
t123 = t161 * t254 + (-t162 + t310) * t253 + t278;
t122 = t169 * t223 + (-t150 + t301) * t259 + t281;
t121 = t151 * t259 - t169 * t224 + t276;
t120 = t150 * t224 - t151 * t223 + t277;
t119 = -t135 * t185 + t143 * t171 + t170 * t223 + (-t153 + t301) * t259 + t281;
t118 = t136 * t185 - t143 * t172 + t154 * t259 - t170 * t224 + t276;
t117 = t135 * t172 - t136 * t171 + t153 * t224 - t154 * t223 + t277;
t1 = m(3) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + t185 * ((t129 * t171 + t130 * t172 + t140 * t185) * t216 + ((-t132 * t268 + t134 * t271) * t172 + (-t131 * t268 + t133 * t271) * t171 + (-t141 * t268 + t142 * t271) * t185) * t217) / 0.2e1 + m(6) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(5) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(7) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t223 * ((t145 * t273 + t187 * t147 + t188 * t149) * t224 + (t144 * t273 + t187 * t146 + t188 * t148) * t223 + (t166 * t273 + t187 * t167 + t188 * t168) * t259) / 0.2e1 + t224 * ((-t145 * t270 - t147 * t189 + t149 * t190) * t224 + (-t144 * t270 - t146 * t189 + t148 * t190) * t223 + (-t166 * t270 - t167 * t189 + t168 * t190) * t259) / 0.2e1 + m(2) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + t171 * ((-t130 * t187 + t132 * t173 + t134 * t174) * t172 + (-t129 * t187 + t131 * t173 + t133 * t174) * t171 + (-t140 * t187 + t141 * t173 + t142 * t174) * t185) / 0.2e1 + t172 * ((t130 * t189 + t132 * t175 + t134 * t176) * t172 + (t129 * t189 + t131 * t175 + t133 * t176) * t171 + (t140 * t189 + t141 * t175 + t142 * t176) * t185) / 0.2e1 + ((-t270 * t240 + t244 * t273 + Icges(1,4)) * V_base(5) + (-t270 * t241 + t245 * t273 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t240 * t273 + t270 * t244 + Icges(1,2)) * V_base(5) + (t241 * t273 + t270 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t210 * t158 + t211 * t160) * t254 + (t210 * t157 + t211 * t159) * t253 + (t210 * t178 + t211 * t179) * t259 + t327 * t273 + t328 * t270) * t253 / 0.2e1 + ((t158 * t212 + t160 * t213) * t254 + (t157 * t212 + t159 * t213) * t253 + (t178 * t212 + t179 * t213) * t259 + t328 * t273 - t327 * t270) * t254 / 0.2e1 + ((-t147 * t216 + t149 * t217) * t224 + (-t146 * t216 + t148 * t217) * t223 + (-t158 * t289 + t160 * t222 + t333 * t269 - t335 * t272) * t254 + (-t157 * t289 + t159 * t222 + t334 * t269 - t336 * t272) * t253 + (-t167 * t216 + t168 * t217 - t178 * t289 + t179 * t222 + t331 * t269 - t332 * t272 + Icges(2,3)) * t259) * t259 / 0.2e1 + t259 * V_base(4) * (Icges(2,5) * t273 - Icges(2,6) * t270) + t259 * V_base(5) * (Icges(2,5) * t270 + Icges(2,6) * t273) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
