% Calculate kinetic energy for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:49:58
% EndTime: 2019-03-09 02:50:02
% DurationCPUTime: 4.46s
% Computational Cost: add. (1696->345), mult. (1883->476), div. (0->0), fcn. (1721->10), ass. (0->176)
t347 = Icges(4,4) + Icges(5,6);
t346 = Icges(4,1) + Icges(5,2);
t345 = Icges(4,2) + Icges(5,3);
t250 = pkin(9) + qJ(3);
t243 = cos(t250);
t344 = t347 * t243;
t241 = sin(t250);
t343 = t347 * t241;
t342 = Icges(5,4) - Icges(4,5);
t341 = Icges(5,5) - Icges(4,6);
t340 = t345 * t241 - t344;
t339 = t346 * t243 - t343;
t338 = Icges(5,1) + Icges(4,3);
t257 = sin(qJ(1));
t258 = cos(qJ(1));
t337 = t340 * t257 - t341 * t258;
t336 = t341 * t257 + t340 * t258;
t335 = -t342 * t257 + t339 * t258;
t334 = t339 * t257 + t342 * t258;
t333 = -t345 * t243 - t343;
t332 = t346 * t241 + t344;
t331 = t341 * t241 - t342 * t243;
t234 = -qJD(3) * t258 + V_base(5);
t235 = qJD(3) * t257 + V_base(4);
t244 = V_base(6) + qJD(1);
t330 = (t333 * t241 + t332 * t243) * t244 + (t336 * t241 + t335 * t243) * t235 + (t337 * t241 + t334 * t243) * t234;
t329 = (-t342 * t241 - t341 * t243) * t244 + (t338 * t257 + t331 * t258) * t235 + (t331 * t257 - t338 * t258) * t234;
t252 = sin(pkin(9));
t325 = pkin(2) * t252;
t251 = sin(pkin(10));
t324 = pkin(5) * t251;
t254 = cos(pkin(9));
t323 = pkin(2) * t254;
t253 = cos(pkin(10));
t322 = pkin(5) * t253;
t321 = Icges(2,4) * t257;
t320 = Icges(3,4) * t252;
t319 = Icges(3,4) * t254;
t314 = qJ(5) * t241;
t249 = pkin(10) + qJ(6);
t240 = sin(t249);
t313 = t240 * t258;
t242 = cos(t249);
t312 = t242 * t258;
t311 = t243 * t257;
t310 = t243 * t258;
t309 = t251 * t258;
t308 = t253 * t258;
t307 = t257 * t240;
t306 = t257 * t242;
t305 = t257 * t251;
t304 = t257 * t253;
t156 = -pkin(7) * t258 + t323 * t257;
t230 = t257 * pkin(1) - qJ(2) * t258;
t301 = -t156 - t230;
t284 = pkin(3) * t243 + qJ(4) * t241;
t191 = t284 * t258;
t200 = t257 * pkin(4) + qJ(5) * t310;
t300 = -t191 - t200;
t299 = qJD(4) * t241;
t298 = qJD(5) * t243;
t297 = qJD(6) * t243;
t296 = V_base(4) * t230 + V_base(3);
t295 = V_base(5) * pkin(6) + V_base(1);
t190 = t284 * t257;
t292 = -t190 + t301;
t291 = qJD(2) * t257 + t295;
t208 = pkin(3) * t241 - qJ(4) * t243;
t290 = -t208 - t314;
t201 = -pkin(4) * t258 + qJ(5) * t311;
t289 = -t201 + t292;
t288 = V_base(5) * t325 + t291;
t287 = rSges(3,1) * t254 - rSges(3,2) * t252;
t286 = rSges(4,1) * t243 - rSges(4,2) * t241;
t285 = -rSges(5,2) * t243 + rSges(5,3) * t241;
t283 = Icges(3,1) * t254 - t320;
t281 = -Icges(3,2) * t252 + t319;
t278 = Icges(3,5) * t254 - Icges(3,6) * t252;
t232 = pkin(1) * t258 + t257 * qJ(2);
t274 = -qJD(2) * t258 + t244 * t232 + V_base(2);
t273 = t234 * t208 + t258 * t299 + t288;
t272 = pkin(8) * t243 + t241 * t324;
t269 = t234 * t314 + t258 * t298 + t273;
t157 = pkin(7) * t257 + t323 * t258;
t268 = V_base(4) * t156 + (-t157 - t232) * V_base(5) + t296;
t267 = (-Icges(3,3) * t258 + t278 * t257) * V_base(5) + (Icges(3,3) * t257 + t278 * t258) * V_base(4) + (Icges(3,5) * t252 + Icges(3,6) * t254) * t244;
t266 = t244 * t157 + (-pkin(6) - t325) * V_base(4) + t274;
t265 = -qJD(4) * t243 + t235 * t190 + t268;
t264 = t244 * t191 + t257 * t299 + t266;
t263 = qJD(5) * t241 + t235 * t201 + t265;
t262 = t244 * t200 + t257 * t298 + t264;
t183 = -Icges(3,6) * t258 + t281 * t257;
t184 = Icges(3,6) * t257 + t281 * t258;
t185 = -Icges(3,5) * t258 + t283 * t257;
t186 = Icges(3,5) * t257 + t283 * t258;
t217 = Icges(3,2) * t254 + t320;
t218 = Icges(3,1) * t252 + t319;
t259 = (-t184 * t252 + t186 * t254) * V_base(4) + (-t183 * t252 + t185 * t254) * V_base(5) + (-t217 * t252 + t218 * t254) * t244;
t247 = Icges(2,4) * t258;
t233 = rSges(2,1) * t258 - t257 * rSges(2,2);
t231 = t257 * rSges(2,1) + rSges(2,2) * t258;
t225 = Icges(2,1) * t258 - t321;
t224 = Icges(2,1) * t257 + t247;
t223 = -Icges(2,2) * t257 + t247;
t222 = Icges(2,2) * t258 + t321;
t221 = Icges(2,5) * t258 - Icges(2,6) * t257;
t220 = Icges(2,5) * t257 + Icges(2,6) * t258;
t219 = rSges(3,1) * t252 + rSges(3,2) * t254;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t212 = qJD(6) * t241 + t244;
t210 = rSges(4,1) * t241 + rSges(4,2) * t243;
t209 = -rSges(5,2) * t241 - rSges(5,3) * t243;
t197 = t258 * t297 + t235;
t196 = t257 * t297 + t234;
t195 = t241 * t305 - t308;
t194 = t241 * t304 + t309;
t193 = t241 * t309 + t304;
t192 = t241 * t308 - t305;
t188 = t257 * rSges(3,3) + t287 * t258;
t187 = -rSges(3,3) * t258 + t287 * t257;
t180 = t241 * t307 - t312;
t179 = t241 * t306 + t313;
t178 = t241 * t313 + t306;
t177 = t241 * t312 - t307;
t174 = -rSges(5,1) * t258 + t285 * t257;
t173 = t257 * rSges(5,1) + t285 * t258;
t172 = t257 * rSges(4,3) + t286 * t258;
t171 = -rSges(4,3) * t258 + t286 * t257;
t155 = pkin(8) * t241 - t243 * t324;
t154 = V_base(5) * rSges(2,3) - t231 * t244 + t295;
t153 = t233 * t244 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = rSges(6,3) * t241 + (-rSges(6,1) * t251 - rSges(6,2) * t253) * t243;
t151 = Icges(6,5) * t241 + (-Icges(6,1) * t251 - Icges(6,4) * t253) * t243;
t150 = Icges(6,6) * t241 + (-Icges(6,4) * t251 - Icges(6,2) * t253) * t243;
t149 = Icges(6,3) * t241 + (-Icges(6,5) * t251 - Icges(6,6) * t253) * t243;
t147 = t231 * V_base(4) - t233 * V_base(5) + V_base(3);
t144 = rSges(7,3) * t241 + (-rSges(7,1) * t240 - rSges(7,2) * t242) * t243;
t143 = Icges(7,5) * t241 + (-Icges(7,1) * t240 - Icges(7,4) * t242) * t243;
t142 = Icges(7,6) * t241 + (-Icges(7,4) * t240 - Icges(7,2) * t242) * t243;
t141 = Icges(7,3) * t241 + (-Icges(7,5) * t240 - Icges(7,6) * t242) * t243;
t140 = t272 * t257 - t322 * t258;
t139 = t322 * t257 + t272 * t258;
t138 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t311;
t137 = t193 * rSges(6,1) + t192 * rSges(6,2) + rSges(6,3) * t310;
t136 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t311;
t135 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t310;
t134 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t311;
t133 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t310;
t132 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t311;
t131 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t310;
t130 = rSges(7,1) * t180 + rSges(7,2) * t179 + rSges(7,3) * t311;
t129 = t178 * rSges(7,1) + t177 * rSges(7,2) + rSges(7,3) * t310;
t128 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t311;
t127 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t310;
t126 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t311;
t125 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t310;
t124 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t311;
t123 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t310;
t122 = t219 * V_base(5) + (-t187 - t230) * t244 + t291;
t121 = t244 * t188 + (-pkin(6) - t219) * V_base(4) + t274;
t120 = t187 * V_base(4) + (-t188 - t232) * V_base(5) + t296;
t119 = t210 * t234 + (-t171 + t301) * t244 + t288;
t118 = t244 * t172 - t235 * t210 + t266;
t117 = t171 * t235 - t172 * t234 + t268;
t116 = t209 * t234 + (-t174 + t292) * t244 + t273;
t115 = t244 * t173 + (-t208 - t209) * t235 + t264;
t114 = t174 * t235 + (-t173 - t191) * t234 + t265;
t113 = t152 * t234 + (-t138 + t289) * t244 + t269;
t112 = t244 * t137 + (-t152 + t290) * t235 + t262;
t111 = t138 * t235 + (-t137 + t300) * t234 + t263;
t110 = -t130 * t212 + t144 * t196 + t155 * t234 + (-t140 + t289) * t244 + t269;
t109 = t212 * t129 + t244 * t139 - t197 * t144 + (-t155 + t290) * t235 + t262;
t108 = -t129 * t196 + t130 * t197 + t140 * t235 + (-t139 + t300) * t234 + t263;
t1 = t197 * ((t123 * t310 + t177 * t125 + t178 * t127) * t197 + (t124 * t310 + t177 * t126 + t178 * t128) * t196 + (t141 * t310 + t177 * t142 + t178 * t143) * t212) / 0.2e1 + t196 * ((t123 * t311 + t125 * t179 + t127 * t180) * t197 + (t124 * t311 + t179 * t126 + t180 * t128) * t196 + (t141 * t311 + t142 * t179 + t143 * t180) * t212) / 0.2e1 + t212 * ((t123 * t197 + t124 * t196 + t141 * t212) * t241 + ((-t125 * t242 - t127 * t240) * t197 + (-t126 * t242 - t128 * t240) * t196 + (-t142 * t242 - t143 * t240) * t212) * t243) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + ((t131 * t311 + t133 * t194 + t135 * t195) * t235 + (t132 * t311 + t134 * t194 + t136 * t195) * t234 + (t149 * t311 + t150 * t194 + t151 * t195) * t244 - t329 * t258 + t330 * t257) * t234 / 0.2e1 + ((t131 * t310 + t192 * t133 + t193 * t135) * t235 + (t132 * t310 + t192 * t134 + t193 * t136) * t234 + (t149 * t310 + t192 * t150 + t193 * t151) * t244 + t330 * t258 + t329 * t257) * t235 / 0.2e1 + (t221 * t244 + t267 * t257 + t259 * t258 + (-t257 * t222 + t224 * t258 + Icges(1,4)) * V_base(5) + (-t257 * t223 + t225 * t258 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t220 * t244 + t259 * t257 - t267 * t258 + (t222 * t258 + t257 * t224 + Icges(1,2)) * V_base(5) + (t223 * t258 + t257 * t225 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t183 * t254 + t185 * t252 + t220) * V_base(5) + (t184 * t254 + t186 * t252 + t221) * V_base(4) + ((-t133 * t253 - t135 * t251 - t336) * t243 + (t131 + t335) * t241) * t235 + ((-t134 * t253 - t136 * t251 - t337) * t243 + (t132 + t334) * t241) * t234 + (t217 * t254 + t218 * t252 + Icges(2,3) + (-t150 * t253 - t151 * t251 - t333) * t243 + (t149 + t332) * t241) * t244) * t244 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
