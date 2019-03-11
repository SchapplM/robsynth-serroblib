% Calculate kinetic energy for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:50
% EndTime: 2019-03-09 08:09:54
% DurationCPUTime: 4.00s
% Computational Cost: add. (1722->338), mult. (1909->462), div. (0->0), fcn. (1747->10), ass. (0->172)
t348 = Icges(4,4) + Icges(5,6);
t347 = Icges(4,1) + Icges(5,2);
t346 = Icges(4,2) + Icges(5,3);
t250 = qJ(2) + pkin(9);
t243 = cos(t250);
t345 = t348 * t243;
t241 = sin(t250);
t344 = t348 * t241;
t343 = Icges(5,4) - Icges(4,5);
t342 = Icges(5,5) - Icges(4,6);
t341 = t346 * t241 - t345;
t340 = t347 * t243 - t344;
t256 = sin(qJ(1));
t258 = cos(qJ(1));
t339 = t256 * t341 - t258 * t342;
t338 = t256 * t342 + t258 * t341;
t337 = -t256 * t343 + t340 * t258;
t336 = t340 * t256 + t258 * t343;
t335 = -t346 * t243 - t344;
t334 = t347 * t241 + t345;
t333 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t255 = sin(qJ(2));
t257 = cos(qJ(2));
t332 = Icges(3,5) * t257 - Icges(3,6) * t255 + t342 * t241 - t243 * t343;
t319 = Icges(3,4) * t257;
t283 = -Icges(3,2) * t255 + t319;
t184 = -Icges(3,6) * t258 + t256 * t283;
t185 = Icges(3,6) * t256 + t258 * t283;
t320 = Icges(3,4) * t255;
t285 = Icges(3,1) * t257 - t320;
t186 = -Icges(3,5) * t258 + t256 * t285;
t187 = Icges(3,5) * t256 + t258 * t285;
t220 = Icges(3,2) * t257 + t320;
t223 = Icges(3,1) * t255 + t319;
t235 = -qJD(2) * t258 + V_base(5);
t236 = qJD(2) * t256 + V_base(4);
t244 = V_base(6) + qJD(1);
t331 = (-t220 * t255 + t223 * t257 + t241 * t335 + t243 * t334) * t244 + (-t185 * t255 + t187 * t257 + t241 * t338 + t243 * t337) * t236 + (-t184 * t255 + t186 * t257 + t241 * t339 + t336 * t243) * t235;
t330 = (Icges(3,5) * t255 + Icges(3,6) * t257 - t241 * t343 - t342 * t243) * t244 + (t256 * t333 + t258 * t332) * t236 + (t256 * t332 - t258 * t333) * t235;
t326 = pkin(2) * t255;
t251 = sin(pkin(10));
t325 = pkin(5) * t251;
t324 = pkin(2) * t257;
t252 = cos(pkin(10));
t323 = pkin(5) * t252;
t321 = Icges(2,4) * t256;
t314 = qJ(5) * t241;
t249 = pkin(10) + qJ(6);
t240 = sin(t249);
t313 = t240 * t258;
t242 = cos(t249);
t312 = t242 * t258;
t311 = t243 * t256;
t310 = t243 * t258;
t309 = t251 * t258;
t308 = t252 * t258;
t307 = t256 * t240;
t306 = t256 * t242;
t305 = t256 * t251;
t304 = t256 * t252;
t157 = -qJ(3) * t258 + t256 * t324;
t233 = t256 * pkin(1) - pkin(7) * t258;
t302 = -t157 - t233;
t158 = qJ(3) * t256 + t258 * t324;
t286 = pkin(3) * t243 + qJ(4) * t241;
t189 = t286 * t258;
t301 = -t158 - t189;
t300 = qJD(4) * t241;
t299 = qJD(5) * t243;
t298 = qJD(6) * t243;
t297 = V_base(5) * pkin(6) + V_base(1);
t188 = t286 * t256;
t294 = -t188 + t302;
t200 = t256 * pkin(4) + qJ(5) * t310;
t293 = -t200 + t301;
t208 = pkin(3) * t241 - qJ(4) * t243;
t292 = -t208 - t326;
t201 = -pkin(4) * t258 + qJ(5) * t311;
t291 = -t201 + t294;
t290 = qJD(3) * t256 + t235 * t326 + t297;
t289 = rSges(3,1) * t257 - rSges(3,2) * t255;
t288 = rSges(4,1) * t243 - rSges(4,2) * t241;
t287 = -rSges(5,2) * t243 + rSges(5,3) * t241;
t276 = t292 - t314;
t234 = pkin(1) * t258 + t256 * pkin(7);
t275 = -V_base(4) * pkin(6) + t244 * t234 + V_base(2);
t274 = V_base(4) * t233 - t234 * V_base(5) + V_base(3);
t273 = t235 * t208 + t258 * t300 + t290;
t272 = t236 * t157 + t274;
t271 = pkin(8) * t243 + t241 * t325;
t267 = t235 * t314 + t258 * t299 + t273;
t266 = -qJD(3) * t258 + t244 * t158 + t275;
t265 = -qJD(4) * t243 + t236 * t188 + t272;
t264 = t244 * t189 + t256 * t300 + t266;
t263 = qJD(5) * t241 + t236 * t201 + t265;
t262 = t244 * t200 + t256 * t299 + t264;
t247 = Icges(2,4) * t258;
t232 = rSges(2,1) * t258 - t256 * rSges(2,2);
t231 = t256 * rSges(2,1) + rSges(2,2) * t258;
t230 = rSges(3,1) * t255 + rSges(3,2) * t257;
t225 = Icges(2,1) * t258 - t321;
t224 = Icges(2,1) * t256 + t247;
t222 = -Icges(2,2) * t256 + t247;
t221 = Icges(2,2) * t258 + t321;
t216 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t215 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t214 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = qJD(6) * t241 + t244;
t210 = rSges(4,1) * t241 + rSges(4,2) * t243;
t209 = -rSges(5,2) * t241 - rSges(5,3) * t243;
t197 = t258 * t298 + t236;
t196 = t256 * t298 + t235;
t195 = t241 * t305 - t308;
t194 = t241 * t304 + t309;
t193 = t241 * t309 + t304;
t192 = t241 * t308 - t305;
t191 = t256 * rSges(3,3) + t258 * t289;
t190 = -rSges(3,3) * t258 + t256 * t289;
t180 = t241 * t307 - t312;
t179 = t241 * t306 + t313;
t178 = t241 * t313 + t306;
t177 = t241 * t312 - t307;
t174 = -rSges(5,1) * t258 + t256 * t287;
t173 = t256 * rSges(5,1) + t258 * t287;
t172 = t256 * rSges(4,3) + t258 * t288;
t171 = -rSges(4,3) * t258 + t256 * t288;
t155 = pkin(8) * t241 - t243 * t325;
t154 = V_base(5) * rSges(2,3) - t231 * t244 + t297;
t153 = t232 * t244 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = rSges(6,3) * t241 + (-rSges(6,1) * t251 - rSges(6,2) * t252) * t243;
t151 = Icges(6,5) * t241 + (-Icges(6,1) * t251 - Icges(6,4) * t252) * t243;
t150 = Icges(6,6) * t241 + (-Icges(6,4) * t251 - Icges(6,2) * t252) * t243;
t149 = Icges(6,3) * t241 + (-Icges(6,5) * t251 - Icges(6,6) * t252) * t243;
t147 = t231 * V_base(4) - t232 * V_base(5) + V_base(3);
t145 = rSges(7,3) * t241 + (-rSges(7,1) * t240 - rSges(7,2) * t242) * t243;
t144 = Icges(7,5) * t241 + (-Icges(7,1) * t240 - Icges(7,4) * t242) * t243;
t143 = Icges(7,6) * t241 + (-Icges(7,4) * t240 - Icges(7,2) * t242) * t243;
t142 = Icges(7,3) * t241 + (-Icges(7,5) * t240 - Icges(7,6) * t242) * t243;
t140 = t256 * t271 - t258 * t323;
t139 = t256 * t323 + t258 * t271;
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
t122 = t230 * t235 + (-t190 - t233) * t244 + t297;
t121 = t191 * t244 - t230 * t236 + t275;
t120 = t190 * t236 - t191 * t235 + t274;
t119 = t210 * t235 + (-t171 + t302) * t244 + t290;
t118 = t244 * t172 + (-t210 - t326) * t236 + t266;
t117 = t171 * t236 + (-t158 - t172) * t235 + t272;
t116 = t209 * t235 + (-t174 + t294) * t244 + t273;
t115 = t244 * t173 + (-t209 + t292) * t236 + t264;
t114 = t174 * t236 + (-t173 + t301) * t235 + t265;
t113 = t152 * t235 + (-t138 + t291) * t244 + t267;
t112 = t244 * t137 + (-t152 + t276) * t236 + t262;
t111 = t138 * t236 + (-t137 + t293) * t235 + t263;
t110 = -t130 * t213 + t145 * t196 + t155 * t235 + (-t140 + t291) * t244 + t267;
t109 = t213 * t129 + t244 * t139 - t197 * t145 + (-t155 + t276) * t236 + t262;
t108 = -t129 * t196 + t130 * t197 + t140 * t236 + (-t139 + t293) * t235 + t263;
t1 = m(2) * (t147 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + t196 * ((t123 * t311 + t125 * t179 + t127 * t180) * t197 + (t124 * t311 + t179 * t126 + t180 * t128) * t196 + (t142 * t311 + t143 * t179 + t144 * t180) * t213) / 0.2e1 + t197 * ((t123 * t310 + t177 * t125 + t178 * t127) * t197 + (t124 * t310 + t177 * t126 + t178 * t128) * t196 + (t142 * t310 + t177 * t143 + t178 * t144) * t213) / 0.2e1 + t213 * ((t123 * t197 + t124 * t196 + t142 * t213) * t241 + ((-t125 * t242 - t127 * t240) * t197 + (-t126 * t242 - t128 * t240) * t196 + (-t143 * t242 - t144 * t240) * t213) * t243) / 0.2e1 + m(1) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + ((-t256 * t221 + t224 * t258 + Icges(1,4)) * V_base(5) + (-t256 * t222 + t258 * t225 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t258 * t221 + t256 * t224 + Icges(1,2)) * V_base(5) + (t222 * t258 + t256 * t225 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t131 * t311 + t133 * t194 + t135 * t195) * t236 + (t132 * t311 + t194 * t134 + t195 * t136) * t235 + (t149 * t311 + t150 * t194 + t151 * t195) * t244 - t330 * t258 + t331 * t256) * t235 / 0.2e1 + ((t131 * t310 + t192 * t133 + t193 * t135) * t236 + (t132 * t310 + t192 * t134 + t193 * t136) * t235 + (t149 * t310 + t192 * t150 + t193 * t151) * t244 + t331 * t258 + t330 * t256) * t236 / 0.2e1 + ((t185 * t257 + t187 * t255 + (-t133 * t252 - t135 * t251 - t338) * t243 + (t131 + t337) * t241) * t236 + (t184 * t257 + t186 * t255 + (-t134 * t252 - t136 * t251 - t339) * t243 + (t132 + t336) * t241) * t235 + (t257 * t220 + t255 * t223 + Icges(2,3) + (-t150 * t252 - t151 * t251 - t335) * t243 + (t149 + t334) * t241) * t244) * t244 / 0.2e1 + t244 * V_base(4) * (Icges(2,5) * t258 - Icges(2,6) * t256) + t244 * V_base(5) * (Icges(2,5) * t256 + Icges(2,6) * t258) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
