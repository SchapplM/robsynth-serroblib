% Calculate kinetic energy for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:46
% EndTime: 2019-03-09 08:13:50
% DurationCPUTime: 3.65s
% Computational Cost: add. (1162->313), mult. (1904->432), div. (0->0), fcn. (1742->8), ass. (0->156)
t332 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t331 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t330 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t244 = sin(qJ(2));
t329 = t332 * t244;
t246 = cos(qJ(2));
t328 = t332 * t246;
t327 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t326 = Icges(5,5) + Icges(3,6) - Icges(4,6);
t325 = t330 * t244 - t328;
t324 = t331 * t246 - t329;
t323 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t245 = sin(qJ(1));
t247 = cos(qJ(1));
t322 = t325 * t245 + t326 * t247;
t321 = -t326 * t245 + t325 * t247;
t320 = t324 * t245 - t327 * t247;
t319 = t327 * t245 + t324 * t247;
t318 = -t330 * t246 - t329;
t317 = t331 * t244 + t328;
t316 = -t326 * t244 + t327 * t246;
t225 = -qJD(2) * t247 + V_base(5);
t226 = qJD(2) * t245 + V_base(4);
t234 = V_base(6) + qJD(1);
t315 = (t318 * t244 + t317 * t246) * t234 + (t321 * t244 + t319 * t246) * t226 + (t322 * t244 + t320 * t246) * t225;
t314 = (t327 * t244 + t326 * t246) * t234 + (t323 * t245 + t316 * t247) * t226 + (t316 * t245 - t323 * t247) * t225;
t310 = pkin(3) * t244;
t242 = cos(pkin(9));
t309 = pkin(5) * t242;
t308 = Icges(2,4) * t245;
t241 = sin(pkin(9));
t301 = t241 * t247;
t300 = t244 * t247;
t240 = pkin(9) + qJ(6);
t232 = sin(t240);
t299 = t245 * t232;
t233 = cos(t240);
t298 = t245 * t233;
t297 = t245 * t241;
t296 = t245 * t242;
t295 = t245 * t246;
t294 = t246 * t247;
t275 = pkin(2) * t246 + qJ(3) * t244;
t189 = t275 * t245;
t223 = t245 * pkin(1) - pkin(7) * t247;
t292 = -t189 - t223;
t191 = t275 * t247;
t196 = pkin(3) * t294 - t245 * qJ(4);
t291 = -t191 - t196;
t290 = qJD(3) * t244;
t289 = qJD(5) * t246;
t288 = qJD(6) * t246;
t287 = V_base(5) * pkin(6) + V_base(1);
t195 = pkin(3) * t295 + qJ(4) * t247;
t284 = -t195 + t292;
t274 = pkin(4) * t244 + qJ(5) * t246;
t190 = t274 * t247;
t283 = -t190 + t291;
t216 = pkin(2) * t244 - qJ(3) * t246;
t282 = -t216 - t310;
t188 = t274 * t245;
t281 = -t188 + t284;
t220 = -pkin(4) * t246 + qJ(5) * t244;
t280 = -t220 + t282;
t279 = t225 * t216 + t247 * t290 + t287;
t278 = rSges(3,1) * t246 - rSges(3,2) * t244;
t277 = rSges(4,1) * t246 + rSges(4,3) * t244;
t276 = rSges(5,1) * t244 - rSges(5,2) * t246;
t224 = pkin(1) * t247 + t245 * pkin(7);
t264 = -V_base(4) * pkin(6) + t234 * t224 + V_base(2);
t263 = V_base(4) * t223 - t224 * V_base(5) + V_base(3);
t259 = pkin(8) * t246 + t244 * t309;
t258 = t234 * t191 + t245 * t290 + t264;
t257 = -qJD(4) * t245 + t225 * t310 + t279;
t256 = -qJD(3) * t246 + t226 * t189 + t263;
t255 = qJD(4) * t247 + t234 * t196 + t258;
t254 = t225 * t220 + t247 * t289 + t257;
t253 = t226 * t195 + t256;
t252 = t234 * t190 + t245 * t289 + t255;
t251 = qJD(5) * t244 + t226 * t188 + t253;
t238 = Icges(2,4) * t247;
t222 = rSges(2,1) * t247 - t245 * rSges(2,2);
t221 = -rSges(5,1) * t246 - rSges(5,2) * t244;
t219 = t245 * rSges(2,1) + rSges(2,2) * t247;
t218 = rSges(3,1) * t244 + rSges(3,2) * t246;
t217 = rSges(4,1) * t244 - rSges(4,3) * t246;
t215 = qJD(6) * t244 + t234;
t214 = Icges(2,1) * t247 - t308;
t213 = Icges(2,1) * t245 + t238;
t209 = -Icges(2,2) * t245 + t238;
t208 = Icges(2,2) * t247 + t308;
t199 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t198 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t197 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t187 = t247 * t288 + t226;
t186 = t245 * t288 + t225;
t185 = t242 * t300 - t297;
t184 = -t241 * t300 - t296;
t183 = t244 * t296 + t301;
t182 = t242 * t247 - t244 * t297;
t178 = t245 * rSges(3,3) + t247 * t278;
t177 = t245 * rSges(4,2) + t247 * t277;
t176 = -t245 * rSges(5,3) + t247 * t276;
t175 = -rSges(3,3) * t247 + t245 * t278;
t174 = -rSges(4,2) * t247 + t245 * t277;
t173 = rSges(5,3) * t247 + t245 * t276;
t154 = t233 * t300 - t299;
t153 = -t232 * t300 - t298;
t152 = t232 * t247 + t244 * t298;
t151 = t233 * t247 - t244 * t299;
t150 = rSges(6,3) * t244 + (-rSges(6,1) * t242 + rSges(6,2) * t241) * t246;
t147 = Icges(6,5) * t244 + (-Icges(6,1) * t242 + Icges(6,4) * t241) * t246;
t146 = Icges(6,6) * t244 + (-Icges(6,4) * t242 + Icges(6,2) * t241) * t246;
t145 = Icges(6,3) * t244 + (-Icges(6,5) * t242 + Icges(6,6) * t241) * t246;
t141 = rSges(7,3) * t244 + (-rSges(7,1) * t233 + rSges(7,2) * t232) * t246;
t140 = Icges(7,5) * t244 + (-Icges(7,1) * t233 + Icges(7,4) * t232) * t246;
t139 = Icges(7,6) * t244 + (-Icges(7,4) * t233 + Icges(7,2) * t232) * t246;
t138 = Icges(7,3) * t244 + (-Icges(7,5) * t233 + Icges(7,6) * t232) * t246;
t137 = pkin(8) * t244 - t246 * t309;
t136 = V_base(5) * rSges(2,3) - t219 * t234 + t287;
t135 = t222 * t234 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t134 = t219 * V_base(4) - t222 * V_base(5) + V_base(3);
t133 = -pkin(5) * t297 + t247 * t259;
t132 = pkin(5) * t301 + t245 * t259;
t131 = t185 * rSges(6,1) + t184 * rSges(6,2) + rSges(6,3) * t294;
t130 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t295;
t129 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t294;
t128 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t295;
t127 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t294;
t126 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t295;
t125 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t294;
t124 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t295;
t123 = t154 * rSges(7,1) + t153 * rSges(7,2) + rSges(7,3) * t294;
t122 = rSges(7,1) * t152 + rSges(7,2) * t151 + rSges(7,3) * t295;
t121 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t294;
t120 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t295;
t119 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t294;
t118 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t295;
t117 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t294;
t116 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t295;
t115 = t218 * t225 + (-t175 - t223) * t234 + t287;
t114 = t178 * t234 - t218 * t226 + t264;
t113 = t175 * t226 - t178 * t225 + t263;
t112 = t217 * t225 + (-t174 + t292) * t234 + t279;
t111 = t177 * t234 + (-t216 - t217) * t226 + t258;
t110 = t174 * t226 + (-t177 - t191) * t225 + t256;
t109 = t221 * t225 + (-t173 + t284) * t234 + t257;
t108 = t176 * t234 + (-t221 + t282) * t226 + t255;
t107 = t173 * t226 + (-t176 + t291) * t225 + t253;
t106 = t150 * t225 + (-t130 + t281) * t234 + t254;
t105 = t131 * t234 + (-t150 + t280) * t226 + t252;
t104 = t130 * t226 + (-t131 + t283) * t225 + t251;
t103 = -t122 * t215 + t137 * t225 + t141 * t186 + (-t132 + t281) * t234 + t254;
t102 = t123 * t215 + t133 * t234 - t141 * t187 + (-t137 + t280) * t226 + t252;
t101 = t122 * t187 - t123 * t186 + t132 * t226 + (-t133 + t283) * t225 + t251;
t1 = t186 * ((t117 * t295 + t119 * t151 + t121 * t152) * t187 + (t116 * t295 + t151 * t118 + t152 * t120) * t186 + (t138 * t295 + t139 * t151 + t140 * t152) * t215) / 0.2e1 + t187 * ((t117 * t294 + t153 * t119 + t154 * t121) * t187 + (t116 * t294 + t153 * t118 + t154 * t120) * t186 + (t138 * t294 + t153 * t139 + t154 * t140) * t215) / 0.2e1 + m(1) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + m(2) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t215 * ((t116 * t186 + t117 * t187 + t138 * t215) * t244 + ((t119 * t232 - t121 * t233) * t187 + (t118 * t232 - t120 * t233) * t186 + (t139 * t232 - t140 * t233) * t215) * t246) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(3) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + ((-t245 * t208 + t213 * t247 + Icges(1,4)) * V_base(5) + (-t245 * t209 + t247 * t214 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t247 * t208 + t245 * t213 + Icges(1,2)) * V_base(5) + (t209 * t247 + t245 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t125 * t295 + t127 * t182 + t129 * t183) * t226 + (t124 * t295 + t182 * t126 + t183 * t128) * t225 + (t145 * t295 + t146 * t182 + t147 * t183) * t234 - t314 * t247 + t315 * t245) * t225 / 0.2e1 + ((t125 * t294 + t184 * t127 + t185 * t129) * t226 + (t124 * t294 + t184 * t126 + t185 * t128) * t225 + (t145 * t294 + t184 * t146 + t185 * t147) * t234 + t315 * t247 + t314 * t245) * t226 / 0.2e1 + (((t127 * t241 - t129 * t242 - t321) * t246 + (t125 + t319) * t244) * t226 + ((t126 * t241 - t128 * t242 - t322) * t246 + (t124 + t320) * t244) * t225 + (Icges(2,3) + (t146 * t241 - t147 * t242 - t318) * t246 + (t145 + t317) * t244) * t234) * t234 / 0.2e1 + V_base(4) * t234 * (Icges(2,5) * t247 - Icges(2,6) * t245) + V_base(5) * t234 * (Icges(2,5) * t245 + Icges(2,6) * t247) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
