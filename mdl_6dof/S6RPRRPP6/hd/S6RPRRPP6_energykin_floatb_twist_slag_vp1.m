% Calculate kinetic energy for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:17
% EndTime: 2019-03-09 04:46:20
% DurationCPUTime: 3.17s
% Computational Cost: add. (1398->299), mult. (2046->402), div. (0->0), fcn. (1950->8), ass. (0->148)
t326 = Icges(2,4) + Icges(3,6);
t325 = Icges(2,1) + Icges(3,2);
t324 = Icges(6,1) + Icges(7,1);
t323 = -Icges(3,4) + Icges(2,5);
t322 = Icges(6,4) - Icges(7,5);
t321 = Icges(7,4) + Icges(6,5);
t320 = Icges(3,5) - Icges(2,6);
t319 = Icges(2,2) + Icges(3,3);
t318 = Icges(6,2) + Icges(7,3);
t317 = Icges(7,6) - Icges(6,6);
t316 = Icges(6,3) + Icges(7,2) + Icges(5,3);
t315 = rSges(7,1) + pkin(5);
t314 = rSges(7,3) + qJ(6);
t241 = cos(qJ(1));
t313 = t326 * t241;
t238 = sin(qJ(1));
t312 = t326 * t238;
t237 = sin(qJ(3));
t234 = qJ(4) + pkin(9);
t225 = cos(t234);
t272 = t241 * t225;
t224 = sin(t234);
t279 = t238 * t224;
t161 = t237 * t279 - t272;
t278 = t238 * t225;
t162 = t224 * t241 + t237 * t278;
t240 = cos(qJ(3));
t275 = t238 * t240;
t311 = t318 * t161 - t322 * t162 - t317 * t275;
t280 = t237 * t241;
t163 = t224 * t280 + t278;
t164 = -t237 * t272 + t279;
t273 = t240 * t241;
t310 = -t318 * t163 - t322 * t164 + t317 * t273;
t309 = -t322 * t161 + t324 * t162 - t321 * t275;
t308 = t322 * t163 + t324 * t164 + t321 * t273;
t307 = (t318 * t224 - t322 * t225) * t240 + t317 * t237;
t306 = (-t322 * t224 + t324 * t225) * t240 + t321 * t237;
t305 = -t319 * t241 - t312;
t304 = t319 * t238 - t313;
t303 = t325 * t238 + t313;
t302 = t325 * t241 - t312;
t239 = cos(qJ(4));
t274 = t239 * t241;
t236 = sin(qJ(4));
t277 = t238 * t236;
t182 = -t237 * t277 + t274;
t276 = t238 * t239;
t281 = t236 * t241;
t183 = t237 * t276 + t281;
t299 = Icges(5,5) * t183 + Icges(5,6) * t182 + t317 * t161 + t321 * t162 - t316 * t275;
t184 = t236 * t280 + t276;
t185 = -t237 * t274 + t277;
t298 = Icges(5,5) * t185 + Icges(5,6) * t184 - t317 * t163 + t321 * t164 + t316 * t273;
t297 = (Icges(5,5) * t239 - Icges(5,6) * t236 + t317 * t224 + t321 * t225) * t240 + t316 * t237;
t288 = pkin(4) * t239;
t296 = -qJ(5) * t240 + t237 * t288;
t285 = Icges(4,4) * t237;
t256 = Icges(4,2) * t240 + t285;
t170 = Icges(4,6) * t241 + t238 * t256;
t171 = Icges(4,6) * t238 - t241 * t256;
t284 = Icges(4,4) * t240;
t257 = Icges(4,1) * t237 + t284;
t173 = Icges(4,5) * t241 + t238 * t257;
t174 = Icges(4,5) * t238 - t241 * t257;
t200 = -Icges(4,2) * t237 + t284;
t205 = Icges(4,1) * t240 - t285;
t219 = qJD(3) * t238 + V_base(5);
t220 = qJD(3) * t241 + V_base(4);
t226 = V_base(6) + qJD(1);
t295 = (t170 * t240 + t173 * t237) * t220 + (t171 * t240 + t174 * t237) * t219 + (t200 * t240 + t205 * t237) * t226;
t290 = pkin(7) * t238;
t289 = pkin(7) * t241;
t271 = -rSges(7,2) * t275 + t314 * t161 + t162 * t315;
t270 = rSges(7,2) * t273 - t314 * t163 + t164 * t315;
t269 = rSges(7,2) * t237 + (t314 * t224 + t225 * t315) * t240;
t268 = qJD(4) * t240;
t267 = qJD(5) * t240;
t209 = t238 * pkin(1) - qJ(2) * t241;
t266 = V_base(4) * t209 + V_base(3);
t265 = V_base(5) * pkin(6) + V_base(1);
t262 = -t209 - t290;
t261 = qJD(2) * t238 + t265;
t260 = V_base(5) * pkin(2) + t261;
t259 = pkin(3) * t237 - pkin(8) * t240;
t258 = rSges(4,1) * t237 + rSges(4,2) * t240;
t255 = Icges(4,5) * t237 + Icges(4,6) * t240;
t213 = pkin(1) * t241 + t238 * qJ(2);
t251 = -qJD(2) * t241 + t226 * t213 + V_base(2);
t250 = (Icges(4,3) * t241 + t238 * t255) * t220 + (Icges(4,3) * t238 - t241 * t255) * t219 + (Icges(4,5) * t240 - Icges(4,6) * t237) * t226;
t249 = V_base(4) * t290 + (-t213 - t289) * V_base(5) + t266;
t248 = t226 * t289 + (-pkin(2) - pkin(6)) * V_base(4) + t251;
t188 = t259 * t241;
t216 = pkin(3) * t240 + pkin(8) * t237;
t247 = t219 * t216 + (t188 + t262) * t226 + t260;
t187 = t259 * t238;
t246 = -t219 * t187 - t220 * t188 + t249;
t245 = t226 * t187 - t220 * t216 + t248;
t144 = pkin(4) * t277 - t241 * t296;
t181 = -t238 * t268 + t220;
t244 = qJD(5) * t237 + t181 * t144 + t246;
t143 = pkin(4) * t281 + t238 * t296;
t208 = qJD(4) * t237 + t226;
t243 = t208 * t143 + t241 * t267 + t245;
t151 = qJ(5) * t237 + t240 * t288;
t180 = t241 * t268 + t219;
t242 = t180 * t151 - t238 * t267 + t247;
t215 = rSges(2,1) * t241 - t238 * rSges(2,2);
t214 = -rSges(3,2) * t241 + t238 * rSges(3,3);
t212 = rSges(4,1) * t240 - rSges(4,2) * t237;
t211 = t238 * rSges(2,1) + rSges(2,2) * t241;
t210 = -t238 * rSges(3,2) - rSges(3,3) * t241;
t192 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t191 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t190 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t178 = t238 * rSges(4,3) - t241 * t258;
t177 = rSges(5,3) * t237 + (rSges(5,1) * t239 - rSges(5,2) * t236) * t240;
t176 = rSges(4,3) * t241 + t238 * t258;
t172 = Icges(5,5) * t237 + (Icges(5,1) * t239 - Icges(5,4) * t236) * t240;
t169 = Icges(5,6) * t237 + (Icges(5,4) * t239 - Icges(5,2) * t236) * t240;
t159 = rSges(6,3) * t237 + (rSges(6,1) * t225 - rSges(6,2) * t224) * t240;
t150 = V_base(5) * rSges(2,3) - t211 * t226 + t265;
t149 = t215 * t226 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t148 = t211 * V_base(4) - t215 * V_base(5) + V_base(3);
t146 = t185 * rSges(5,1) + t184 * rSges(5,2) + rSges(5,3) * t273;
t145 = rSges(5,1) * t183 + rSges(5,2) * t182 - rSges(5,3) * t275;
t142 = Icges(5,1) * t185 + Icges(5,4) * t184 + Icges(5,5) * t273;
t141 = Icges(5,1) * t183 + Icges(5,4) * t182 - Icges(5,5) * t275;
t140 = Icges(5,4) * t185 + Icges(5,2) * t184 + Icges(5,6) * t273;
t139 = Icges(5,4) * t183 + Icges(5,2) * t182 - Icges(5,6) * t275;
t134 = V_base(5) * rSges(3,1) + (-t209 - t210) * t226 + t261;
t133 = t226 * t214 + (-rSges(3,1) - pkin(6)) * V_base(4) + t251;
t132 = t164 * rSges(6,1) + t163 * rSges(6,2) + rSges(6,3) * t273;
t130 = rSges(6,1) * t162 - rSges(6,2) * t161 - rSges(6,3) * t275;
t115 = t210 * V_base(4) + (-t213 - t214) * V_base(5) + t266;
t113 = t212 * t219 + (-t178 + t262) * t226 + t260;
t112 = t226 * t176 - t220 * t212 + t248;
t111 = -t219 * t176 + t220 * t178 + t249;
t110 = -t146 * t208 + t177 * t180 + t247;
t109 = t208 * t145 - t181 * t177 + t245;
t108 = -t180 * t145 + t181 * t146 + t246;
t107 = t159 * t180 + (-t132 - t144) * t208 + t242;
t106 = t208 * t130 + (-t151 - t159) * t181 + t243;
t105 = t181 * t132 + (-t130 - t143) * t180 + t244;
t104 = qJD(6) * t161 + t269 * t180 + (-t144 - t270) * t208 + t242;
t103 = -qJD(6) * t163 + t271 * t208 + (-t151 - t269) * t181 + t243;
t102 = qJD(6) * t240 * t224 + t270 * t181 + (-t143 - t271) * t180 + t244;
t1 = t220 * (t295 * t238 + t250 * t241) / 0.2e1 + t219 * (t250 * t238 - t295 * t241) / 0.2e1 + m(1) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(2) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(6) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + ((-t163 * t307 + t164 * t306 + t184 * t169 + t185 * t172 + t273 * t297) * t208 + (t184 * t139 + t185 * t141 - t163 * t311 + t309 * t164 + t299 * t273) * t181 + (t184 * t140 + t185 * t142 - t310 * t163 + t308 * t164 + t298 * t273) * t180) * t180 / 0.2e1 + ((t161 * t307 + t162 * t306 + t169 * t182 + t172 * t183 - t275 * t297) * t208 + (t182 * t139 + t183 * t141 + t311 * t161 + t309 * t162 - t299 * t275) * t181 + (t140 * t182 + t142 * t183 + t161 * t310 + t162 * t308 - t275 * t298) * t180) * t181 / 0.2e1 + (((-t169 * t236 + t172 * t239 + t224 * t307 + t225 * t306) * t208 + (-t139 * t236 + t141 * t239 + t224 * t311 + t309 * t225) * t181 + (-t140 * t236 + t142 * t239 + t224 * t310 + t225 * t308) * t180) * t240 + (t180 * t298 + t181 * t299 + t208 * t297) * t237) * t208 / 0.2e1 + ((-t170 * t237 + t173 * t240) * t220 + (-t171 * t237 + t174 * t240) * t219 + (-t237 * t200 + t240 * t205 + Icges(3,1) + Icges(2,3)) * t226) * t226 / 0.2e1 + ((t238 * t305 + t241 * t303 + Icges(1,4)) * V_base(5) + (t238 * t304 + t241 * t302 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t238 * t303 - t241 * t305 + Icges(1,2)) * V_base(5) + (t238 * t302 - t241 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t226 * (t320 * t238 + t323 * t241) + V_base(5) * t226 * (t323 * t238 - t320 * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
