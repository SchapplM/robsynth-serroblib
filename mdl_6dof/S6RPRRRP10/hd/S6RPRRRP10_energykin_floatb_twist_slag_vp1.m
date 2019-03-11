% Calculate kinetic energy for
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:13
% EndTime: 2019-03-09 06:30:16
% DurationCPUTime: 3.44s
% Computational Cost: add. (1434->306), mult. (2106->435), div. (0->0), fcn. (2010->8), ass. (0->152)
t325 = Icges(2,4) + Icges(3,6);
t324 = Icges(2,1) + Icges(3,2);
t323 = Icges(6,1) + Icges(7,1);
t322 = -Icges(3,4) + Icges(2,5);
t321 = Icges(6,4) - Icges(7,5);
t320 = Icges(7,4) + Icges(6,5);
t319 = Icges(3,5) - Icges(2,6);
t318 = Icges(2,2) + Icges(3,3);
t317 = Icges(6,2) + Icges(7,3);
t316 = Icges(7,6) - Icges(6,6);
t315 = Icges(6,3) + Icges(7,2);
t314 = rSges(7,1) + pkin(5);
t313 = rSges(7,3) + qJ(6);
t241 = cos(qJ(1));
t312 = t325 * t241;
t238 = sin(qJ(1));
t311 = t325 * t238;
t235 = qJ(4) + qJ(5);
t230 = sin(t235);
t231 = cos(t235);
t272 = t241 * t231;
t237 = sin(qJ(3));
t278 = t237 * t238;
t174 = t230 * t278 - t272;
t175 = t230 * t241 + t231 * t278;
t240 = cos(qJ(3));
t275 = t238 * t240;
t310 = t174 * t317 - t175 * t321 - t275 * t316;
t277 = t237 * t241;
t176 = t230 * t277 + t231 * t238;
t177 = t230 * t238 - t237 * t272;
t273 = t240 * t241;
t309 = -t176 * t317 - t177 * t321 + t273 * t316;
t308 = t174 * t316 + t175 * t320 - t275 * t315;
t307 = -t176 * t316 + t177 * t320 + t273 * t315;
t306 = -t174 * t321 + t175 * t323 - t275 * t320;
t305 = t176 * t321 + t177 * t323 + t273 * t320;
t304 = (t230 * t317 - t231 * t321) * t240 + t316 * t237;
t303 = (t230 * t316 + t231 * t320) * t240 + t315 * t237;
t302 = (-t230 * t321 + t231 * t323) * t240 + t320 * t237;
t301 = -t241 * t318 - t311;
t300 = t238 * t318 - t312;
t299 = t238 * t324 + t312;
t298 = t241 * t324 - t311;
t239 = cos(qJ(4));
t287 = pkin(4) * t239;
t295 = -pkin(9) * t240 + t237 * t287;
t284 = Icges(4,4) * t237;
t257 = Icges(4,2) * t240 + t284;
t166 = Icges(4,6) * t241 + t238 * t257;
t167 = Icges(4,6) * t238 - t241 * t257;
t283 = Icges(4,4) * t240;
t258 = Icges(4,1) * t237 + t283;
t169 = Icges(4,5) * t241 + t238 * t258;
t170 = Icges(4,5) * t238 - t241 * t258;
t201 = -Icges(4,2) * t237 + t283;
t206 = Icges(4,1) * t240 - t284;
t219 = qJD(3) * t238 + V_base(5);
t220 = qJD(3) * t241 + V_base(4);
t225 = V_base(6) + qJD(1);
t294 = (t166 * t240 + t169 * t237) * t220 + (t167 * t240 + t170 * t237) * t219 + (t201 * t240 + t206 * t237) * t225;
t289 = pkin(7) * t238;
t288 = pkin(7) * t241;
t236 = sin(qJ(4));
t280 = t236 * t238;
t279 = t236 * t241;
t276 = t238 * t239;
t274 = t239 * t241;
t271 = -rSges(7,2) * t275 + t313 * t174 + t314 * t175;
t270 = rSges(7,2) * t273 - t313 * t176 + t314 * t177;
t269 = rSges(7,2) * t237 + (t313 * t230 + t314 * t231) * t240;
t268 = qJD(4) * t240;
t210 = pkin(1) * t238 - qJ(2) * t241;
t267 = V_base(4) * t210 + V_base(3);
t266 = V_base(5) * pkin(6) + V_base(1);
t263 = -t210 - t289;
t262 = qJD(2) * t238 + t266;
t180 = t241 * t268 + t219;
t209 = qJD(4) * t237 + t225;
t261 = V_base(5) * pkin(2) + t262;
t260 = pkin(3) * t237 - pkin(8) * t240;
t259 = rSges(4,1) * t237 + rSges(4,2) * t240;
t256 = Icges(4,5) * t237 + Icges(4,6) * t240;
t214 = pkin(1) * t241 + qJ(2) * t238;
t252 = -qJD(2) * t241 + t225 * t214 + V_base(2);
t251 = (Icges(4,3) * t241 + t238 * t256) * t220 + (Icges(4,3) * t238 - t241 * t256) * t219 + (Icges(4,5) * t240 - Icges(4,6) * t237) * t225;
t250 = V_base(4) * t289 + (-t214 - t288) * V_base(5) + t267;
t249 = t225 * t288 + (-pkin(2) - pkin(6)) * V_base(4) + t252;
t188 = t260 * t241;
t217 = t240 * pkin(3) + t237 * pkin(8);
t248 = t219 * t217 + (t188 + t263) * t225 + t261;
t187 = t260 * t238;
t247 = -t187 * t219 - t220 * t188 + t250;
t246 = t225 * t187 - t217 * t220 + t249;
t144 = pkin(4) * t280 - t241 * t295;
t149 = pkin(9) * t237 + t240 * t287;
t245 = -t144 * t209 + t180 * t149 + t248;
t143 = pkin(4) * t279 + t238 * t295;
t181 = -t238 * t268 + t220;
t244 = -t143 * t180 + t181 * t144 + t247;
t243 = t209 * t143 - t149 * t181 + t246;
t216 = rSges(2,1) * t241 - rSges(2,2) * t238;
t215 = -rSges(3,2) * t241 + rSges(3,3) * t238;
t213 = rSges(4,1) * t240 - rSges(4,2) * t237;
t212 = rSges(2,1) * t238 + rSges(2,2) * t241;
t211 = -rSges(3,2) * t238 - rSges(3,3) * t241;
t193 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t192 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t191 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t189 = qJD(5) * t237 + t209;
t185 = -t237 * t274 + t280;
t184 = t236 * t277 + t276;
t183 = t237 * t276 + t279;
t182 = -t236 * t278 + t274;
t173 = rSges(4,3) * t238 - t241 * t259;
t172 = rSges(5,3) * t237 + (rSges(5,1) * t239 - rSges(5,2) * t236) * t240;
t171 = rSges(4,3) * t241 + t238 * t259;
t168 = Icges(5,5) * t237 + (Icges(5,1) * t239 - Icges(5,4) * t236) * t240;
t165 = Icges(5,6) * t237 + (Icges(5,4) * t239 - Icges(5,2) * t236) * t240;
t162 = Icges(5,3) * t237 + (Icges(5,5) * t239 - Icges(5,6) * t236) * t240;
t160 = (-qJD(4) - qJD(5)) * t275 + t220;
t159 = qJD(5) * t273 + t180;
t158 = rSges(6,3) * t237 + (rSges(6,1) * t231 - rSges(6,2) * t230) * t240;
t148 = V_base(5) * rSges(2,3) - t212 * t225 + t266;
t147 = t216 * t225 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t146 = t212 * V_base(4) - t216 * V_base(5) + V_base(3);
t142 = rSges(5,1) * t185 + rSges(5,2) * t184 + rSges(5,3) * t273;
t141 = rSges(5,1) * t183 + rSges(5,2) * t182 - rSges(5,3) * t275;
t140 = Icges(5,1) * t185 + Icges(5,4) * t184 + Icges(5,5) * t273;
t139 = Icges(5,1) * t183 + Icges(5,4) * t182 - Icges(5,5) * t275;
t138 = Icges(5,4) * t185 + Icges(5,2) * t184 + Icges(5,6) * t273;
t137 = Icges(5,4) * t183 + Icges(5,2) * t182 - Icges(5,6) * t275;
t136 = Icges(5,5) * t185 + Icges(5,6) * t184 + Icges(5,3) * t273;
t135 = Icges(5,5) * t183 + Icges(5,6) * t182 - Icges(5,3) * t275;
t132 = V_base(5) * rSges(3,1) + (-t210 - t211) * t225 + t262;
t131 = t215 * t225 + (-rSges(3,1) - pkin(6)) * V_base(4) + t252;
t130 = rSges(6,1) * t177 + rSges(6,2) * t176 + rSges(6,3) * t273;
t128 = rSges(6,1) * t175 - rSges(6,2) * t174 - rSges(6,3) * t275;
t113 = t211 * V_base(4) + (-t214 - t215) * V_base(5) + t267;
t111 = t213 * t219 + (-t173 + t263) * t225 + t261;
t110 = t171 * t225 - t213 * t220 + t249;
t109 = -t171 * t219 + t173 * t220 + t250;
t108 = -t142 * t209 + t172 * t180 + t248;
t107 = t141 * t209 - t172 * t181 + t246;
t106 = -t141 * t180 + t142 * t181 + t247;
t105 = -t130 * t189 + t158 * t159 + t245;
t104 = t128 * t189 - t158 * t160 + t243;
t103 = -t128 * t159 + t130 * t160 + t244;
t102 = qJD(6) * t174 + t159 * t269 - t189 * t270 + t245;
t101 = -qJD(6) * t176 - t160 * t269 + t189 * t271 + t243;
t100 = qJD(6) * t230 * t240 - t159 * t271 + t160 * t270 + t244;
t1 = t181 * ((-t135 * t275 + t182 * t137 + t183 * t139) * t181 + (-t136 * t275 + t138 * t182 + t140 * t183) * t180 + (-t162 * t275 + t165 * t182 + t168 * t183) * t209) / 0.2e1 + t180 * ((t135 * t273 + t137 * t184 + t139 * t185) * t181 + (t136 * t273 + t184 * t138 + t185 * t140) * t180 + (t162 * t273 + t165 * t184 + t185 * t168) * t209) / 0.2e1 + t220 * (t238 * t294 + t251 * t241) / 0.2e1 + t219 * (t251 * t238 - t241 * t294) / 0.2e1 + t209 * ((t135 * t181 + t136 * t180 + t162 * t209) * t237 + ((-t137 * t236 + t139 * t239) * t181 + (-t138 * t236 + t140 * t239) * t180 + (-t165 * t236 + t168 * t239) * t209) * t240) / 0.2e1 + m(1) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(2) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(3) * (t113 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + ((-t304 * t176 + t302 * t177 + t303 * t273) * t189 + (-t310 * t176 + t306 * t177 + t308 * t273) * t160 + (-t309 * t176 + t305 * t177 + t307 * t273) * t159) * t159 / 0.2e1 + ((t304 * t174 + t302 * t175 - t303 * t275) * t189 + (t310 * t174 + t306 * t175 - t308 * t275) * t160 + (t309 * t174 + t305 * t175 - t307 * t275) * t159) * t160 / 0.2e1 + (((t304 * t230 + t302 * t231) * t189 + (t310 * t230 + t306 * t231) * t160 + (t309 * t230 + t305 * t231) * t159) * t240 + (t307 * t159 + t308 * t160 + t303 * t189) * t237) * t189 / 0.2e1 + ((-t166 * t237 + t169 * t240) * t220 + (-t167 * t237 + t170 * t240) * t219 + (-t237 * t201 + t240 * t206 + Icges(3,1) + Icges(2,3)) * t225) * t225 / 0.2e1 + ((t301 * t238 + t299 * t241 + Icges(1,4)) * V_base(5) + (t300 * t238 + t298 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t299 * t238 - t301 * t241 + Icges(1,2)) * V_base(5) + (t298 * t238 - t300 * t241 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t225 * (t238 * t322 - t319 * t241) + V_base(4) * t225 * (t238 * t319 + t322 * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
