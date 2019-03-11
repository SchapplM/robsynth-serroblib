% Calculate kinetic energy for
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:33
% EndTime: 2019-03-09 04:42:35
% DurationCPUTime: 2.96s
% Computational Cost: add. (1785->287), mult. (2201->399), div. (0->0), fcn. (2119->8), ass. (0->146)
t313 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t312 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t311 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t310 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t309 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t308 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t307 = rSges(7,1) + pkin(5);
t306 = rSges(7,3) + qJ(6);
t237 = pkin(9) + qJ(3);
t231 = cos(t237);
t243 = cos(qJ(4));
t244 = cos(qJ(1));
t280 = t243 * t244;
t241 = sin(qJ(4));
t242 = sin(qJ(1));
t282 = t242 * t241;
t194 = t231 * t282 + t280;
t281 = t242 * t243;
t283 = t241 * t244;
t195 = t231 * t281 - t283;
t230 = sin(t237);
t285 = t230 * t242;
t305 = t309 * t194 + t311 * t195 - t308 * t285;
t196 = t231 * t283 - t281;
t197 = t231 * t280 + t282;
t284 = t230 * t244;
t304 = t309 * t196 + t311 * t197 - t308 * t284;
t303 = t310 * t194 + t312 * t195 - t309 * t285;
t302 = t310 * t196 + t312 * t197 - t309 * t284;
t301 = t312 * t194 + t313 * t195 - t311 * t285;
t300 = t312 * t196 + t313 * t197 - t311 * t284;
t299 = t308 * t231 + (t309 * t241 + t311 * t243) * t230;
t298 = t309 * t231 + (t310 * t241 + t312 * t243) * t230;
t297 = t311 * t231 + (t312 * t241 + t313 * t243) * t230;
t238 = sin(pkin(9));
t292 = pkin(2) * t238;
t239 = cos(pkin(9));
t291 = pkin(2) * t239;
t290 = Icges(2,4) * t242;
t289 = Icges(3,4) * t238;
t288 = Icges(3,4) * t239;
t287 = Icges(4,4) * t230;
t286 = Icges(4,4) * t231;
t278 = rSges(7,2) * t194 + t307 * t195 - t306 * t285;
t277 = t196 * rSges(7,2) + t307 * t197 - t306 * t284;
t276 = t306 * t231 + (rSges(7,2) * t241 + t307 * t243) * t230;
t167 = -pkin(7) * t244 + t242 * t291;
t221 = t242 * pkin(1) - qJ(2) * t244;
t275 = -t167 - t221;
t274 = qJD(4) * t230;
t273 = qJD(6) * t230;
t272 = V_base(4) * t221 + V_base(3);
t271 = V_base(5) * pkin(6) + V_base(1);
t226 = qJD(3) * t242 + V_base(4);
t232 = V_base(6) + qJD(1);
t268 = qJD(2) * t242 + t271;
t267 = V_base(5) * t292 + t268;
t266 = pkin(3) * t231 + pkin(8) * t230;
t225 = -qJD(3) * t244 + V_base(5);
t265 = rSges(3,1) * t239 - rSges(3,2) * t238;
t264 = rSges(4,1) * t231 - rSges(4,2) * t230;
t263 = Icges(3,1) * t239 - t289;
t262 = Icges(4,1) * t231 - t287;
t261 = -Icges(3,2) * t238 + t288;
t260 = -Icges(4,2) * t230 + t286;
t259 = Icges(3,5) * t239 - Icges(3,6) * t238;
t258 = Icges(4,5) * t231 - Icges(4,6) * t230;
t223 = pkin(1) * t244 + t242 * qJ(2);
t257 = -qJD(2) * t244 + t232 * t223 + V_base(2);
t256 = (-Icges(4,3) * t244 + t242 * t258) * t225 + (Icges(4,3) * t242 + t244 * t258) * t226 + (Icges(4,5) * t230 + Icges(4,6) * t231) * t232;
t168 = pkin(7) * t242 + t244 * t291;
t255 = V_base(4) * t167 + (-t168 - t223) * V_base(5) + t272;
t254 = (-Icges(3,3) * t244 + t242 * t259) * V_base(5) + (Icges(3,3) * t242 + t244 * t259) * V_base(4) + (Icges(3,5) * t238 + Icges(3,6) * t239) * t232;
t190 = t266 * t242;
t204 = pkin(3) * t230 - pkin(8) * t231;
t253 = t225 * t204 + (-t190 + t275) * t232 + t267;
t191 = t266 * t244;
t252 = t226 * t190 - t191 * t225 + t255;
t251 = t232 * t168 + (-pkin(6) - t292) * V_base(4) + t257;
t189 = (pkin(4) * t243 + qJ(5) * t241) * t230;
t192 = t242 * t274 + t225;
t250 = qJD(5) * t196 + t192 * t189 + t253;
t144 = pkin(4) * t195 + qJ(5) * t194;
t193 = t244 * t274 + t226;
t249 = qJD(5) * t230 * t241 + t193 * t144 + t252;
t248 = t232 * t191 - t226 * t204 + t251;
t145 = pkin(4) * t197 + qJ(5) * t196;
t206 = -qJD(4) * t231 + t232;
t247 = qJD(5) * t194 + t206 * t145 + t248;
t172 = -Icges(4,6) * t244 + t242 * t260;
t173 = Icges(4,6) * t242 + t244 * t260;
t174 = -Icges(4,5) * t244 + t242 * t262;
t175 = Icges(4,5) * t242 + t244 * t262;
t201 = Icges(4,2) * t231 + t287;
t202 = Icges(4,1) * t230 + t286;
t246 = (-t173 * t230 + t175 * t231) * t226 + (-t172 * t230 + t174 * t231) * t225 + (-t201 * t230 + t202 * t231) * t232;
t181 = -Icges(3,6) * t244 + t242 * t261;
t182 = Icges(3,6) * t242 + t244 * t261;
t183 = -Icges(3,5) * t244 + t242 * t263;
t184 = Icges(3,5) * t242 + t244 * t263;
t211 = Icges(3,2) * t239 + t289;
t212 = Icges(3,1) * t238 + t288;
t245 = (-t182 * t238 + t184 * t239) * V_base(4) + (-t181 * t238 + t183 * t239) * V_base(5) + (-t211 * t238 + t212 * t239) * t232;
t235 = Icges(2,4) * t244;
t224 = rSges(2,1) * t244 - t242 * rSges(2,2);
t222 = t242 * rSges(2,1) + rSges(2,2) * t244;
t219 = Icges(2,1) * t244 - t290;
t218 = Icges(2,1) * t242 + t235;
t217 = -Icges(2,2) * t242 + t235;
t216 = Icges(2,2) * t244 + t290;
t215 = Icges(2,5) * t244 - Icges(2,6) * t242;
t214 = Icges(2,5) * t242 + Icges(2,6) * t244;
t213 = rSges(3,1) * t238 + rSges(3,2) * t239;
t209 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t208 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t207 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = rSges(4,1) * t230 + rSges(4,2) * t231;
t186 = t242 * rSges(3,3) + t244 * t265;
t185 = -rSges(3,3) * t244 + t242 * t265;
t177 = t242 * rSges(4,3) + t244 * t264;
t176 = -rSges(4,3) * t244 + t242 * t264;
t166 = -rSges(5,3) * t231 + (rSges(5,1) * t243 - rSges(5,2) * t241) * t230;
t165 = -rSges(6,2) * t231 + (rSges(6,1) * t243 + rSges(6,3) * t241) * t230;
t163 = V_base(5) * rSges(2,3) - t222 * t232 + t271;
t162 = t224 * t232 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t151 = t222 * V_base(4) - t224 * V_base(5) + V_base(3);
t143 = t197 * rSges(5,1) - t196 * rSges(5,2) + rSges(5,3) * t284;
t142 = t197 * rSges(6,1) + rSges(6,2) * t284 + t196 * rSges(6,3);
t140 = rSges(5,1) * t195 - rSges(5,2) * t194 + rSges(5,3) * t285;
t139 = rSges(6,1) * t195 + rSges(6,2) * t285 + rSges(6,3) * t194;
t117 = t213 * V_base(5) + (-t185 - t221) * t232 + t268;
t116 = t232 * t186 + (-pkin(6) - t213) * V_base(4) + t257;
t115 = t185 * V_base(4) + (-t186 - t223) * V_base(5) + t272;
t114 = t203 * t225 + (-t176 + t275) * t232 + t267;
t113 = t232 * t177 - t226 * t203 + t251;
t112 = t176 * t226 - t177 * t225 + t255;
t111 = -t140 * t206 + t166 * t192 + t253;
t110 = t206 * t143 - t193 * t166 + t248;
t109 = t140 * t193 - t143 * t192 + t252;
t108 = t165 * t192 + (-t139 - t144) * t206 + t250;
t107 = t206 * t142 + (-t165 - t189) * t193 + t247;
t106 = t139 * t193 + (-t142 - t145) * t192 + t249;
t105 = -t244 * t273 + t276 * t192 + (-t144 - t278) * t206 + t250;
t104 = -t242 * t273 + t277 * t206 + (-t189 - t276) * t193 + t247;
t103 = qJD(6) * t231 + t278 * t193 + (-t145 - t277) * t192 + t249;
t1 = t226 * (t242 * t256 + t244 * t246) / 0.2e1 + t225 * (t242 * t246 - t256 * t244) / 0.2e1 + m(1) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(2) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + ((t298 * t194 + t297 * t195 - t299 * t285) * t206 + (t302 * t194 + t300 * t195 - t304 * t285) * t193 + (t303 * t194 + t301 * t195 - t305 * t285) * t192) * t192 / 0.2e1 + ((t298 * t196 + t297 * t197 - t299 * t284) * t206 + (t302 * t196 + t300 * t197 - t304 * t284) * t193 + (t303 * t196 + t301 * t197 - t305 * t284) * t192) * t193 / 0.2e1 + ((t305 * t192 + t304 * t193 + t299 * t206) * t231 + ((t298 * t241 + t297 * t243) * t206 + (t302 * t241 + t300 * t243) * t193 + (t303 * t241 + t301 * t243) * t192) * t230) * t206 / 0.2e1 + ((t173 * t231 + t175 * t230) * t226 + (t172 * t231 + t174 * t230) * t225 + (t181 * t239 + t183 * t238 + t214) * V_base(5) + (t182 * t239 + t184 * t238 + t215) * V_base(4) + (t201 * t231 + t202 * t230 + t211 * t239 + t212 * t238 + Icges(2,3)) * t232) * t232 / 0.2e1 + (t215 * t232 + t242 * t254 + t244 * t245 + (-t242 * t216 + t218 * t244 + Icges(1,4)) * V_base(5) + (-t242 * t217 + t219 * t244 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t214 * t232 + t242 * t245 - t244 * t254 + (t216 * t244 + t242 * t218 + Icges(1,2)) * V_base(5) + (t217 * t244 + t242 * t219 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
