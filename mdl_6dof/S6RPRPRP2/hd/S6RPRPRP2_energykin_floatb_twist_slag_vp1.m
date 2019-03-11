% Calculate kinetic energy for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:26
% EndTime: 2019-03-09 03:04:28
% DurationCPUTime: 2.40s
% Computational Cost: add. (2076->290), mult. (1724->388), div. (0->0), fcn. (1570->10), ass. (0->149)
t315 = Icges(6,1) + Icges(7,1);
t314 = -Icges(6,4) + Icges(7,5);
t313 = Icges(7,4) + Icges(6,5);
t312 = Icges(6,2) + Icges(7,3);
t311 = -Icges(7,6) + Icges(6,6);
t310 = Icges(4,3) + Icges(5,3);
t309 = -Icges(6,3) - Icges(7,2);
t227 = qJ(3) + pkin(10);
t219 = sin(t227);
t221 = cos(t227);
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t308 = Icges(4,5) * t234 + Icges(5,5) * t221 - Icges(4,6) * t231 - Icges(5,6) * t219;
t307 = rSges(7,1) + pkin(5);
t306 = rSges(7,3) + qJ(6);
t228 = qJ(1) + pkin(9);
t222 = cos(t228);
t233 = cos(qJ(5));
t269 = t222 * t233;
t220 = sin(t228);
t230 = sin(qJ(5));
t272 = t220 * t230;
t171 = t221 * t272 + t269;
t270 = t222 * t230;
t271 = t220 * t233;
t172 = t221 * t271 - t270;
t274 = t219 * t220;
t305 = t312 * t171 + t314 * t172 - t311 * t274;
t173 = t221 * t270 - t271;
t174 = t221 * t269 + t272;
t273 = t219 * t222;
t304 = t312 * t173 + t314 * t174 - t311 * t273;
t303 = -t311 * t171 + t313 * t172 - t309 * t274;
t302 = -t311 * t173 + t313 * t174 - t309 * t273;
t301 = t314 * t171 + t315 * t172 + t313 * t274;
t300 = t314 * t173 + t315 * t174 + t313 * t273;
t299 = t311 * t221 + (t312 * t230 + t314 * t233) * t219;
t298 = t309 * t221 + (-t311 * t230 + t313 * t233) * t219;
t297 = -t313 * t221 + (t314 * t230 + t315 * t233) * t219;
t275 = Icges(5,4) * t221;
t250 = -Icges(5,2) * t219 + t275;
t140 = -Icges(5,6) * t222 + t220 * t250;
t141 = Icges(5,6) * t220 + t222 * t250;
t276 = Icges(5,4) * t219;
t252 = Icges(5,1) * t221 - t276;
t142 = -Icges(5,5) * t222 + t220 * t252;
t143 = Icges(5,5) * t220 + t222 * t252;
t277 = Icges(4,4) * t234;
t251 = -Icges(4,2) * t231 + t277;
t154 = -Icges(4,6) * t222 + t220 * t251;
t155 = Icges(4,6) * t220 + t222 * t251;
t278 = Icges(4,4) * t231;
t253 = Icges(4,1) * t234 - t278;
t158 = -Icges(4,5) * t222 + t220 * t253;
t159 = Icges(4,5) * t220 + t222 * t253;
t182 = Icges(5,2) * t221 + t276;
t185 = Icges(5,1) * t219 + t275;
t198 = -qJD(3) * t222 + V_base(5);
t199 = qJD(3) * t220 + V_base(4);
t203 = Icges(4,2) * t234 + t278;
t206 = Icges(4,1) * t231 + t277;
t223 = V_base(6) + qJD(1);
t294 = (-t182 * t219 + t185 * t221 - t203 * t231 + t206 * t234) * t223 + (-t141 * t219 + t143 * t221 - t155 * t231 + t159 * t234) * t199 + (-t140 * t219 + t142 * t221 - t154 * t231 + t158 * t234) * t198;
t293 = (Icges(4,5) * t231 + Icges(5,5) * t219 + Icges(4,6) * t234 + Icges(5,6) * t221) * t223 + (t310 * t220 + t308 * t222) * t199 + (t308 * t220 - t310 * t222) * t198;
t232 = sin(qJ(1));
t286 = pkin(1) * t232;
t235 = cos(qJ(1));
t285 = pkin(1) * t235;
t284 = pkin(3) * t231;
t283 = pkin(3) * t234;
t282 = -pkin(6) - qJ(2);
t280 = Icges(2,4) * t232;
t279 = Icges(3,4) * t220;
t268 = rSges(7,2) * t274 + t306 * t171 + t307 * t172;
t267 = rSges(7,2) * t273 + t306 * t173 + t307 * t174;
t266 = -rSges(7,2) * t221 + (t306 * t230 + t307 * t233) * t219;
t265 = qJD(5) * t219;
t264 = t223 * t285 + V_base(2);
t263 = V_base(5) * pkin(6) + V_base(1);
t192 = pkin(2) * t220 - pkin(7) * t222;
t260 = -t192 - t286;
t259 = V_base(5) * qJ(2) + t263;
t258 = V_base(4) * t286 + qJD(2) + V_base(3);
t136 = -qJ(4) * t222 + t220 * t283;
t257 = -t136 + t260;
t256 = pkin(4) * t221 + pkin(8) * t219;
t255 = rSges(4,1) * t234 - rSges(4,2) * t231;
t254 = rSges(5,1) * t221 - rSges(5,2) * t219;
t247 = qJD(4) * t220 + t198 * t284 + t259;
t193 = pkin(2) * t222 + pkin(7) * t220;
t244 = t223 * t193 + t282 * V_base(4) + t264;
t243 = V_base(4) * t192 + (-t193 - t285) * V_base(5) + t258;
t242 = t199 * t136 + t243;
t137 = qJ(4) * t220 + t222 * t283;
t241 = -qJD(4) * t222 + t223 * t137 + t244;
t167 = t256 * t220;
t191 = pkin(4) * t219 - pkin(8) * t221;
t240 = t198 * t191 + (-t167 + t257) * t223 + t247;
t168 = t256 * t222;
t239 = t199 * t167 + (-t137 - t168) * t198 + t242;
t238 = t223 * t168 + (-t191 - t284) * t199 + t241;
t225 = Icges(2,4) * t235;
t217 = Icges(3,4) * t222;
t211 = rSges(2,1) * t235 - t232 * rSges(2,2);
t210 = t232 * rSges(2,1) + rSges(2,2) * t235;
t209 = rSges(4,1) * t231 + rSges(4,2) * t234;
t208 = Icges(2,1) * t235 - t280;
t207 = Icges(2,1) * t232 + t225;
t205 = -Icges(2,2) * t232 + t225;
t204 = Icges(2,2) * t235 + t280;
t197 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t196 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t195 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t194 = -qJD(5) * t221 + t223;
t190 = rSges(3,1) * t222 - rSges(3,2) * t220;
t189 = rSges(3,1) * t220 + rSges(3,2) * t222;
t188 = rSges(5,1) * t219 + rSges(5,2) * t221;
t187 = Icges(3,1) * t222 - t279;
t186 = Icges(3,1) * t220 + t217;
t184 = -Icges(3,2) * t220 + t217;
t183 = Icges(3,2) * t222 + t279;
t170 = t222 * t265 + t199;
t169 = t220 * t265 + t198;
t165 = rSges(4,3) * t220 + t222 * t255;
t164 = -rSges(4,3) * t222 + t220 * t255;
t163 = -rSges(6,3) * t221 + (rSges(6,1) * t233 - rSges(6,2) * t230) * t219;
t161 = V_base(5) * rSges(2,3) - t210 * t223 + t263;
t160 = t211 * t223 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t147 = t210 * V_base(4) - t211 * V_base(5) + V_base(3);
t145 = rSges(5,3) * t220 + t222 * t254;
t144 = -rSges(5,3) * t222 + t220 * t254;
t133 = V_base(5) * rSges(3,3) + (-t189 - t286) * t223 + t259;
t132 = t190 * t223 + (-rSges(3,3) + t282) * V_base(4) + t264;
t128 = V_base(4) * t189 + (-t190 - t285) * V_base(5) + t258;
t127 = rSges(6,1) * t174 - rSges(6,2) * t173 + rSges(6,3) * t273;
t125 = rSges(6,1) * t172 - rSges(6,2) * t171 + rSges(6,3) * t274;
t111 = t198 * t209 + (-t164 + t260) * t223 + t259;
t110 = t165 * t223 - t199 * t209 + t244;
t109 = t199 * t164 - t198 * t165 + t243;
t108 = t188 * t198 + (-t144 + t257) * t223 + t247;
t107 = t145 * t223 + (-t188 - t284) * t199 + t241;
t106 = t199 * t144 + (-t137 - t145) * t198 + t242;
t105 = -t125 * t194 + t163 * t169 + t240;
t104 = t127 * t194 - t163 * t170 + t238;
t103 = t170 * t125 - t169 * t127 + t239;
t102 = qJD(6) * t173 + t169 * t266 - t194 * t268 + t240;
t101 = qJD(6) * t171 - t170 * t266 + t194 * t267 + t238;
t100 = qJD(6) * t219 * t230 - t169 * t267 + t170 * t268 + t239;
t1 = m(1) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(4) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + ((t171 * t299 + t172 * t297 + t274 * t298) * t194 + (t171 * t304 + t172 * t300 + t274 * t302) * t170 + (t305 * t171 + t301 * t172 + t303 * t274) * t169) * t169 / 0.2e1 + ((t173 * t299 + t174 * t297 + t273 * t298) * t194 + (t304 * t173 + t300 * t174 + t302 * t273) * t170 + (t173 * t305 + t301 * t174 + t303 * t273) * t169) * t170 / 0.2e1 + ((-t169 * t303 - t170 * t302 - t194 * t298) * t221 + ((t230 * t299 + t233 * t297) * t194 + (t230 * t304 + t233 * t300) * t170 + (t230 * t305 + t301 * t233) * t169) * t219) * t194 / 0.2e1 + (t294 * t220 - t293 * t222) * t198 / 0.2e1 + (t293 * t220 + t294 * t222) * t199 / 0.2e1 + ((-t183 * t220 + t186 * t222 - t232 * t204 + t207 * t235 + Icges(1,4)) * V_base(5) + (-t184 * t220 + t187 * t222 - t232 * t205 + t208 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t183 * t222 + t186 * t220 + t204 * t235 + t232 * t207 + Icges(1,2)) * V_base(5) + (t184 * t222 + t187 * t220 + t205 * t235 + t232 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t141 * t221 + t143 * t219 + t155 * t234 + t231 * t159) * t199 + (t140 * t221 + t142 * t219 + t154 * t234 + t158 * t231) * t198 + (t182 * t221 + t185 * t219 + t203 * t234 + t206 * t231 + Icges(2,3) + Icges(3,3)) * t223) * t223 / 0.2e1 + t223 * V_base(5) * (Icges(2,5) * t232 + Icges(3,5) * t220 + Icges(2,6) * t235 + Icges(3,6) * t222) + t223 * V_base(4) * (Icges(2,5) * t235 + Icges(3,5) * t222 - Icges(2,6) * t232 - Icges(3,6) * t220) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
