% Calculate kinetic energy for
% S6RPRPRP1
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:15
% EndTime: 2019-03-09 03:01:18
% DurationCPUTime: 2.54s
% Computational Cost: add. (2108->295), mult. (1733->395), div. (0->0), fcn. (1569->10), ass. (0->150)
t316 = Icges(6,1) + Icges(7,1);
t315 = Icges(6,4) + Icges(7,4);
t314 = -Icges(7,5) - Icges(6,5);
t313 = Icges(6,2) + Icges(7,2);
t312 = -Icges(7,6) - Icges(6,6);
t311 = Icges(4,3) + Icges(5,3);
t310 = -Icges(7,3) - Icges(6,3);
t225 = qJ(3) + pkin(10);
t217 = sin(t225);
t219 = cos(t225);
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t309 = Icges(4,5) * t233 + Icges(5,5) * t219 - Icges(4,6) * t230 - Icges(5,6) * t217;
t226 = qJ(1) + pkin(9);
t220 = cos(t226);
t232 = cos(qJ(5));
t270 = t220 * t232;
t218 = sin(t226);
t229 = sin(qJ(5));
t273 = t218 * t229;
t170 = -t219 * t273 - t270;
t271 = t220 * t229;
t272 = t218 * t232;
t171 = t219 * t272 - t271;
t275 = t217 * t218;
t308 = -t312 * t170 - t314 * t171 - t310 * t275;
t172 = -t219 * t271 + t272;
t173 = t219 * t270 + t273;
t274 = t217 * t220;
t307 = -t312 * t172 - t314 * t173 - t310 * t274;
t306 = t313 * t170 + t315 * t171 - t312 * t275;
t305 = t313 * t172 + t315 * t173 - t312 * t274;
t304 = t315 * t170 + t316 * t171 - t314 * t275;
t303 = t315 * t172 + t316 * t173 - t314 * t274;
t302 = t310 * t219 + (t312 * t229 - t314 * t232) * t217;
t301 = t312 * t219 + (-t313 * t229 + t315 * t232) * t217;
t300 = t314 * t219 + (-t315 * t229 + t316 * t232) * t217;
t276 = Icges(5,4) * t219;
t250 = -Icges(5,2) * t217 + t276;
t139 = -Icges(5,6) * t220 + t218 * t250;
t140 = Icges(5,6) * t218 + t220 * t250;
t277 = Icges(5,4) * t217;
t252 = Icges(5,1) * t219 - t277;
t141 = -Icges(5,5) * t220 + t218 * t252;
t142 = Icges(5,5) * t218 + t220 * t252;
t278 = Icges(4,4) * t233;
t251 = -Icges(4,2) * t230 + t278;
t153 = -Icges(4,6) * t220 + t218 * t251;
t154 = Icges(4,6) * t218 + t220 * t251;
t279 = Icges(4,4) * t230;
t253 = Icges(4,1) * t233 - t279;
t157 = -Icges(4,5) * t220 + t218 * t253;
t158 = Icges(4,5) * t218 + t220 * t253;
t180 = Icges(5,2) * t219 + t277;
t183 = Icges(5,1) * t217 + t276;
t196 = -qJD(3) * t220 + V_base(5);
t197 = qJD(3) * t218 + V_base(4);
t201 = Icges(4,2) * t233 + t279;
t204 = Icges(4,1) * t230 + t278;
t221 = V_base(6) + qJD(1);
t297 = (-t180 * t217 + t183 * t219 - t201 * t230 + t204 * t233) * t221 + (-t140 * t217 + t142 * t219 - t154 * t230 + t158 * t233) * t197 + (-t139 * t217 + t141 * t219 - t153 * t230 + t157 * t233) * t196;
t296 = (Icges(4,5) * t230 + Icges(5,5) * t217 + Icges(4,6) * t233 + Icges(5,6) * t219) * t221 + (t311 * t218 + t309 * t220) * t197 + (t309 * t218 - t311 * t220) * t196;
t231 = sin(qJ(1));
t289 = pkin(1) * t231;
t234 = cos(qJ(1));
t288 = pkin(1) * t234;
t287 = pkin(3) * t230;
t286 = pkin(3) * t233;
t285 = pkin(5) * t232;
t284 = -pkin(6) - qJ(2);
t281 = Icges(2,4) * t231;
t280 = Icges(3,4) * t218;
t246 = qJ(6) * t217 + t219 * t285;
t269 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t275 - pkin(5) * t271 + t218 * t246;
t268 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t274 + pkin(5) * t273 + t220 * t246;
t267 = (-qJ(6) - rSges(7,3)) * t219 + (rSges(7,1) * t232 - rSges(7,2) * t229 + t285) * t217;
t266 = qJD(5) * t217;
t265 = qJD(6) * t217;
t264 = t221 * t288 + V_base(2);
t263 = V_base(5) * pkin(6) + V_base(1);
t190 = pkin(2) * t218 - pkin(7) * t220;
t260 = -t190 - t289;
t259 = V_base(5) * qJ(2) + t263;
t258 = V_base(4) * t289 + qJD(2) + V_base(3);
t135 = -qJ(4) * t220 + t218 * t286;
t257 = -t135 + t260;
t256 = pkin(4) * t219 + pkin(8) * t217;
t255 = rSges(4,1) * t233 - rSges(4,2) * t230;
t254 = rSges(5,1) * t219 - rSges(5,2) * t217;
t247 = qJD(4) * t218 + t196 * t287 + t259;
t191 = pkin(2) * t220 + pkin(7) * t218;
t243 = t221 * t191 + t284 * V_base(4) + t264;
t242 = V_base(4) * t190 + (-t191 - t288) * V_base(5) + t258;
t241 = t197 * t135 + t242;
t136 = qJ(4) * t218 + t220 * t286;
t240 = -qJD(4) * t220 + t221 * t136 + t243;
t166 = t256 * t218;
t189 = pkin(4) * t217 - pkin(8) * t219;
t239 = t196 * t189 + (-t166 + t257) * t221 + t247;
t167 = t256 * t220;
t238 = t197 * t166 + (-t136 - t167) * t196 + t241;
t237 = t221 * t167 + (-t189 - t287) * t197 + t240;
t223 = Icges(2,4) * t234;
t214 = Icges(3,4) * t220;
t209 = rSges(2,1) * t234 - t231 * rSges(2,2);
t208 = t231 * rSges(2,1) + rSges(2,2) * t234;
t207 = rSges(4,1) * t230 + rSges(4,2) * t233;
t206 = Icges(2,1) * t234 - t281;
t205 = Icges(2,1) * t231 + t223;
t203 = -Icges(2,2) * t231 + t223;
t202 = Icges(2,2) * t234 + t281;
t195 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t194 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t193 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t192 = -qJD(5) * t219 + t221;
t188 = rSges(3,1) * t220 - rSges(3,2) * t218;
t187 = rSges(3,1) * t218 + rSges(3,2) * t220;
t186 = rSges(5,1) * t217 + rSges(5,2) * t219;
t185 = Icges(3,1) * t220 - t280;
t184 = Icges(3,1) * t218 + t214;
t182 = -Icges(3,2) * t218 + t214;
t181 = Icges(3,2) * t220 + t280;
t169 = t220 * t266 + t197;
t168 = t218 * t266 + t196;
t164 = rSges(4,3) * t218 + t220 * t255;
t163 = -rSges(4,3) * t220 + t218 * t255;
t162 = -rSges(6,3) * t219 + (rSges(6,1) * t232 - rSges(6,2) * t229) * t217;
t160 = V_base(5) * rSges(2,3) - t208 * t221 + t263;
t159 = t209 * t221 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t146 = t208 * V_base(4) - t209 * V_base(5) + V_base(3);
t144 = rSges(5,3) * t218 + t220 * t254;
t143 = -rSges(5,3) * t220 + t218 * t254;
t131 = V_base(5) * rSges(3,3) + (-t187 - t289) * t221 + t259;
t130 = t188 * t221 + (-rSges(3,3) + t284) * V_base(4) + t264;
t128 = V_base(4) * t187 + (-t188 - t288) * V_base(5) + t258;
t127 = rSges(6,1) * t173 + rSges(6,2) * t172 + rSges(6,3) * t274;
t125 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t275;
t109 = t196 * t207 + (-t163 + t260) * t221 + t259;
t108 = t164 * t221 - t197 * t207 + t243;
t107 = t197 * t163 - t196 * t164 + t242;
t106 = t186 * t196 + (-t143 + t257) * t221 + t247;
t105 = t144 * t221 + (-t186 - t287) * t197 + t240;
t104 = t197 * t143 + (-t136 - t144) * t196 + t241;
t103 = -t125 * t192 + t162 * t168 + t239;
t102 = t127 * t192 - t162 * t169 + t237;
t101 = t169 * t125 - t168 * t127 + t238;
t100 = t168 * t267 - t192 * t269 + t220 * t265 + t239;
t99 = -t169 * t267 + t192 * t268 + t218 * t265 + t237;
t98 = -qJD(6) * t219 - t168 * t268 + t169 * t269 + t238;
t1 = m(2) * (t146 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(3) * (t128 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(1) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + ((t170 * t301 + t171 * t300 + t275 * t302) * t192 + (t170 * t305 + t171 * t303 + t275 * t307) * t169 + (t306 * t170 + t304 * t171 + t308 * t275) * t168) * t168 / 0.2e1 + ((t172 * t301 + t173 * t300 + t274 * t302) * t192 + (t305 * t172 + t303 * t173 + t307 * t274) * t169 + (t306 * t172 + t304 * t173 + t274 * t308) * t168) * t169 / 0.2e1 + ((-t168 * t308 - t307 * t169 - t302 * t192) * t219 + ((-t229 * t301 + t232 * t300) * t192 + (-t229 * t305 + t232 * t303) * t169 + (-t229 * t306 + t304 * t232) * t168) * t217) * t192 / 0.2e1 + (t297 * t218 - t296 * t220) * t196 / 0.2e1 + (t296 * t218 + t297 * t220) * t197 / 0.2e1 + ((-t181 * t218 + t184 * t220 - t231 * t202 + t205 * t234 + Icges(1,4)) * V_base(5) + (-t182 * t218 + t185 * t220 - t231 * t203 + t206 * t234 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t181 * t220 + t184 * t218 + t202 * t234 + t231 * t205 + Icges(1,2)) * V_base(5) + (t182 * t220 + t185 * t218 + t203 * t234 + t231 * t206 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t140 * t219 + t142 * t217 + t154 * t233 + t158 * t230) * t197 + (t139 * t219 + t141 * t217 + t153 * t233 + t157 * t230) * t196 + (t180 * t219 + t183 * t217 + t201 * t233 + t204 * t230 + Icges(2,3) + Icges(3,3)) * t221) * t221 / 0.2e1 + t221 * V_base(5) * (Icges(2,5) * t231 + Icges(3,5) * t218 + Icges(2,6) * t234 + Icges(3,6) * t220) + t221 * V_base(4) * (Icges(2,5) * t234 + Icges(3,5) * t220 - Icges(2,6) * t231 - Icges(3,6) * t218) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
