% Calculate kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:20
% EndTime: 2019-12-05 16:02:24
% DurationCPUTime: 3.38s
% Computational Cost: add. (1513->314), mult. (3496->451), div. (0->0), fcn. (4016->10), ass. (0->139)
t294 = Icges(3,1) + Icges(4,2);
t293 = Icges(3,4) + Icges(4,6);
t292 = Icges(3,5) - Icges(4,4);
t291 = Icges(3,2) + Icges(4,3);
t290 = Icges(3,6) - Icges(4,5);
t289 = Icges(3,3) + Icges(4,1);
t240 = sin(pkin(9));
t242 = cos(pkin(9));
t246 = sin(qJ(2));
t243 = cos(pkin(5));
t248 = cos(qJ(2));
t265 = t243 * t248;
t207 = t240 * t265 + t242 * t246;
t266 = t243 * t246;
t208 = -t240 * t266 + t242 * t248;
t241 = sin(pkin(5));
t271 = t240 * t241;
t286 = t291 * t207 - t293 * t208 - t290 * t271;
t205 = t240 * t246 - t242 * t265;
t206 = t240 * t248 + t242 * t266;
t270 = t241 * t242;
t285 = t291 * t205 - t293 * t206 + t290 * t270;
t284 = -t293 * t207 + t294 * t208 + t292 * t271;
t283 = -t293 * t205 + t294 * t206 - t292 * t270;
t282 = -t290 * t207 + t292 * t208 + t289 * t271;
t281 = -t290 * t205 + t292 * t206 - t289 * t270;
t280 = t289 * t243 + (t292 * t246 + t290 * t248) * t241;
t279 = t290 * t243 + (t293 * t246 + t291 * t248) * t241;
t278 = t292 * t243 + (t294 * t246 + t293 * t248) * t241;
t274 = cos(qJ(4));
t273 = pkin(6) * t243;
t272 = Icges(2,4) * t240;
t245 = sin(qJ(4));
t269 = t241 * t245;
t268 = t241 * t246;
t267 = t241 * t248;
t264 = qJD(2) * t241;
t263 = V_base(5) * qJ(1) + V_base(1);
t259 = qJD(1) + V_base(3);
t258 = t241 * t274;
t222 = t240 * t264 + V_base(4);
t233 = qJD(2) * t243 + V_base(6);
t179 = qJD(4) * t208 + t222;
t209 = qJD(4) * t268 + t233;
t221 = -t242 * t264 + V_base(5);
t178 = qJD(4) * t206 + t221;
t215 = pkin(1) * t240 - pkin(6) * t270;
t257 = -t215 * V_base(6) + V_base(5) * t273 + t263;
t216 = pkin(1) * t242 + pkin(6) * t271;
t256 = V_base(4) * t215 - V_base(5) * t216 + t259;
t255 = V_base(6) * t216 + V_base(2) + (-qJ(1) - t273) * V_base(4);
t214 = (pkin(2) * t246 - qJ(3) * t248) * t241;
t254 = qJD(3) * t207 + t221 * t214 + t257;
t173 = pkin(2) * t208 + qJ(3) * t207;
t253 = qJD(3) * t205 + t233 * t173 + t255;
t172 = pkin(2) * t206 + qJ(3) * t205;
t252 = -qJD(3) * t267 + t222 * t172 + t256;
t188 = -pkin(3) * t270 + pkin(7) * t206;
t217 = pkin(3) * t243 + pkin(7) * t268;
t251 = t221 * t217 + (-t172 - t188) * t233 + t254;
t187 = pkin(3) * t271 + pkin(7) * t208;
t250 = t233 * t187 + (-t214 - t217) * t222 + t253;
t249 = t222 * t188 + (-t173 - t187) * t221 + t252;
t247 = cos(qJ(5));
t244 = sin(qJ(5));
t238 = Icges(2,4) * t242;
t230 = rSges(2,1) * t242 - rSges(2,2) * t240;
t229 = rSges(2,1) * t240 + rSges(2,2) * t242;
t228 = Icges(2,1) * t242 - t272;
t227 = Icges(2,1) * t240 + t238;
t226 = -Icges(2,2) * t240 + t238;
t225 = Icges(2,2) * t242 + t272;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = t243 * t274 - t245 * t267;
t212 = t243 * t245 + t248 * t258;
t199 = t243 * rSges(4,1) + (-rSges(4,2) * t246 - rSges(4,3) * t248) * t241;
t198 = t243 * rSges(3,3) + (rSges(3,1) * t246 + rSges(3,2) * t248) * t241;
t190 = V_base(5) * rSges(2,3) - t229 * V_base(6) + t263;
t189 = t230 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t185 = t213 * t247 + t244 * t268;
t184 = -t213 * t244 + t247 * t268;
t183 = t205 * t245 - t242 * t258;
t182 = t205 * t274 + t242 * t269;
t181 = t207 * t245 + t240 * t258;
t180 = -t207 * t274 + t240 * t269;
t177 = t229 * V_base(4) - t230 * V_base(5) + t259;
t176 = qJD(5) * t212 + t209;
t174 = pkin(4) * t213 + pkin(8) * t212;
t170 = rSges(5,1) * t213 - rSges(5,2) * t212 + rSges(5,3) * t268;
t169 = Icges(5,1) * t213 - Icges(5,4) * t212 + Icges(5,5) * t268;
t168 = Icges(5,4) * t213 - Icges(5,2) * t212 + Icges(5,6) * t268;
t167 = Icges(5,5) * t213 - Icges(5,6) * t212 + Icges(5,3) * t268;
t166 = rSges(3,1) * t208 - rSges(3,2) * t207 + rSges(3,3) * t271;
t165 = rSges(3,1) * t206 - rSges(3,2) * t205 - rSges(3,3) * t270;
t164 = -rSges(4,1) * t270 - rSges(4,2) * t206 + rSges(4,3) * t205;
t163 = rSges(4,1) * t271 - rSges(4,2) * t208 + rSges(4,3) * t207;
t148 = t183 * t247 + t206 * t244;
t147 = -t183 * t244 + t206 * t247;
t146 = t181 * t247 + t208 * t244;
t145 = -t181 * t244 + t208 * t247;
t144 = qJD(5) * t180 + t179;
t143 = -qJD(5) * t182 + t178;
t142 = pkin(4) * t183 - pkin(8) * t182;
t141 = pkin(4) * t181 + pkin(8) * t180;
t140 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t212;
t139 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t212;
t138 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t212;
t137 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t212;
t136 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t206;
t135 = rSges(5,1) * t181 - rSges(5,2) * t180 + rSges(5,3) * t208;
t134 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t206;
t133 = Icges(5,1) * t181 - Icges(5,4) * t180 + Icges(5,5) * t208;
t132 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t206;
t131 = Icges(5,4) * t181 - Icges(5,2) * t180 + Icges(5,6) * t208;
t130 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t206;
t129 = Icges(5,5) * t181 - Icges(5,6) * t180 + Icges(5,3) * t208;
t128 = -t165 * t233 + t198 * t221 + t257;
t127 = t166 * t233 - t198 * t222 + t255;
t126 = rSges(6,1) * t148 + rSges(6,2) * t147 - rSges(6,3) * t182;
t125 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t180;
t124 = Icges(6,1) * t148 + Icges(6,4) * t147 - Icges(6,5) * t182;
t123 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t180;
t122 = Icges(6,4) * t148 + Icges(6,2) * t147 - Icges(6,6) * t182;
t121 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t180;
t120 = Icges(6,5) * t148 + Icges(6,6) * t147 - Icges(6,3) * t182;
t119 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t180;
t118 = t165 * t222 - t166 * t221 + t256;
t117 = t199 * t221 + (-t164 - t172) * t233 + t254;
t116 = t163 * t233 + (-t199 - t214) * t222 + t253;
t115 = t222 * t164 + (-t163 - t173) * t221 + t252;
t114 = -t136 * t209 + t170 * t178 + t251;
t113 = t135 * t209 - t170 * t179 + t250;
t112 = -t178 * t135 + t179 * t136 + t249;
t111 = -t126 * t176 + t140 * t143 - t142 * t209 + t174 * t178 + t251;
t110 = t125 * t176 - t140 * t144 + t141 * t209 - t174 * t179 + t250;
t109 = -t143 * t125 + t144 * t126 - t178 * t141 + t179 * t142 + t249;
t1 = m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(2) * (t177 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(3) * (t118 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t179 * ((t208 * t129 - t180 * t131 + t181 * t133) * t179 + (t130 * t208 - t132 * t180 + t134 * t181) * t178 + (t167 * t208 - t168 * t180 + t169 * t181) * t209) / 0.2e1 + t178 * ((t129 * t206 + t131 * t182 + t133 * t183) * t179 + (t206 * t130 + t182 * t132 + t183 * t134) * t178 + (t167 * t206 + t168 * t182 + t169 * t183) * t209) / 0.2e1 + t209 * ((t129 * t268 - t131 * t212 + t133 * t213) * t179 + (t130 * t268 - t132 * t212 + t134 * t213) * t178 + (t167 * t268 - t212 * t168 + t213 * t169) * t209) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t144 * ((t180 * t119 + t145 * t121 + t146 * t123) * t144 + (t120 * t180 + t122 * t145 + t124 * t146) * t143 + (t137 * t180 + t138 * t145 + t139 * t146) * t176) / 0.2e1 + t143 * ((-t119 * t182 + t121 * t147 + t123 * t148) * t144 + (-t182 * t120 + t147 * t122 + t148 * t124) * t143 + (-t137 * t182 + t138 * t147 + t139 * t148) * t176) / 0.2e1 + t176 * ((t119 * t212 + t121 * t184 + t123 * t185) * t144 + (t120 * t212 + t122 * t184 + t124 * t185) * t143 + (t212 * t137 + t184 * t138 + t185 * t139) * t176) / 0.2e1 + ((-t205 * t279 + t206 * t278 - t270 * t280) * t233 + (t205 * t286 + t206 * t284 - t270 * t282) * t222 + (t205 * t285 + t206 * t283 - t270 * t281) * t221) * t221 / 0.2e1 + ((-t207 * t279 + t208 * t278 + t271 * t280) * t233 + (t207 * t286 + t208 * t284 + t271 * t282) * t222 + (t207 * t285 + t208 * t283 + t271 * t281) * t221) * t222 / 0.2e1 + ((t221 * t281 + t222 * t282 + t233 * t280) * t243 + ((t246 * t278 + t248 * t279) * t233 + (t246 * t284 - t248 * t286) * t222 + (t246 * t283 - t248 * t285) * t221) * t241) * t233 / 0.2e1 + ((-t225 * t240 + t227 * t242 + Icges(1,4)) * V_base(5) + (-t226 * t240 + t228 * t242 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t225 * t242 + t227 * t240 + Icges(1,2)) * V_base(5) + (t226 * t242 + t228 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t240 + Icges(2,6) * t242 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t242 - Icges(2,6) * t240 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
