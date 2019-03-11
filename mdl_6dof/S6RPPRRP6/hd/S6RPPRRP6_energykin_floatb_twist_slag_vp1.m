% Calculate kinetic energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:08
% EndTime: 2019-03-09 02:10:10
% DurationCPUTime: 2.36s
% Computational Cost: add. (847->251), mult. (1544->333), div. (0->0), fcn. (1392->6), ass. (0->129)
t292 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t291 = Icges(6,1) + Icges(7,1);
t290 = -Icges(6,4) + Icges(7,5);
t289 = Icges(7,4) + Icges(6,5);
t288 = Icges(6,2) + Icges(7,3);
t287 = Icges(7,6) - Icges(6,6);
t286 = Icges(6,3) + Icges(7,2);
t285 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t284 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t283 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t282 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t281 = rSges(7,1) + pkin(5);
t280 = rSges(7,3) + qJ(6);
t209 = sin(qJ(1));
t279 = t292 * t209;
t212 = cos(qJ(1));
t278 = t292 * t212;
t208 = sin(qJ(4));
t210 = cos(qJ(5));
t242 = t212 * t210;
t207 = sin(qJ(5));
t246 = t209 * t207;
t146 = t208 * t246 - t242;
t245 = t209 * t210;
t247 = t207 * t212;
t147 = t208 * t245 + t247;
t211 = cos(qJ(4));
t244 = t209 * t211;
t277 = t146 * t288 + t147 * t290 - t244 * t287;
t148 = t208 * t247 + t245;
t149 = t208 * t242 - t246;
t243 = t211 * t212;
t276 = t148 * t288 + t149 * t290 - t243 * t287;
t275 = t146 * t287 + t147 * t289 - t244 * t286;
t274 = t148 * t287 + t149 * t289 - t243 * t286;
t273 = t290 * t146 + t147 * t291 - t289 * t244;
t272 = t290 * t148 + t149 * t291 - t289 * t243;
t271 = (t207 * t288 + t210 * t290) * t211 + t287 * t208;
t270 = (t207 * t287 + t210 * t289) * t211 + t286 * t208;
t269 = (t290 * t207 + t210 * t291) * t211 + t289 * t208;
t268 = t212 * t285 - t279;
t267 = t209 * t285 + t278;
t266 = -t212 * t282 - t279;
t265 = t209 * t282 - t278;
t262 = -pkin(2) - pkin(6);
t257 = pkin(7) * t209;
t256 = pkin(7) * t212;
t254 = Icges(5,4) * t208;
t253 = Icges(5,4) * t211;
t249 = qJ(3) * t209;
t248 = qJ(3) * t212;
t241 = -rSges(7,2) * t244 + t280 * t146 + t281 * t147;
t240 = -rSges(7,2) * t243 + t280 * t148 + t281 * t149;
t239 = rSges(7,2) * t208 + (t280 * t207 + t281 * t210) * t211;
t238 = qJD(5) * t211;
t180 = t209 * pkin(1) - qJ(2) * t212;
t237 = V_base(4) * t180 + V_base(3);
t236 = V_base(5) * pkin(6) + V_base(1);
t191 = qJD(4) * t212 + V_base(5);
t197 = V_base(6) + qJD(1);
t233 = V_base(4) * t249 + t237;
t232 = qJD(2) * t209 + t236;
t231 = -t180 - t249;
t185 = pkin(1) * t212 + t209 * qJ(2);
t230 = -t185 - t248;
t229 = pkin(4) * t208 - pkin(8) * t211;
t192 = -qJD(4) * t209 + V_base(4);
t228 = rSges(5,1) * t208 + rSges(5,2) * t211;
t227 = Icges(5,1) * t208 + t253;
t226 = Icges(5,2) * t211 + t254;
t225 = Icges(5,5) * t208 + Icges(5,6) * t211;
t224 = -qJD(2) * t212 + t197 * t185 + V_base(2);
t223 = V_base(5) * pkin(2) + qJD(3) * t212 + t232;
t222 = t231 - t256;
t221 = V_base(5) * pkin(3) + t223;
t220 = qJD(3) * t209 + t197 * t248 + t224;
t219 = (Icges(5,3) * t212 + t209 * t225) * t191 + (-Icges(5,3) * t209 + t212 * t225) * t192 + (Icges(5,5) * t211 - Icges(5,6) * t208) * t197;
t218 = V_base(4) * t256 + t233 + (t230 + t257) * V_base(5);
t217 = (-pkin(3) + t262) * V_base(4) + t220;
t152 = t229 * t209;
t153 = t229 * t212;
t216 = t192 * t152 - t191 * t153 + t218;
t190 = pkin(4) * t211 + pkin(8) * t208;
t215 = t191 * t190 + (-t152 + t222) * t197 + t221;
t214 = -t192 * t190 + t217 + (t153 - t257) * t197;
t133 = Icges(5,6) * t212 + t209 * t226;
t134 = -Icges(5,6) * t209 + t212 * t226;
t137 = Icges(5,5) * t212 + t209 * t227;
t138 = -Icges(5,5) * t209 + t212 * t227;
t169 = -Icges(5,2) * t208 + t253;
t176 = Icges(5,1) * t211 - t254;
t213 = (t134 * t211 + t138 * t208) * t192 + (t133 * t211 + t137 * t208) * t191 + (t169 * t211 + t176 * t208) * t197;
t188 = rSges(2,1) * t212 - t209 * rSges(2,2);
t187 = -rSges(3,2) * t212 + t209 * rSges(3,3);
t186 = -rSges(4,2) * t212 + t209 * rSges(4,3);
t184 = rSges(5,1) * t211 - rSges(5,2) * t208;
t183 = t209 * rSges(2,1) + rSges(2,2) * t212;
t182 = -t209 * rSges(3,2) - rSges(3,3) * t212;
t181 = t209 * rSges(4,2) + rSges(4,3) * t212;
t179 = qJD(5) * t208 + t197;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t145 = -t212 * t238 + t192;
t144 = -t209 * t238 + t191;
t142 = -t209 * rSges(5,3) + t212 * t228;
t141 = rSges(6,3) * t208 + (rSges(6,1) * t210 - rSges(6,2) * t207) * t211;
t139 = rSges(5,3) * t212 + t209 * t228;
t124 = V_base(5) * rSges(2,3) - t183 * t197 + t236;
t123 = t188 * t197 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t122 = t183 * V_base(4) - t188 * V_base(5) + V_base(3);
t119 = t149 * rSges(6,1) - t148 * rSges(6,2) - rSges(6,3) * t243;
t117 = rSges(6,1) * t147 - rSges(6,2) * t146 - rSges(6,3) * t244;
t103 = V_base(5) * rSges(3,1) + (-t180 - t182) * t197 + t232;
t102 = t197 * t187 + (-rSges(3,1) - pkin(6)) * V_base(4) + t224;
t101 = t182 * V_base(4) + (-t185 - t187) * V_base(5) + t237;
t100 = V_base(5) * rSges(4,1) + (-t186 + t231) * t197 + t223;
t99 = t197 * t181 + (-rSges(4,1) + t262) * V_base(4) + t220;
t98 = V_base(4) * t186 + (-t181 + t230) * V_base(5) + t233;
t97 = t191 * t184 + (-t139 + t222) * t197 + t221;
t96 = -t192 * t184 + (t142 - t257) * t197 + t217;
t95 = t192 * t139 - t191 * t142 + t218;
t94 = -t179 * t117 + t144 * t141 + t215;
t93 = t179 * t119 - t145 * t141 + t214;
t92 = t145 * t117 - t144 * t119 + t216;
t91 = qJD(6) * t148 + t144 * t239 - t179 * t241 + t215;
t90 = qJD(6) * t146 - t145 * t239 + t179 * t240 + t214;
t89 = qJD(6) * t211 * t207 - t144 * t240 + t145 * t241 + t216;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(7) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(6) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t191 * (t213 * t209 + t219 * t212) / 0.2e1 + t192 * (-t219 * t209 + t213 * t212) / 0.2e1 + ((t271 * t146 + t269 * t147 - t270 * t244) * t179 + (t276 * t146 + t272 * t147 - t274 * t244) * t145 + (t277 * t146 + t273 * t147 - t275 * t244) * t144) * t144 / 0.2e1 + ((t271 * t148 + t269 * t149 - t270 * t243) * t179 + (t276 * t148 + t272 * t149 - t274 * t243) * t145 + (t277 * t148 + t273 * t149 - t275 * t243) * t144) * t145 / 0.2e1 + (((t271 * t207 + t269 * t210) * t179 + (t276 * t207 + t272 * t210) * t145 + (t277 * t207 + t273 * t210) * t144) * t211 + (t275 * t144 + t274 * t145 + t270 * t179) * t208) * t179 / 0.2e1 + ((-t134 * t208 + t138 * t211) * t192 + (-t133 * t208 + t137 * t211) * t191 + (-t208 * t169 + t211 * t176 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t197) * t197 / 0.2e1 + ((t266 * t209 + t267 * t212 + Icges(1,4)) * V_base(5) + (t265 * t209 + t268 * t212 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t209 - t266 * t212 + Icges(1,2)) * V_base(5) + (t268 * t209 - t265 * t212 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t197 * (t209 * t284 - t212 * t283) + V_base(4) * t197 * (t209 * t283 + t212 * t284) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
