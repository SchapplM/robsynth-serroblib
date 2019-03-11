% Calculate kinetic energy for
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:28
% EndTime: 2019-03-09 01:41:32
% DurationCPUTime: 3.84s
% Computational Cost: add. (1763->303), mult. (1442->408), div. (0->0), fcn. (1222->10), ass. (0->159)
t310 = Icges(5,4) + Icges(6,6);
t309 = Icges(5,1) + Icges(6,2);
t308 = -Icges(5,2) - Icges(6,3);
t217 = pkin(10) + qJ(4);
t211 = cos(t217);
t307 = t310 * t211;
t209 = sin(t217);
t306 = t310 * t209;
t305 = Icges(6,4) - Icges(5,5);
t304 = Icges(6,5) - Icges(5,6);
t303 = t308 * t209 + t307;
t302 = t309 * t211 - t306;
t301 = Icges(6,1) + Icges(5,3);
t218 = qJ(1) + pkin(9);
t210 = sin(t218);
t212 = cos(t218);
t300 = t303 * t210 + t304 * t212;
t299 = -t304 * t210 + t303 * t212;
t298 = t302 * t210 + t305 * t212;
t297 = -t305 * t210 + t302 * t212;
t296 = t308 * t211 - t306;
t295 = t309 * t209 + t307;
t294 = t304 * t209 - t305 * t211;
t223 = sin(qJ(1));
t225 = cos(qJ(1));
t293 = Icges(2,5) * t223 + Icges(3,5) * t210 + Icges(2,6) * t225 + Icges(3,6) * t212;
t292 = Icges(2,5) * t225 + Icges(3,5) * t212 - Icges(2,6) * t223 - Icges(3,6) * t210;
t191 = -qJD(4) * t212 + V_base(5);
t192 = qJD(4) * t210 + V_base(4);
t213 = V_base(6) + qJD(1);
t291 = (t296 * t209 + t295 * t211) * t213 + (-t299 * t209 + t297 * t211) * t192 + (-t300 * t209 + t298 * t211) * t191;
t290 = (-t305 * t209 - t304 * t211) * t213 + (t301 * t210 + t294 * t212) * t192 + (t294 * t210 - t301 * t212) * t191;
t286 = pkin(1) * t223;
t285 = pkin(1) * t225;
t219 = sin(pkin(10));
t284 = pkin(3) * t219;
t283 = pkin(8) * t209;
t220 = cos(pkin(10));
t282 = pkin(3) * t220;
t281 = -pkin(6) - qJ(2);
t280 = Icges(2,4) * t223;
t279 = Icges(3,4) * t210;
t278 = Icges(4,4) * t219;
t277 = Icges(4,4) * t220;
t272 = t210 * t211;
t222 = sin(qJ(6));
t271 = t210 * t222;
t224 = cos(qJ(6));
t270 = t210 * t224;
t269 = t211 * t212;
t268 = t212 * t222;
t267 = t212 * t224;
t265 = qJD(5) * t209;
t264 = qJD(6) * t211;
t263 = t213 * t285 + V_base(2);
t262 = V_base(5) * pkin(6) + V_base(1);
t178 = pkin(2) * t210 - qJ(3) * t212;
t259 = -t178 - t286;
t180 = pkin(2) * t212 + qJ(3) * t210;
t258 = -t180 - t285;
t257 = V_base(5) * qJ(2) + t262;
t256 = V_base(4) * t286 + qJD(2) + V_base(3);
t116 = -pkin(7) * t212 + t210 * t282;
t255 = -t116 + t259;
t254 = qJD(3) * t210 + t257;
t253 = V_base(4) * t178 + t256;
t252 = rSges(4,1) * t220 - rSges(4,2) * t219;
t251 = rSges(5,1) * t211 - rSges(5,2) * t209;
t250 = -rSges(6,2) * t211 + rSges(6,3) * t209;
t249 = pkin(4) * t211 + qJ(5) * t209;
t248 = Icges(4,1) * t220 - t278;
t246 = -Icges(4,2) * t219 + t277;
t243 = Icges(4,5) * t220 - Icges(4,6) * t219;
t151 = t249 * t210;
t239 = -t151 + t255;
t238 = V_base(5) * t284 + t254;
t237 = -qJD(3) * t212 + t213 * t180 + t263;
t175 = pkin(4) * t209 - qJ(5) * t211;
t236 = t191 * t175 + t212 * t265 + t238;
t233 = (-Icges(4,3) * t212 + t210 * t243) * V_base(5) + (Icges(4,3) * t210 + t212 * t243) * V_base(4) + (Icges(4,5) * t219 + Icges(4,6) * t220) * t213;
t117 = pkin(7) * t210 + t212 * t282;
t232 = V_base(4) * t116 + (-t117 + t258) * V_base(5) + t253;
t231 = t213 * t117 + (t281 - t284) * V_base(4) + t237;
t152 = t249 * t212;
t230 = t213 * t152 + t210 * t265 + t231;
t229 = -qJD(5) * t211 + t192 * t151 + t232;
t138 = -Icges(4,6) * t212 + t210 * t246;
t139 = Icges(4,6) * t210 + t212 * t246;
t140 = -Icges(4,5) * t212 + t210 * t248;
t141 = Icges(4,5) * t210 + t212 * t248;
t189 = Icges(4,2) * t220 + t278;
t190 = Icges(4,1) * t219 + t277;
t226 = (-t139 * t219 + t141 * t220) * V_base(4) + (-t138 * t219 + t140 * t220) * V_base(5) + (-t189 * t219 + t190 * t220) * t213;
t215 = Icges(2,4) * t225;
t207 = Icges(3,4) * t212;
t201 = rSges(2,1) * t225 - t223 * rSges(2,2);
t200 = t223 * rSges(2,1) + rSges(2,2) * t225;
t199 = Icges(2,1) * t225 - t280;
t198 = Icges(2,1) * t223 + t215;
t197 = -Icges(2,2) * t223 + t215;
t196 = Icges(2,2) * t225 + t280;
t193 = rSges(4,1) * t219 + rSges(4,2) * t220;
t187 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t186 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t185 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t184 = qJD(6) * t209 + t213;
t181 = rSges(3,1) * t212 - rSges(3,2) * t210;
t179 = rSges(3,1) * t210 + rSges(3,2) * t212;
t177 = rSges(5,1) * t209 + rSges(5,2) * t211;
t176 = -rSges(6,2) * t209 - rSges(6,3) * t211;
t174 = Icges(3,1) * t212 - t279;
t173 = Icges(3,1) * t210 + t207;
t171 = -Icges(3,2) * t210 + t207;
t170 = Icges(3,2) * t212 + t279;
t161 = -pkin(5) * t212 + pkin(8) * t272;
t160 = pkin(5) * t210 + pkin(8) * t269;
t158 = t209 * t271 - t267;
t157 = t209 * t270 + t268;
t156 = t209 * t268 + t270;
t155 = t209 * t267 - t271;
t154 = t212 * t264 + t192;
t153 = t210 * t264 + t191;
t149 = rSges(7,3) * t209 + (-rSges(7,1) * t222 - rSges(7,2) * t224) * t211;
t148 = V_base(5) * rSges(2,3) - t200 * t213 + t262;
t147 = t201 * t213 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t146 = Icges(7,5) * t209 + (-Icges(7,1) * t222 - Icges(7,4) * t224) * t211;
t145 = Icges(7,6) * t209 + (-Icges(7,4) * t222 - Icges(7,2) * t224) * t211;
t144 = Icges(7,3) * t209 + (-Icges(7,5) * t222 - Icges(7,6) * t224) * t211;
t143 = rSges(4,3) * t210 + t212 * t252;
t142 = -rSges(4,3) * t212 + t210 * t252;
t135 = t200 * V_base(4) - t201 * V_base(5) + V_base(3);
t133 = -rSges(6,1) * t212 + t210 * t250;
t132 = rSges(6,1) * t210 + t212 * t250;
t131 = rSges(5,3) * t210 + t212 * t251;
t130 = -rSges(5,3) * t212 + t210 * t251;
t112 = V_base(5) * rSges(3,3) + (-t179 - t286) * t213 + t257;
t111 = t181 * t213 + (-rSges(3,3) + t281) * V_base(4) + t263;
t110 = V_base(4) * t179 + (-t181 - t285) * V_base(5) + t256;
t109 = rSges(7,1) * t158 + rSges(7,2) * t157 + rSges(7,3) * t272;
t108 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t269;
t107 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t272;
t106 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t269;
t105 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t272;
t104 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t269;
t103 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t272;
t102 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t269;
t101 = t193 * V_base(5) + (-t142 + t259) * t213 + t254;
t100 = t143 * t213 + (-t193 + t281) * V_base(4) + t237;
t99 = V_base(4) * t142 + (-t143 + t258) * V_base(5) + t253;
t98 = t177 * t191 + (-t130 + t255) * t213 + t238;
t97 = t131 * t213 - t177 * t192 + t231;
t96 = t192 * t130 - t191 * t131 + t232;
t95 = t176 * t191 + (-t133 + t239) * t213 + t236;
t94 = t132 * t213 + (-t175 - t176) * t192 + t230;
t93 = t192 * t133 + (-t132 - t152) * t191 + t229;
t92 = t191 * t283 - t109 * t184 + t149 * t153 + (-t161 + t239) * t213 + t236;
t91 = t108 * t184 - t149 * t154 + t160 * t213 + (-t175 - t283) * t192 + t230;
t90 = -t153 * t108 + t154 * t109 + t192 * t161 + (-t152 - t160) * t191 + t229;
t1 = m(2) * (t135 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + t154 * ((t102 * t269 + t155 * t104 + t156 * t106) * t154 + (t103 * t269 + t105 * t155 + t107 * t156) * t153 + (t144 * t269 + t145 * t155 + t146 * t156) * t184) / 0.2e1 + t153 * ((t102 * t272 + t104 * t157 + t106 * t158) * t154 + (t103 * t272 + t157 * t105 + t158 * t107) * t153 + (t144 * t272 + t145 * t157 + t146 * t158) * t184) / 0.2e1 + m(1) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + t184 * ((t102 * t154 + t103 * t153 + t144 * t184) * t209 + ((-t104 * t224 - t106 * t222) * t154 + (-t105 * t224 - t107 * t222) * t153 + (-t145 * t224 - t146 * t222) * t184) * t211) / 0.2e1 + (t291 * t210 - t290 * t212) * t191 / 0.2e1 + (t290 * t210 + t291 * t212) * t192 / 0.2e1 + (t233 * t210 + t226 * t212 + t292 * t213 + (-t170 * t210 + t173 * t212 - t223 * t196 + t198 * t225 + Icges(1,4)) * V_base(5) + (-t210 * t171 + t174 * t212 - t223 * t197 + t199 * t225 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t226 * t210 - t233 * t212 + t293 * t213 + (t170 * t212 + t173 * t210 + t196 * t225 + t223 * t198 + Icges(1,2)) * V_base(5) + (t171 * t212 + t174 * t210 + t197 * t225 + t223 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t297 * t209 + t299 * t211) * t192 + (t298 * t209 + t300 * t211) * t191 + (t138 * t220 + t140 * t219 + t293) * V_base(5) + (t139 * t220 + t141 * t219 + t292) * V_base(4) + (t189 * t220 + t190 * t219 + t295 * t209 - t296 * t211 + Icges(2,3) + Icges(3,3)) * t213) * t213 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
