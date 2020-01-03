% Calculate kinetic energy for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR16_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR16_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:38
% EndTime: 2019-12-31 20:44:41
% DurationCPUTime: 3.28s
% Computational Cost: add. (1558->314), mult. (3496->452), div. (0->0), fcn. (4016->10), ass. (0->140)
t294 = Icges(3,1) + Icges(4,2);
t293 = Icges(3,4) + Icges(4,6);
t292 = Icges(3,5) - Icges(4,4);
t291 = Icges(3,2) + Icges(4,3);
t290 = Icges(3,6) - Icges(4,5);
t289 = Icges(3,3) + Icges(4,1);
t244 = cos(pkin(5));
t248 = sin(qJ(1));
t250 = cos(qJ(2));
t267 = t248 * t250;
t247 = sin(qJ(2));
t251 = cos(qJ(1));
t268 = t247 * t251;
t213 = t244 * t267 + t268;
t266 = t250 * t251;
t269 = t247 * t248;
t214 = -t244 * t269 + t266;
t243 = sin(pkin(5));
t272 = t243 * t248;
t288 = t291 * t213 - t293 * t214 - t290 * t272;
t211 = -t244 * t266 + t269;
t212 = t244 * t268 + t267;
t270 = t243 * t251;
t287 = t291 * t211 - t293 * t212 + t290 * t270;
t286 = -t293 * t213 + t294 * t214 + t292 * t272;
t285 = -t293 * t211 + t294 * t212 - t292 * t270;
t284 = -t290 * t213 + t292 * t214 + t289 * t272;
t283 = -t290 * t211 + t292 * t212 - t289 * t270;
t282 = t289 * t244 + (t292 * t247 + t290 * t250) * t243;
t281 = t290 * t244 + (t293 * t247 + t291 * t250) * t243;
t280 = t292 * t244 + (t294 * t247 + t293 * t250) * t243;
t276 = cos(qJ(4));
t275 = pkin(7) * t244;
t274 = Icges(2,4) * t248;
t273 = t243 * t247;
t271 = t243 * t250;
t265 = qJD(2) * t243;
t264 = V_base(5) * pkin(6) + V_base(1);
t261 = t243 * t276;
t224 = t248 * t265 + V_base(4);
t240 = V_base(6) + qJD(1);
t182 = qJD(4) * t214 + t224;
t225 = qJD(2) * t244 + t240;
t207 = qJD(4) * t273 + t225;
t223 = -t251 * t265 + V_base(5);
t218 = t248 * pkin(1) - pkin(7) * t270;
t260 = -t218 * t240 + V_base(5) * t275 + t264;
t219 = pkin(1) * t251 + pkin(7) * t272;
t259 = V_base(4) * t218 - t219 * V_base(5) + V_base(3);
t181 = qJD(4) * t212 + t223;
t215 = (pkin(2) * t247 - qJ(3) * t250) * t243;
t258 = qJD(3) * t213 + t223 * t215 + t260;
t257 = t240 * t219 + V_base(2) + (-pkin(6) - t275) * V_base(4);
t178 = pkin(2) * t214 + qJ(3) * t213;
t256 = qJD(3) * t211 + t225 * t178 + t257;
t177 = pkin(2) * t212 + qJ(3) * t211;
t255 = -qJD(3) * t271 + t224 * t177 + t259;
t192 = -pkin(3) * t270 + t212 * pkin(8);
t217 = pkin(3) * t244 + pkin(8) * t273;
t254 = t223 * t217 + (-t177 - t192) * t225 + t258;
t191 = pkin(3) * t272 + pkin(8) * t214;
t253 = t225 * t191 + (-t215 - t217) * t224 + t256;
t252 = t224 * t192 + (-t178 - t191) * t223 + t255;
t249 = cos(qJ(5));
t246 = sin(qJ(4));
t245 = sin(qJ(5));
t241 = Icges(2,4) * t251;
t233 = rSges(2,1) * t251 - t248 * rSges(2,2);
t232 = t248 * rSges(2,1) + rSges(2,2) * t251;
t231 = Icges(2,1) * t251 - t274;
t230 = Icges(2,1) * t248 + t241;
t229 = -Icges(2,2) * t248 + t241;
t228 = Icges(2,2) * t251 + t274;
t222 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t221 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t220 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t210 = t244 * t276 - t246 * t271;
t209 = t244 * t246 + t250 * t261;
t201 = rSges(4,1) * t244 + (-rSges(4,2) * t247 - rSges(4,3) * t250) * t243;
t200 = rSges(3,3) * t244 + (rSges(3,1) * t247 + rSges(3,2) * t250) * t243;
t190 = V_base(5) * rSges(2,3) - t232 * t240 + t264;
t189 = t233 * t240 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t187 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t186 = t211 * t246 - t251 * t261;
t185 = t211 * t276 + t246 * t270;
t184 = t213 * t246 + t248 * t261;
t183 = -t213 * t276 + t246 * t272;
t180 = t210 * t249 + t245 * t273;
t179 = -t210 * t245 + t249 * t273;
t175 = qJD(5) * t209 + t207;
t174 = pkin(4) * t210 + pkin(9) * t209;
t172 = rSges(3,1) * t214 - rSges(3,2) * t213 + rSges(3,3) * t272;
t171 = t212 * rSges(3,1) - t211 * rSges(3,2) - rSges(3,3) * t270;
t170 = -rSges(4,1) * t270 - t212 * rSges(4,2) + t211 * rSges(4,3);
t169 = rSges(4,1) * t272 - rSges(4,2) * t214 + rSges(4,3) * t213;
t156 = rSges(5,1) * t210 - rSges(5,2) * t209 + rSges(5,3) * t273;
t155 = Icges(5,1) * t210 - Icges(5,4) * t209 + Icges(5,5) * t273;
t154 = Icges(5,4) * t210 - Icges(5,2) * t209 + Icges(5,6) * t273;
t153 = Icges(5,5) * t210 - Icges(5,6) * t209 + Icges(5,3) * t273;
t150 = t186 * t249 + t212 * t245;
t149 = -t186 * t245 + t212 * t249;
t148 = t184 * t249 + t214 * t245;
t147 = -t184 * t245 + t214 * t249;
t146 = qJD(5) * t183 + t182;
t145 = -qJD(5) * t185 + t181;
t144 = pkin(4) * t186 - pkin(9) * t185;
t143 = pkin(4) * t184 + pkin(9) * t183;
t142 = rSges(5,1) * t186 + rSges(5,2) * t185 + rSges(5,3) * t212;
t141 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t214;
t140 = Icges(5,1) * t186 + Icges(5,4) * t185 + Icges(5,5) * t212;
t139 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t214;
t138 = Icges(5,4) * t186 + Icges(5,2) * t185 + Icges(5,6) * t212;
t137 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t214;
t136 = Icges(5,5) * t186 + Icges(5,6) * t185 + Icges(5,3) * t212;
t135 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t214;
t134 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t209;
t133 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t209;
t132 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t209;
t131 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t209;
t130 = -t171 * t225 + t200 * t223 + t260;
t129 = t172 * t225 - t200 * t224 + t257;
t128 = rSges(6,1) * t150 + rSges(6,2) * t149 - rSges(6,3) * t185;
t127 = rSges(6,1) * t148 + rSges(6,2) * t147 + rSges(6,3) * t183;
t126 = Icges(6,1) * t150 + Icges(6,4) * t149 - Icges(6,5) * t185;
t125 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t183;
t124 = Icges(6,4) * t150 + Icges(6,2) * t149 - Icges(6,6) * t185;
t123 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t183;
t122 = Icges(6,5) * t150 + Icges(6,6) * t149 - Icges(6,3) * t185;
t121 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t183;
t120 = t171 * t224 - t172 * t223 + t259;
t119 = t201 * t223 + (-t170 - t177) * t225 + t258;
t118 = t169 * t225 + (-t201 - t215) * t224 + t256;
t117 = t170 * t224 + (-t169 - t178) * t223 + t255;
t116 = -t142 * t207 + t156 * t181 + t254;
t115 = t141 * t207 - t156 * t182 + t253;
t114 = -t141 * t181 + t142 * t182 + t252;
t113 = -t128 * t175 + t134 * t145 - t144 * t207 + t174 * t181 + t254;
t112 = t127 * t175 - t134 * t146 + t143 * t207 - t174 * t182 + t253;
t111 = -t127 * t145 + t128 * t146 - t143 * t181 + t144 * t182 + t252;
t1 = m(1) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(2) * (t187 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + t182 * ((t135 * t214 - t137 * t183 + t139 * t184) * t182 + (t136 * t214 - t138 * t183 + t140 * t184) * t181 + (t153 * t214 - t154 * t183 + t155 * t184) * t207) / 0.2e1 + t181 * ((t135 * t212 + t137 * t185 + t139 * t186) * t182 + (t136 * t212 + t138 * t185 + t140 * t186) * t181 + (t153 * t212 + t154 * t185 + t155 * t186) * t207) / 0.2e1 + t207 * ((t135 * t273 - t137 * t209 + t139 * t210) * t182 + (t136 * t273 - t138 * t209 + t140 * t210) * t181 + (t153 * t273 - t154 * t209 + t155 * t210) * t207) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + t146 * ((t183 * t121 + t147 * t123 + t148 * t125) * t146 + (t122 * t183 + t124 * t147 + t126 * t148) * t145 + (t131 * t183 + t132 * t147 + t133 * t148) * t175) / 0.2e1 + t145 * ((-t121 * t185 + t123 * t149 + t125 * t150) * t146 + (-t185 * t122 + t149 * t124 + t150 * t126) * t145 + (-t131 * t185 + t132 * t149 + t133 * t150) * t175) / 0.2e1 + t175 * ((t121 * t209 + t123 * t179 + t125 * t180) * t146 + (t122 * t209 + t124 * t179 + t126 * t180) * t145 + (t131 * t209 + t132 * t179 + t133 * t180) * t175) / 0.2e1 + ((-t211 * t281 + t212 * t280 - t270 * t282) * t225 + (t211 * t288 + t286 * t212 - t284 * t270) * t224 + (t287 * t211 + t285 * t212 - t283 * t270) * t223) * t223 / 0.2e1 + ((-t213 * t281 + t214 * t280 + t272 * t282) * t225 + (t213 * t288 + t286 * t214 + t284 * t272) * t224 + (t213 * t287 + t214 * t285 + t272 * t283) * t223) * t224 / 0.2e1 + ((t223 * t283 + t224 * t284 + t225 * t282) * t244 + ((t247 * t280 + t250 * t281) * t225 + (t286 * t247 - t250 * t288) * t224 + (t247 * t285 - t250 * t287) * t223) * t243) * t225 / 0.2e1 + ((-t248 * t228 + t230 * t251 + Icges(1,4)) * V_base(5) + (-t248 * t229 + t231 * t251 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t251 + t248 * t230 + Icges(1,2)) * V_base(5) + (t229 * t251 + t248 * t231 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t248 + Icges(2,6) * t251) * V_base(5) + (Icges(2,5) * t251 - Icges(2,6) * t248) * V_base(4) + Icges(2,3) * t240 / 0.2e1) * t240;
T = t1;
