% Calculate kinetic energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:07
% EndTime: 2024-09-27 21:45:09
% DurationCPUTime: 0.39s
% Computational Cost: add. (2750->340), mult. (1978->475), div. (0->0), fcn. (1752->22), ass. (0->166)
t271 = pkin(7) + pkin(8);
t266 = cos(pkin(10));
t302 = pkin(1) * t266;
t267 = cos(pkin(5));
t301 = t267 * pkin(7);
t270 = cos(qJ(3));
t300 = t270 * pkin(3);
t299 = -pkin(6) - qJ(1);
t264 = sin(pkin(10));
t298 = Icges(2,4) * t264;
t261 = pkin(10) + qJ(2);
t249 = sin(t261);
t297 = Icges(3,4) * t249;
t265 = sin(pkin(5));
t296 = t249 * t265;
t250 = cos(t261);
t295 = t250 * t265;
t269 = sin(qJ(3));
t294 = t267 * t269;
t293 = t267 * t270;
t263 = qJ(3) + qJ(4);
t209 = pkin(3) * t294 - t265 * t271;
t245 = cos(qJ(4)) * pkin(4) + pkin(3);
t262 = pkin(9) + t271;
t286 = pkin(4) * sin(qJ(4)) * t270;
t292 = -t262 * t265 + (t245 * t269 + t286) * t267 - t209;
t256 = cos(t263);
t291 = pkin(4) * t256;
t290 = qJD(3) * t265;
t289 = -qJD(3) - qJD(4);
t281 = pkin(1) * V_base(6);
t288 = t266 * t281 + V_base(2);
t287 = V_base(5) * qJ(1) + V_base(1);
t282 = qJD(1) + V_base(3);
t213 = t249 * t290 + V_base(4);
t254 = V_base(6) + qJD(2);
t253 = pkin(5) - t263;
t252 = pkin(5) + t263;
t280 = pkin(7) * t265 + t209;
t279 = V_base(4) * t264 * pkin(1) + t282;
t179 = qJD(4) * t296 + t213;
t228 = qJD(3) * t267 + t254;
t198 = qJD(4) * t267 + t228;
t278 = V_base(5) * pkin(6) - t264 * t281 + t287;
t197 = pkin(2) * t250 + pkin(7) * t296;
t277 = t254 * t197 + (t299 - t301) * V_base(4) + t288;
t196 = pkin(2) * t249 - pkin(7) * t295;
t276 = V_base(4) * t196 + (-t197 - t302) * V_base(5) + t279;
t275 = -t196 * t254 + V_base(5) * t301 + t278;
t151 = -t249 * t280 + t250 * t300;
t191 = t265 * t269 * pkin(3) + (-pkin(7) + t271) * t267;
t274 = t228 * t151 - t191 * t213 + t277;
t150 = t249 * t300 + t250 * t280;
t212 = -t250 * t290 + V_base(5);
t273 = t213 * t150 - t151 * t212 + t276;
t272 = -t150 * t228 + t212 * t191 + t275;
t258 = qJ(5) + t263;
t255 = sin(t263);
t251 = Icges(2,4) * t266;
t244 = cos(t258);
t243 = sin(t258);
t242 = -qJ(5) + t253;
t241 = qJ(5) + t252;
t240 = cos(t252);
t239 = sin(t253);
t238 = Icges(3,4) * t250;
t234 = cos(t241);
t233 = sin(t242);
t232 = cos(t253) / 0.2e1;
t231 = sin(t252) / 0.2e1;
t230 = cos(t242) / 0.2e1;
t229 = sin(t241) / 0.2e1;
t225 = rSges(2,1) * t266 - rSges(2,2) * t264;
t224 = rSges(2,1) * t264 + rSges(2,2) * t266;
t223 = Icges(2,1) * t266 - t298;
t222 = Icges(2,1) * t264 + t251;
t221 = -Icges(2,2) * t264 + t251;
t220 = Icges(2,2) * t266 + t298;
t217 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t216 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t215 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t211 = rSges(3,1) * t250 - rSges(3,2) * t249;
t210 = rSges(3,1) * t249 + rSges(3,2) * t250;
t208 = Icges(3,1) * t250 - t297;
t207 = Icges(3,1) * t249 + t238;
t206 = -Icges(3,2) * t249 + t238;
t205 = Icges(3,2) * t250 + t297;
t202 = t232 - t240 / 0.2e1;
t201 = t232 + t240 / 0.2e1;
t200 = t231 - t239 / 0.2e1;
t199 = t231 + t239 / 0.2e1;
t195 = t230 - t234 / 0.2e1;
t194 = t230 + t234 / 0.2e1;
t193 = t229 - t233 / 0.2e1;
t192 = t229 + t233 / 0.2e1;
t189 = qJD(5) * t267 + t198;
t188 = -t249 * t294 + t250 * t270;
t187 = -t249 * t293 - t250 * t269;
t186 = t249 * t270 + t250 * t294;
t185 = -t249 * t269 + t250 * t293;
t184 = rSges(4,3) * t267 + (rSges(4,1) * t269 + rSges(4,2) * t270) * t265;
t183 = Icges(4,5) * t267 + (Icges(4,1) * t269 + Icges(4,4) * t270) * t265;
t182 = Icges(4,6) * t267 + (Icges(4,4) * t269 + Icges(4,2) * t270) * t265;
t181 = Icges(4,3) * t267 + (Icges(4,5) * t269 + Icges(4,6) * t270) * t265;
t178 = t289 * t295 + V_base(5);
t177 = V_base(5) * rSges(2,3) - t224 * V_base(6) + t287;
t176 = t225 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t174 = -t200 * t249 + t250 * t256;
t173 = -t201 * t249 - t250 * t255;
t172 = t200 * t250 + t249 * t256;
t171 = t201 * t250 - t249 * t255;
t170 = t224 * V_base(4) - t225 * V_base(5) + t282;
t168 = qJD(5) * t296 + t179;
t167 = V_base(5) + (-qJD(5) + t289) * t295;
t166 = -t193 * t249 + t244 * t250;
t165 = -t194 * t249 - t243 * t250;
t164 = t193 * t250 + t244 * t249;
t163 = t194 * t250 - t243 * t249;
t162 = rSges(5,1) * t202 + rSges(5,2) * t199 + rSges(5,3) * t267;
t161 = Icges(5,1) * t202 + Icges(5,4) * t199 + Icges(5,5) * t267;
t160 = Icges(5,4) * t202 + Icges(5,2) * t199 + Icges(5,6) * t267;
t159 = Icges(5,5) * t202 + Icges(5,6) * t199 + Icges(5,3) * t267;
t158 = V_base(5) * rSges(3,3) - t210 * t254 + t278;
t157 = t211 * t254 + (-rSges(3,3) + t299) * V_base(4) + t288;
t156 = rSges(6,1) * t195 + rSges(6,2) * t192 + rSges(6,3) * t267;
t155 = (t262 - t271) * t267 + (t286 + (-pkin(3) + t245) * t269) * t265;
t154 = Icges(6,1) * t195 + Icges(6,4) * t192 + Icges(6,5) * t267;
t153 = Icges(6,4) * t195 + Icges(6,2) * t192 + Icges(6,6) * t267;
t152 = Icges(6,5) * t195 + Icges(6,6) * t192 + Icges(6,3) * t267;
t149 = t210 * V_base(4) + (-t211 - t302) * V_base(5) + t279;
t148 = rSges(4,1) * t188 + rSges(4,2) * t187 + rSges(4,3) * t296;
t147 = rSges(4,1) * t186 + rSges(4,2) * t185 - rSges(4,3) * t295;
t146 = Icges(4,1) * t188 + Icges(4,4) * t187 + Icges(4,5) * t296;
t145 = Icges(4,1) * t186 + Icges(4,4) * t185 - Icges(4,5) * t295;
t144 = Icges(4,4) * t188 + Icges(4,2) * t187 + Icges(4,6) * t296;
t143 = Icges(4,4) * t186 + Icges(4,2) * t185 - Icges(4,6) * t295;
t142 = Icges(4,5) * t188 + Icges(4,6) * t187 + Icges(4,3) * t296;
t141 = Icges(4,5) * t186 + Icges(4,6) * t185 - Icges(4,3) * t295;
t138 = rSges(5,1) * t174 + rSges(5,2) * t173 + rSges(5,3) * t296;
t137 = rSges(5,1) * t172 + rSges(5,2) * t171 - rSges(5,3) * t295;
t136 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t296;
t135 = Icges(5,1) * t172 + Icges(5,4) * t171 - Icges(5,5) * t295;
t134 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t296;
t133 = Icges(5,4) * t172 + Icges(5,2) * t171 - Icges(5,6) * t295;
t132 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t296;
t131 = Icges(5,5) * t172 + Icges(5,6) * t171 - Icges(5,3) * t295;
t130 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t296;
t129 = rSges(6,1) * t164 + rSges(6,2) * t163 - rSges(6,3) * t295;
t128 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t296;
t127 = Icges(6,1) * t164 + Icges(6,4) * t163 - Icges(6,5) * t295;
t126 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t296;
t125 = Icges(6,4) * t164 + Icges(6,2) * t163 - Icges(6,6) * t295;
t124 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t296;
t123 = Icges(6,5) * t164 + Icges(6,6) * t163 - Icges(6,3) * t295;
t122 = -t249 * t292 + t250 * t291;
t121 = t249 * t291 + t250 * t292;
t120 = -t147 * t228 + t184 * t212 + t275;
t119 = t148 * t228 - t184 * t213 + t277;
t118 = t147 * t213 - t148 * t212 + t276;
t117 = -t137 * t198 + t162 * t178 + t272;
t116 = t138 * t198 - t162 * t179 + t274;
t115 = t137 * t179 - t138 * t178 + t273;
t114 = -t121 * t198 - t129 * t189 + t155 * t178 + t156 * t167 + t272;
t113 = t122 * t198 + t130 * t189 - t155 * t179 - t156 * t168 + t274;
t112 = t121 * t179 - t122 * t178 + t129 * t168 - t130 * t167 + t273;
t1 = m(1) * (t215 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + m(2) * (t170 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(3) * (t149 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + t213 * ((t142 * t296 + t144 * t187 + t146 * t188) * t213 + (t141 * t296 + t143 * t187 + t145 * t188) * t212 + (t181 * t296 + t182 * t187 + t183 * t188) * t228) / 0.2e1 + t212 * ((-t142 * t295 + t144 * t185 + t146 * t186) * t213 + (-t141 * t295 + t143 * t185 + t145 * t186) * t212 + (-t181 * t295 + t182 * t185 + t183 * t186) * t228) / 0.2e1 + t228 * ((t141 * t212 + t142 * t213 + t181 * t228) * t267 + ((t144 * t270 + t146 * t269) * t213 + (t143 * t270 + t145 * t269) * t212 + (t182 * t270 + t183 * t269) * t228) * t265) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + t179 * ((t132 * t296 + t173 * t134 + t174 * t136) * t179 + (t131 * t296 + t133 * t173 + t135 * t174) * t178 + (t159 * t296 + t160 * t173 + t161 * t174) * t198) / 0.2e1 + t178 * ((-t132 * t295 + t134 * t171 + t136 * t172) * t179 + (-t131 * t295 + t171 * t133 + t172 * t135) * t178 + (-t159 * t295 + t160 * t171 + t161 * t172) * t198) / 0.2e1 + t198 * ((t132 * t267 + t134 * t199 + t136 * t202) * t179 + (t131 * t267 + t133 * t199 + t135 * t202) * t178 + (t159 * t267 + t160 * t199 + t161 * t202) * t198) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t168 * ((t124 * t296 + t165 * t126 + t166 * t128) * t168 + (t123 * t296 + t125 * t165 + t127 * t166) * t167 + (t152 * t296 + t153 * t165 + t154 * t166) * t189) / 0.2e1 + t167 * ((-t124 * t295 + t126 * t163 + t128 * t164) * t168 + (-t123 * t295 + t163 * t125 + t164 * t127) * t167 + (-t152 * t295 + t153 * t163 + t154 * t164) * t189) / 0.2e1 + t189 * ((t124 * t267 + t126 * t192 + t128 * t195) * t168 + (t123 * t267 + t125 * t192 + t127 * t195) * t167 + (t152 * t267 + t153 * t192 + t154 * t195) * t189) / 0.2e1 + ((-t205 * t249 + t207 * t250 - t220 * t264 + t222 * t266 + Icges(1,4)) * V_base(5) + (-t206 * t249 + t208 * t250 - t221 * t264 + t223 * t266 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t205 * t250 + t249 * t207 + t220 * t266 + t222 * t264 + Icges(1,2)) * V_base(5) + (t206 * t250 + t208 * t249 + t221 * t266 + t223 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t264 + Icges(2,6) * t266 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t266 - Icges(2,6) * t264 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t249 + Icges(3,6) * t250) * V_base(5) + (Icges(3,5) * t250 - Icges(3,6) * t249) * V_base(4) + Icges(3,3) * t254 / 0.2e1) * t254;
T = t1;
