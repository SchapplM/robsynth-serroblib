% Calculate kinetic energy for
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:42
% EndTime: 2019-12-05 17:14:46
% DurationCPUTime: 3.44s
% Computational Cost: add. (2294->361), mult. (4144->540), div. (0->0), fcn. (4842->12), ass. (0->159)
t261 = cos(pkin(5));
t301 = pkin(6) * t261;
t266 = cos(qJ(3));
t300 = pkin(3) * t266;
t258 = sin(pkin(10));
t298 = Icges(2,4) * t258;
t259 = sin(pkin(5));
t297 = t258 * t259;
t260 = cos(pkin(10));
t296 = t259 * t260;
t263 = sin(qJ(3));
t295 = t259 * t263;
t264 = sin(qJ(2));
t294 = t259 * t264;
t293 = t259 * t266;
t267 = cos(qJ(2));
t292 = t259 * t267;
t291 = t261 * t263;
t290 = t261 * t264;
t289 = t261 * t267;
t288 = qJ(3) + qJ(4);
t287 = qJD(2) * t259;
t286 = V_base(5) * qJ(1) + V_base(1);
t282 = qJD(1) + V_base(3);
t281 = t258 * t295;
t280 = t260 * t295;
t238 = t258 * t287 + V_base(4);
t249 = qJD(2) * t261 + V_base(6);
t279 = cos(t288);
t224 = t258 * t289 + t260 * t264;
t203 = qJD(3) * t224 + t238;
t278 = t259 * t279;
t176 = qJD(4) * t224 + t203;
t237 = -t260 * t287 + V_base(5);
t222 = t258 * t264 - t260 * t289;
t202 = qJD(3) * t222 + t237;
t232 = pkin(1) * t258 - pkin(6) * t296;
t277 = -t232 * V_base(6) + V_base(5) * t301 + t286;
t233 = pkin(1) * t260 + pkin(6) * t297;
t276 = V_base(4) * t232 - t233 * V_base(5) + t282;
t175 = qJD(4) * t222 + t202;
t211 = (-qJD(3) - qJD(4)) * t292 + t249;
t275 = V_base(6) * t233 + V_base(2) + (-qJ(1) - t301) * V_base(4);
t223 = t258 * t267 + t260 * t290;
t193 = t223 * pkin(2) + t222 * pkin(7);
t231 = (pkin(2) * t264 - pkin(7) * t267) * t259;
t274 = -t193 * t249 + t237 * t231 + t277;
t225 = -t258 * t290 + t260 * t267;
t194 = t225 * pkin(2) + t224 * pkin(7);
t273 = t238 * t193 - t194 * t237 + t276;
t272 = t249 * t194 - t231 * t238 + t275;
t151 = -pkin(3) * t280 + pkin(8) * t222 + t223 * t300;
t192 = pkin(3) * t291 + (-pkin(8) * t267 + t264 * t300) * t259;
t226 = -qJD(3) * t292 + t249;
t271 = -t151 * t226 + t202 * t192 + t274;
t152 = pkin(3) * t281 + pkin(8) * t224 + t225 * t300;
t270 = t203 * t151 - t152 * t202 + t273;
t269 = t226 * t152 - t192 * t203 + t272;
t265 = cos(qJ(5));
t262 = sin(qJ(5));
t256 = sin(t288);
t255 = Icges(2,4) * t260;
t247 = rSges(2,1) * t260 - rSges(2,2) * t258;
t246 = rSges(2,1) * t258 + rSges(2,2) * t260;
t245 = Icges(2,1) * t260 - t298;
t244 = Icges(2,1) * t258 + t255;
t243 = -Icges(2,2) * t258 + t255;
t242 = Icges(2,2) * t260 + t298;
t236 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t235 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t234 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t230 = t264 * t293 + t291;
t229 = t261 * t266 - t263 * t294;
t217 = t261 * t256 + t264 * t278;
t216 = t256 * t294 - t261 * t279;
t215 = rSges(3,3) * t261 + (rSges(3,1) * t264 + rSges(3,2) * t267) * t259;
t214 = Icges(3,5) * t261 + (Icges(3,1) * t264 + Icges(3,4) * t267) * t259;
t213 = Icges(3,6) * t261 + (Icges(3,4) * t264 + Icges(3,2) * t267) * t259;
t212 = Icges(3,3) * t261 + (Icges(3,5) * t264 + Icges(3,6) * t267) * t259;
t210 = V_base(5) * rSges(2,3) - t246 * V_base(6) + t286;
t209 = t247 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t207 = t225 * t266 + t281;
t206 = -t225 * t263 + t258 * t293;
t205 = t223 * t266 - t280;
t204 = -t223 * t263 - t260 * t293;
t201 = t246 * V_base(4) - t247 * V_base(5) + t282;
t200 = t217 * t265 - t262 * t292;
t199 = -t217 * t262 - t265 * t292;
t198 = t225 * t279 + t256 * t297;
t197 = t225 * t256 - t258 * t278;
t196 = t223 * t279 - t256 * t296;
t195 = t223 * t256 + t260 * t278;
t191 = pkin(4) * t217 + pkin(9) * t216;
t190 = rSges(4,1) * t230 + rSges(4,2) * t229 - rSges(4,3) * t292;
t189 = Icges(4,1) * t230 + Icges(4,4) * t229 - Icges(4,5) * t292;
t188 = Icges(4,4) * t230 + Icges(4,2) * t229 - Icges(4,6) * t292;
t187 = Icges(4,5) * t230 + Icges(4,6) * t229 - Icges(4,3) * t292;
t186 = rSges(3,1) * t225 - rSges(3,2) * t224 + rSges(3,3) * t297;
t185 = rSges(3,1) * t223 - rSges(3,2) * t222 - rSges(3,3) * t296;
t184 = Icges(3,1) * t225 - Icges(3,4) * t224 + Icges(3,5) * t297;
t183 = Icges(3,1) * t223 - Icges(3,4) * t222 - Icges(3,5) * t296;
t182 = Icges(3,4) * t225 - Icges(3,2) * t224 + Icges(3,6) * t297;
t181 = Icges(3,4) * t223 - Icges(3,2) * t222 - Icges(3,6) * t296;
t180 = Icges(3,5) * t225 - Icges(3,6) * t224 + Icges(3,3) * t297;
t179 = Icges(3,5) * t223 - Icges(3,6) * t222 - Icges(3,3) * t296;
t177 = qJD(5) * t216 + t211;
t173 = rSges(5,1) * t217 - rSges(5,2) * t216 - rSges(5,3) * t292;
t172 = Icges(5,1) * t217 - Icges(5,4) * t216 - Icges(5,5) * t292;
t171 = Icges(5,4) * t217 - Icges(5,2) * t216 - Icges(5,6) * t292;
t170 = Icges(5,5) * t217 - Icges(5,6) * t216 - Icges(5,3) * t292;
t169 = t198 * t265 + t224 * t262;
t168 = -t198 * t262 + t224 * t265;
t167 = t196 * t265 + t222 * t262;
t166 = -t196 * t262 + t222 * t265;
t165 = pkin(4) * t198 + pkin(9) * t197;
t164 = pkin(4) * t196 + pkin(9) * t195;
t162 = rSges(4,1) * t207 + rSges(4,2) * t206 + rSges(4,3) * t224;
t161 = rSges(4,1) * t205 + rSges(4,2) * t204 + rSges(4,3) * t222;
t160 = Icges(4,1) * t207 + Icges(4,4) * t206 + Icges(4,5) * t224;
t159 = Icges(4,1) * t205 + Icges(4,4) * t204 + Icges(4,5) * t222;
t158 = Icges(4,4) * t207 + Icges(4,2) * t206 + Icges(4,6) * t224;
t157 = Icges(4,4) * t205 + Icges(4,2) * t204 + Icges(4,6) * t222;
t156 = Icges(4,5) * t207 + Icges(4,6) * t206 + Icges(4,3) * t224;
t155 = Icges(4,5) * t205 + Icges(4,6) * t204 + Icges(4,3) * t222;
t154 = qJD(5) * t197 + t176;
t153 = qJD(5) * t195 + t175;
t150 = rSges(5,1) * t198 - rSges(5,2) * t197 + rSges(5,3) * t224;
t149 = rSges(5,1) * t196 - rSges(5,2) * t195 + rSges(5,3) * t222;
t148 = Icges(5,1) * t198 - Icges(5,4) * t197 + Icges(5,5) * t224;
t147 = Icges(5,1) * t196 - Icges(5,4) * t195 + Icges(5,5) * t222;
t146 = Icges(5,4) * t198 - Icges(5,2) * t197 + Icges(5,6) * t224;
t145 = Icges(5,4) * t196 - Icges(5,2) * t195 + Icges(5,6) * t222;
t144 = Icges(5,5) * t198 - Icges(5,6) * t197 + Icges(5,3) * t224;
t143 = Icges(5,5) * t196 - Icges(5,6) * t195 + Icges(5,3) * t222;
t142 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t216;
t141 = Icges(6,1) * t200 + Icges(6,4) * t199 + Icges(6,5) * t216;
t140 = Icges(6,4) * t200 + Icges(6,2) * t199 + Icges(6,6) * t216;
t139 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t216;
t137 = -t185 * t249 + t215 * t237 + t277;
t136 = t186 * t249 - t215 * t238 + t275;
t134 = t185 * t238 - t186 * t237 + t276;
t133 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t197;
t132 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t195;
t131 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t197;
t130 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t195;
t129 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t197;
t128 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t195;
t127 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t197;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t195;
t125 = -t161 * t226 + t190 * t202 + t274;
t124 = t162 * t226 - t190 * t203 + t272;
t123 = t161 * t203 - t162 * t202 + t273;
t122 = -t149 * t211 + t173 * t175 + t271;
t121 = t150 * t211 - t173 * t176 + t269;
t120 = t149 * t176 - t150 * t175 + t270;
t119 = -t132 * t177 + t142 * t153 - t164 * t211 + t175 * t191 + t271;
t118 = t133 * t177 - t142 * t154 + t165 * t211 - t176 * t191 + t269;
t117 = t132 * t154 - t133 * t153 + t164 * t176 - t165 * t175 + t270;
t1 = m(1) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + m(2) * (t201 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t238 * ((t180 * t297 - t182 * t224 + t184 * t225) * t238 + (t179 * t297 - t181 * t224 + t183 * t225) * t237 + (t212 * t297 - t213 * t224 + t214 * t225) * t249) / 0.2e1 + t237 * ((-t180 * t296 - t182 * t222 + t184 * t223) * t238 + (-t179 * t296 - t181 * t222 + t183 * t223) * t237 + (-t212 * t296 - t213 * t222 + t214 * t223) * t249) / 0.2e1 + t249 * ((t179 * t237 + t180 * t238 + t212 * t249) * t261 + ((t182 * t267 + t184 * t264) * t238 + (t181 * t267 + t183 * t264) * t237 + (t213 * t267 + t214 * t264) * t249) * t259) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + t203 * ((t156 * t224 + t158 * t206 + t160 * t207) * t203 + (t155 * t224 + t157 * t206 + t159 * t207) * t202 + (t187 * t224 + t188 * t206 + t189 * t207) * t226) / 0.2e1 + t202 * ((t156 * t222 + t158 * t204 + t160 * t205) * t203 + (t155 * t222 + t157 * t204 + t159 * t205) * t202 + (t187 * t222 + t188 * t204 + t189 * t205) * t226) / 0.2e1 + t226 * ((-t156 * t292 + t158 * t229 + t160 * t230) * t203 + (-t155 * t292 + t157 * t229 + t159 * t230) * t202 + (-t187 * t292 + t229 * t188 + t230 * t189) * t226) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + t176 * ((t144 * t224 - t146 * t197 + t148 * t198) * t176 + (t143 * t224 - t145 * t197 + t147 * t198) * t175 + (t170 * t224 - t171 * t197 + t172 * t198) * t211) / 0.2e1 + t175 * ((t144 * t222 - t146 * t195 + t148 * t196) * t176 + (t143 * t222 - t145 * t195 + t147 * t196) * t175 + (t170 * t222 - t171 * t195 + t172 * t196) * t211) / 0.2e1 + t211 * ((-t144 * t292 - t146 * t216 + t148 * t217) * t176 + (-t143 * t292 - t145 * t216 + t147 * t217) * t175 + (-t170 * t292 - t216 * t171 + t217 * t172) * t211) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t154 * ((t197 * t127 + t129 * t168 + t169 * t131) * t154 + (t126 * t197 + t128 * t168 + t130 * t169) * t153 + (t139 * t197 + t140 * t168 + t141 * t169) * t177) / 0.2e1 + t153 * ((t127 * t195 + t129 * t166 + t131 * t167) * t154 + (t126 * t195 + t128 * t166 + t130 * t167) * t153 + (t139 * t195 + t140 * t166 + t141 * t167) * t177) / 0.2e1 + t177 * ((t127 * t216 + t129 * t199 + t131 * t200) * t154 + (t126 * t216 + t128 * t199 + t130 * t200) * t153 + (t139 * t216 + t140 * t199 + t141 * t200) * t177) / 0.2e1 + ((-t242 * t258 + t244 * t260 + Icges(1,4)) * V_base(5) + (-t243 * t258 + t245 * t260 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t242 * t260 + t244 * t258 + Icges(1,2)) * V_base(5) + (t243 * t260 + t245 * t258 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t258 + Icges(2,6) * t260 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t260 - Icges(2,6) * t258 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
