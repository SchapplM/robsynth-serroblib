% Calculate kinetic energy for
% S5RRRRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:54
% EndTime: 2019-12-31 22:31:57
% DurationCPUTime: 3.18s
% Computational Cost: add. (2339->361), mult. (4144->541), div. (0->0), fcn. (4842->12), ass. (0->160)
t260 = cos(pkin(5));
t301 = pkin(7) * t260;
t266 = cos(qJ(3));
t300 = pkin(3) * t266;
t264 = sin(qJ(1));
t298 = Icges(2,4) * t264;
t259 = sin(pkin(5));
t263 = sin(qJ(2));
t297 = t259 * t263;
t296 = t259 * t264;
t295 = t259 * t266;
t267 = cos(qJ(2));
t294 = t259 * t267;
t268 = cos(qJ(1));
t293 = t259 * t268;
t262 = sin(qJ(3));
t292 = t260 * t262;
t291 = t263 * t264;
t290 = t263 * t268;
t289 = t264 * t267;
t288 = t267 * t268;
t287 = qJ(3) + qJ(4);
t286 = qJD(2) * t259;
t285 = V_base(5) * pkin(6) + V_base(1);
t282 = t262 * t296;
t281 = t262 * t293;
t238 = t264 * t286 + V_base(4);
t280 = cos(t287);
t255 = V_base(6) + qJD(1);
t228 = t260 * t289 + t290;
t202 = qJD(3) * t228 + t238;
t240 = qJD(2) * t260 + t255;
t279 = t259 * t280;
t177 = qJD(4) * t228 + t202;
t237 = -t268 * t286 + V_base(5);
t232 = pkin(1) * t264 - pkin(7) * t293;
t278 = -t232 * t255 + V_base(5) * t301 + t285;
t233 = pkin(1) * t268 + pkin(7) * t296;
t277 = V_base(4) * t232 - t233 * V_base(5) + V_base(3);
t226 = -t260 * t288 + t291;
t201 = qJD(3) * t226 + t237;
t176 = qJD(4) * t226 + t201;
t276 = t255 * t233 + V_base(2) + (-pkin(6) - t301) * V_base(4);
t211 = (-qJD(3) - qJD(4)) * t294 + t240;
t227 = t260 * t290 + t289;
t193 = t227 * pkin(2) + t226 * pkin(8);
t231 = (pkin(2) * t263 - pkin(8) * t267) * t259;
t275 = -t193 * t240 + t237 * t231 + t278;
t229 = -t260 * t291 + t288;
t194 = t229 * pkin(2) + t228 * pkin(8);
t274 = t238 * t193 - t194 * t237 + t277;
t273 = t240 * t194 - t231 * t238 + t276;
t151 = -pkin(3) * t281 + pkin(9) * t226 + t227 * t300;
t190 = pkin(3) * t292 + (-pkin(9) * t267 + t263 * t300) * t259;
t222 = -qJD(3) * t294 + t240;
t272 = -t151 * t222 + t201 * t190 + t275;
t152 = pkin(3) * t282 + pkin(9) * t228 + t229 * t300;
t271 = t202 * t151 - t152 * t201 + t274;
t270 = t222 * t152 - t190 * t202 + t273;
t265 = cos(qJ(5));
t261 = sin(qJ(5));
t257 = Icges(2,4) * t268;
t256 = sin(t287);
t248 = rSges(2,1) * t268 - rSges(2,2) * t264;
t247 = rSges(2,1) * t264 + rSges(2,2) * t268;
t246 = Icges(2,1) * t268 - t298;
t245 = Icges(2,1) * t264 + t257;
t244 = -Icges(2,2) * t264 + t257;
t243 = Icges(2,2) * t268 + t298;
t236 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t235 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t234 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t225 = t263 * t295 + t292;
t224 = t260 * t266 - t262 * t297;
t217 = t260 * t256 + t263 * t279;
t216 = t256 * t297 - t260 * t280;
t215 = rSges(3,3) * t260 + (rSges(3,1) * t263 + rSges(3,2) * t267) * t259;
t214 = Icges(3,5) * t260 + (Icges(3,1) * t263 + Icges(3,4) * t267) * t259;
t213 = Icges(3,6) * t260 + (Icges(3,4) * t263 + Icges(3,2) * t267) * t259;
t212 = Icges(3,3) * t260 + (Icges(3,5) * t263 + Icges(3,6) * t267) * t259;
t210 = V_base(5) * rSges(2,3) - t247 * t255 + t285;
t209 = t248 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t207 = t247 * V_base(4) - t248 * V_base(5) + V_base(3);
t206 = t229 * t266 + t282;
t205 = -t229 * t262 + t264 * t295;
t204 = t227 * t266 - t281;
t203 = -t227 * t262 - t266 * t293;
t200 = t229 * t280 + t256 * t296;
t199 = t229 * t256 - t264 * t279;
t198 = t227 * t280 - t256 * t293;
t197 = t227 * t256 + t268 * t279;
t196 = t217 * t265 - t261 * t294;
t195 = -t217 * t261 - t265 * t294;
t192 = rSges(3,1) * t229 - rSges(3,2) * t228 + rSges(3,3) * t296;
t191 = rSges(3,1) * t227 - rSges(3,2) * t226 - rSges(3,3) * t293;
t189 = Icges(3,1) * t229 - Icges(3,4) * t228 + Icges(3,5) * t296;
t188 = Icges(3,1) * t227 - Icges(3,4) * t226 - Icges(3,5) * t293;
t187 = Icges(3,4) * t229 - Icges(3,2) * t228 + Icges(3,6) * t296;
t186 = Icges(3,4) * t227 - Icges(3,2) * t226 - Icges(3,6) * t293;
t185 = Icges(3,5) * t229 - Icges(3,6) * t228 + Icges(3,3) * t296;
t184 = Icges(3,5) * t227 - Icges(3,6) * t226 - Icges(3,3) * t293;
t183 = pkin(4) * t217 + pkin(10) * t216;
t182 = rSges(4,1) * t225 + rSges(4,2) * t224 - rSges(4,3) * t294;
t181 = Icges(4,1) * t225 + Icges(4,4) * t224 - Icges(4,5) * t294;
t180 = Icges(4,4) * t225 + Icges(4,2) * t224 - Icges(4,6) * t294;
t179 = Icges(4,5) * t225 + Icges(4,6) * t224 - Icges(4,3) * t294;
t174 = qJD(5) * t216 + t211;
t173 = rSges(5,1) * t217 - rSges(5,2) * t216 - rSges(5,3) * t294;
t172 = Icges(5,1) * t217 - Icges(5,4) * t216 - Icges(5,5) * t294;
t171 = Icges(5,4) * t217 - Icges(5,2) * t216 - Icges(5,6) * t294;
t170 = Icges(5,5) * t217 - Icges(5,6) * t216 - Icges(5,3) * t294;
t169 = t200 * t265 + t228 * t261;
t168 = -t200 * t261 + t228 * t265;
t167 = t198 * t265 + t226 * t261;
t166 = -t198 * t261 + t226 * t265;
t165 = pkin(4) * t200 + pkin(10) * t199;
t164 = pkin(4) * t198 + pkin(10) * t197;
t162 = rSges(4,1) * t206 + rSges(4,2) * t205 + rSges(4,3) * t228;
t161 = rSges(4,1) * t204 + rSges(4,2) * t203 + rSges(4,3) * t226;
t160 = Icges(4,1) * t206 + Icges(4,4) * t205 + Icges(4,5) * t228;
t159 = Icges(4,1) * t204 + Icges(4,4) * t203 + Icges(4,5) * t226;
t158 = Icges(4,4) * t206 + Icges(4,2) * t205 + Icges(4,6) * t228;
t157 = Icges(4,4) * t204 + Icges(4,2) * t203 + Icges(4,6) * t226;
t156 = Icges(4,5) * t206 + Icges(4,6) * t205 + Icges(4,3) * t228;
t155 = Icges(4,5) * t204 + Icges(4,6) * t203 + Icges(4,3) * t226;
t154 = qJD(5) * t199 + t177;
t153 = qJD(5) * t197 + t176;
t150 = rSges(5,1) * t200 - rSges(5,2) * t199 + rSges(5,3) * t228;
t149 = rSges(5,1) * t198 - rSges(5,2) * t197 + rSges(5,3) * t226;
t148 = Icges(5,1) * t200 - Icges(5,4) * t199 + Icges(5,5) * t228;
t147 = Icges(5,1) * t198 - Icges(5,4) * t197 + Icges(5,5) * t226;
t146 = Icges(5,4) * t200 - Icges(5,2) * t199 + Icges(5,6) * t228;
t145 = Icges(5,4) * t198 - Icges(5,2) * t197 + Icges(5,6) * t226;
t144 = Icges(5,5) * t200 - Icges(5,6) * t199 + Icges(5,3) * t228;
t143 = Icges(5,5) * t198 - Icges(5,6) * t197 + Icges(5,3) * t226;
t142 = rSges(6,1) * t196 + rSges(6,2) * t195 + rSges(6,3) * t216;
t141 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t216;
t140 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t216;
t139 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t216;
t136 = -t191 * t240 + t215 * t237 + t278;
t135 = t192 * t240 - t215 * t238 + t276;
t134 = t191 * t238 - t192 * t237 + t277;
t133 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t199;
t132 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t197;
t131 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t199;
t130 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t197;
t129 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t199;
t128 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t197;
t127 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t199;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t197;
t125 = -t161 * t222 + t182 * t201 + t275;
t124 = t162 * t222 - t182 * t202 + t273;
t123 = t161 * t202 - t162 * t201 + t274;
t122 = -t149 * t211 + t173 * t176 + t272;
t121 = t150 * t211 - t173 * t177 + t270;
t120 = t149 * t177 - t150 * t176 + t271;
t119 = -t132 * t174 + t142 * t153 - t164 * t211 + t176 * t183 + t272;
t118 = t133 * t174 - t142 * t154 + t165 * t211 - t177 * t183 + t270;
t117 = t132 * t154 - t133 * t153 + t164 * t177 - t165 * t176 + t271;
t1 = m(1) * (t234 ^ 2 + t235 ^ 2 + t236 ^ 2) / 0.2e1 + m(2) * (t207 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t238 * ((t185 * t296 - t187 * t228 + t189 * t229) * t238 + (t184 * t296 - t186 * t228 + t188 * t229) * t237 + (t212 * t296 - t213 * t228 + t214 * t229) * t240) / 0.2e1 + t237 * ((-t185 * t293 - t187 * t226 + t189 * t227) * t238 + (-t184 * t293 - t186 * t226 + t188 * t227) * t237 + (-t212 * t293 - t213 * t226 + t214 * t227) * t240) / 0.2e1 + t240 * ((t184 * t237 + t185 * t238 + t212 * t240) * t260 + ((t187 * t267 + t189 * t263) * t238 + (t186 * t267 + t188 * t263) * t237 + (t213 * t267 + t214 * t263) * t240) * t259) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + t202 * ((t156 * t228 + t158 * t205 + t160 * t206) * t202 + (t155 * t228 + t157 * t205 + t159 * t206) * t201 + (t179 * t228 + t180 * t205 + t181 * t206) * t222) / 0.2e1 + t201 * ((t156 * t226 + t158 * t203 + t160 * t204) * t202 + (t155 * t226 + t157 * t203 + t159 * t204) * t201 + (t179 * t226 + t180 * t203 + t181 * t204) * t222) / 0.2e1 + t222 * ((-t156 * t294 + t158 * t224 + t160 * t225) * t202 + (-t155 * t294 + t157 * t224 + t159 * t225) * t201 + (-t179 * t294 + t180 * t224 + t181 * t225) * t222) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + t177 * ((t144 * t228 - t146 * t199 + t148 * t200) * t177 + (t143 * t228 - t145 * t199 + t147 * t200) * t176 + (t170 * t228 - t171 * t199 + t172 * t200) * t211) / 0.2e1 + t176 * ((t144 * t226 - t146 * t197 + t148 * t198) * t177 + (t143 * t226 - t145 * t197 + t147 * t198) * t176 + (t170 * t226 - t171 * t197 + t172 * t198) * t211) / 0.2e1 + t211 * ((-t144 * t294 - t146 * t216 + t148 * t217) * t177 + (-t143 * t294 - t145 * t216 + t147 * t217) * t176 + (-t170 * t294 - t171 * t216 + t172 * t217) * t211) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t154 * ((t127 * t199 + t129 * t168 + t131 * t169) * t154 + (t126 * t199 + t128 * t168 + t130 * t169) * t153 + (t139 * t199 + t140 * t168 + t141 * t169) * t174) / 0.2e1 + t153 * ((t127 * t197 + t129 * t166 + t131 * t167) * t154 + (t126 * t197 + t128 * t166 + t130 * t167) * t153 + (t139 * t197 + t140 * t166 + t141 * t167) * t174) / 0.2e1 + t174 * ((t127 * t216 + t129 * t195 + t131 * t196) * t154 + (t126 * t216 + t128 * t195 + t130 * t196) * t153 + (t139 * t216 + t140 * t195 + t141 * t196) * t174) / 0.2e1 + ((-t243 * t264 + t245 * t268 + Icges(1,4)) * V_base(5) + (-t244 * t264 + t246 * t268 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t243 * t268 + t245 * t264 + Icges(1,2)) * V_base(5) + (t244 * t268 + t246 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t264 + Icges(2,6) * t268) * V_base(5) + (Icges(2,5) * t268 - Icges(2,6) * t264) * V_base(4) + Icges(2,3) * t255 / 0.2e1) * t255;
T = t1;
