% Calculate kinetic energy for
% S5RRPRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR14_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:41
% EndTime: 2019-12-31 20:35:45
% DurationCPUTime: 3.21s
% Computational Cost: add. (2234->360), mult. (3934->519), div. (0->0), fcn. (4590->12), ass. (0->159)
t259 = cos(pkin(5));
t263 = sin(qJ(1));
t265 = cos(qJ(2));
t287 = t263 * t265;
t262 = sin(qJ(2));
t266 = cos(qJ(1));
t288 = t262 * t266;
t224 = t259 * t288 + t287;
t256 = sin(pkin(10));
t258 = cos(pkin(10));
t257 = sin(pkin(5));
t290 = t257 * t266;
t199 = -t224 * t256 - t258 * t290;
t278 = t256 * t290;
t200 = t224 * t258 - t278;
t286 = t265 * t266;
t289 = t262 * t263;
t223 = -t259 * t286 + t289;
t153 = Icges(4,5) * t200 + Icges(4,6) * t199 + Icges(4,3) * t223;
t184 = Icges(3,4) * t224 - Icges(3,2) * t223 - Icges(3,6) * t290;
t303 = t153 - t184;
t226 = -t259 * t289 + t286;
t292 = t257 * t263;
t201 = -t226 * t256 + t258 * t292;
t279 = t256 * t292;
t202 = t226 * t258 + t279;
t225 = t259 * t287 + t288;
t154 = Icges(4,5) * t202 + Icges(4,6) * t201 + Icges(4,3) * t225;
t185 = Icges(3,4) * t226 - Icges(3,2) * t225 + Icges(3,6) * t292;
t302 = t154 - t185;
t293 = t257 * t262;
t221 = -t256 * t293 + t258 * t259;
t294 = t256 * t259;
t222 = t258 * t293 + t294;
t291 = t257 * t265;
t176 = Icges(4,5) * t222 + Icges(4,6) * t221 - Icges(4,3) * t291;
t210 = Icges(3,6) * t259 + (Icges(3,4) * t262 + Icges(3,2) * t265) * t257;
t301 = t176 - t210;
t297 = pkin(7) * t259;
t296 = pkin(3) * t258;
t295 = Icges(2,4) * t263;
t284 = qJD(2) * t257;
t283 = pkin(10) + qJ(4);
t282 = V_base(5) * pkin(6) + V_base(1);
t236 = t263 * t284 + V_base(4);
t253 = V_base(6) + qJD(1);
t277 = cos(t283);
t204 = qJD(4) * t225 + t236;
t237 = qJD(2) * t259 + t253;
t276 = t257 * t277;
t235 = -t266 * t284 + V_base(5);
t229 = t263 * pkin(1) - pkin(7) * t290;
t275 = -t229 * t253 + V_base(5) * t297 + t282;
t230 = pkin(1) * t266 + pkin(7) * t292;
t274 = V_base(4) * t229 - t230 * V_base(5) + V_base(3);
t203 = qJD(4) * t223 + t235;
t219 = -qJD(4) * t291 + t237;
t227 = (pkin(2) * t262 - qJ(3) * t265) * t257;
t273 = qJD(3) * t225 + t235 * t227 + t275;
t272 = t253 * t230 + V_base(2) + (-pkin(6) - t297) * V_base(4);
t192 = pkin(2) * t226 + qJ(3) * t225;
t271 = qJD(3) * t223 + t237 * t192 + t272;
t191 = t224 * pkin(2) + t223 * qJ(3);
t270 = -qJD(3) * t291 + t236 * t191 + t274;
t151 = -pkin(3) * t278 + pkin(8) * t223 + t224 * t296;
t181 = pkin(3) * t294 + (-pkin(8) * t265 + t262 * t296) * t257;
t269 = t235 * t181 + (-t151 - t191) * t237 + t273;
t152 = pkin(3) * t279 + pkin(8) * t225 + t226 * t296;
t268 = t237 * t152 + (-t181 - t227) * t236 + t271;
t267 = t236 * t151 + (-t152 - t192) * t235 + t270;
t264 = cos(qJ(5));
t261 = sin(qJ(5));
t254 = Icges(2,4) * t266;
t252 = sin(t283);
t245 = rSges(2,1) * t266 - t263 * rSges(2,2);
t244 = t263 * rSges(2,1) + rSges(2,2) * t266;
t243 = Icges(2,1) * t266 - t295;
t242 = Icges(2,1) * t263 + t254;
t241 = -Icges(2,2) * t263 + t254;
t240 = Icges(2,2) * t266 + t295;
t233 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t232 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t231 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = rSges(3,3) * t259 + (rSges(3,1) * t262 + rSges(3,2) * t265) * t257;
t213 = t259 * t252 + t262 * t276;
t212 = t252 * t293 - t259 * t277;
t211 = Icges(3,5) * t259 + (Icges(3,1) * t262 + Icges(3,4) * t265) * t257;
t209 = Icges(3,3) * t259 + (Icges(3,5) * t262 + Icges(3,6) * t265) * t257;
t208 = V_base(5) * rSges(2,3) - t244 * t253 + t282;
t207 = t245 * t253 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t205 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t198 = t226 * t277 + t252 * t292;
t197 = t226 * t252 - t263 * t276;
t196 = t224 * t277 - t252 * t290;
t195 = t224 * t252 + t266 * t276;
t194 = t213 * t264 - t261 * t291;
t193 = -t213 * t261 - t264 * t291;
t190 = qJD(5) * t212 + t219;
t189 = rSges(3,1) * t226 - rSges(3,2) * t225 + rSges(3,3) * t292;
t188 = t224 * rSges(3,1) - t223 * rSges(3,2) - rSges(3,3) * t290;
t187 = Icges(3,1) * t226 - Icges(3,4) * t225 + Icges(3,5) * t292;
t186 = Icges(3,1) * t224 - Icges(3,4) * t223 - Icges(3,5) * t290;
t183 = Icges(3,5) * t226 - Icges(3,6) * t225 + Icges(3,3) * t292;
t182 = Icges(3,5) * t224 - Icges(3,6) * t223 - Icges(3,3) * t290;
t180 = pkin(4) * t213 + pkin(9) * t212;
t179 = rSges(4,1) * t222 + rSges(4,2) * t221 - rSges(4,3) * t291;
t178 = Icges(4,1) * t222 + Icges(4,4) * t221 - Icges(4,5) * t291;
t177 = Icges(4,4) * t222 + Icges(4,2) * t221 - Icges(4,6) * t291;
t173 = rSges(5,1) * t213 - rSges(5,2) * t212 - rSges(5,3) * t291;
t172 = Icges(5,1) * t213 - Icges(5,4) * t212 - Icges(5,5) * t291;
t171 = Icges(5,4) * t213 - Icges(5,2) * t212 - Icges(5,6) * t291;
t170 = Icges(5,5) * t213 - Icges(5,6) * t212 - Icges(5,3) * t291;
t169 = t198 * t264 + t225 * t261;
t168 = -t198 * t261 + t225 * t264;
t167 = t196 * t264 + t223 * t261;
t166 = -t196 * t261 + t223 * t264;
t164 = qJD(5) * t197 + t204;
t163 = qJD(5) * t195 + t203;
t162 = pkin(4) * t198 + pkin(9) * t197;
t161 = pkin(4) * t196 + pkin(9) * t195;
t160 = rSges(4,1) * t202 + rSges(4,2) * t201 + rSges(4,3) * t225;
t159 = rSges(4,1) * t200 + rSges(4,2) * t199 + rSges(4,3) * t223;
t158 = Icges(4,1) * t202 + Icges(4,4) * t201 + Icges(4,5) * t225;
t157 = Icges(4,1) * t200 + Icges(4,4) * t199 + Icges(4,5) * t223;
t156 = Icges(4,4) * t202 + Icges(4,2) * t201 + Icges(4,6) * t225;
t155 = Icges(4,4) * t200 + Icges(4,2) * t199 + Icges(4,6) * t223;
t150 = rSges(5,1) * t198 - rSges(5,2) * t197 + rSges(5,3) * t225;
t149 = rSges(5,1) * t196 - rSges(5,2) * t195 + rSges(5,3) * t223;
t148 = Icges(5,1) * t198 - Icges(5,4) * t197 + Icges(5,5) * t225;
t147 = Icges(5,1) * t196 - Icges(5,4) * t195 + Icges(5,5) * t223;
t146 = Icges(5,4) * t198 - Icges(5,2) * t197 + Icges(5,6) * t225;
t145 = Icges(5,4) * t196 - Icges(5,2) * t195 + Icges(5,6) * t223;
t144 = Icges(5,5) * t198 - Icges(5,6) * t197 + Icges(5,3) * t225;
t143 = Icges(5,5) * t196 - Icges(5,6) * t195 + Icges(5,3) * t223;
t142 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t212;
t141 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t212;
t140 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t212;
t139 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t212;
t136 = -t188 * t237 + t214 * t235 + t275;
t135 = t189 * t237 - t214 * t236 + t272;
t134 = t188 * t236 - t189 * t235 + t274;
t133 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t197;
t132 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t195;
t131 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t197;
t130 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t195;
t129 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t197;
t128 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t195;
t127 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t197;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t195;
t125 = t179 * t235 + (-t159 - t191) * t237 + t273;
t124 = t160 * t237 + (-t179 - t227) * t236 + t271;
t123 = t159 * t236 + (-t160 - t192) * t235 + t270;
t122 = -t149 * t219 + t173 * t203 + t269;
t121 = t150 * t219 - t173 * t204 + t268;
t120 = t149 * t204 - t150 * t203 + t267;
t119 = -t132 * t190 + t142 * t163 - t161 * t219 + t180 * t203 + t269;
t118 = t133 * t190 - t142 * t164 + t162 * t219 - t180 * t204 + t268;
t117 = t132 * t164 - t133 * t163 + t161 * t204 - t162 * t203 + t267;
t1 = m(1) * (t231 ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + m(2) * (t205 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + t204 * ((t144 * t225 - t146 * t197 + t148 * t198) * t204 + (t143 * t225 - t145 * t197 + t147 * t198) * t203 + (t170 * t225 - t171 * t197 + t172 * t198) * t219) / 0.2e1 + t203 * ((t144 * t223 - t146 * t195 + t148 * t196) * t204 + (t143 * t223 - t145 * t195 + t147 * t196) * t203 + (t170 * t223 - t171 * t195 + t172 * t196) * t219) / 0.2e1 + t219 * ((-t144 * t291 - t146 * t212 + t148 * t213) * t204 + (-t143 * t291 - t145 * t212 + t147 * t213) * t203 + (-t170 * t291 - t171 * t212 + t172 * t213) * t219) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t164 * ((t127 * t197 + t129 * t168 + t131 * t169) * t164 + (t126 * t197 + t128 * t168 + t130 * t169) * t163 + (t139 * t197 + t140 * t168 + t141 * t169) * t190) / 0.2e1 + t163 * ((t127 * t195 + t129 * t166 + t131 * t167) * t164 + (t126 * t195 + t128 * t166 + t130 * t167) * t163 + (t139 * t195 + t140 * t166 + t141 * t167) * t190) / 0.2e1 + t190 * ((t127 * t212 + t129 * t193 + t131 * t194) * t164 + (t126 * t212 + t128 * t193 + t130 * t194) * t163 + (t139 * t212 + t140 * t193 + t141 * t194) * t190) / 0.2e1 + ((t177 * t199 + t178 * t200 - t209 * t290 + t224 * t211 + t223 * t301) * t237 + (t156 * t199 + t158 * t200 - t183 * t290 + t224 * t187 + t223 * t302) * t236 + (t155 * t199 + t157 * t200 - t182 * t290 + t224 * t186 + t223 * t303) * t235) * t235 / 0.2e1 + ((t177 * t201 + t178 * t202 + t209 * t292 + t211 * t226 + t225 * t301) * t237 + (t156 * t201 + t158 * t202 + t183 * t292 + t187 * t226 + t225 * t302) * t236 + (t155 * t201 + t157 * t202 + t182 * t292 + t186 * t226 + t225 * t303) * t235) * t236 / 0.2e1 + ((t182 * t235 + t183 * t236 + t209 * t237) * t259 + ((t185 * t265 + t187 * t262) * t236 + (t184 * t265 + t186 * t262) * t235 + (t210 * t265 + t211 * t262) * t237) * t257 + (-t154 * t291 + t156 * t221 + t158 * t222) * t236 + (-t153 * t291 + t155 * t221 + t157 * t222) * t235 + (-t176 * t291 + t177 * t221 + t178 * t222) * t237) * t237 / 0.2e1 + ((-t263 * t240 + t242 * t266 + Icges(1,4)) * V_base(5) + (-t263 * t241 + t243 * t266 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t240 * t266 + t263 * t242 + Icges(1,2)) * V_base(5) + (t241 * t266 + t263 * t243 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t263 + Icges(2,6) * t266) * V_base(5) + (Icges(2,5) * t266 - Icges(2,6) * t263) * V_base(4) + Icges(2,3) * t253 / 0.2e1) * t253;
T = t1;
