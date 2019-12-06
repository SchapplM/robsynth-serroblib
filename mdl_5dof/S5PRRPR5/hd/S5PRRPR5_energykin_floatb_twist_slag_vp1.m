% Calculate kinetic energy for
% S5PRRPR5
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:13
% EndTime: 2019-12-05 16:25:16
% DurationCPUTime: 3.72s
% Computational Cost: add. (2234->355), mult. (4024->510), div. (0->0), fcn. (4698->12), ass. (0->154)
t307 = Icges(4,3) + Icges(5,3);
t255 = sin(pkin(9));
t257 = cos(pkin(9));
t265 = cos(qJ(2));
t258 = cos(pkin(5));
t262 = sin(qJ(2));
t287 = t258 * t262;
t220 = t255 * t265 + t257 * t287;
t284 = qJ(3) + pkin(10);
t252 = sin(t284);
t256 = sin(pkin(5));
t276 = cos(t284);
t275 = t256 * t276;
t193 = t220 * t252 + t257 * t275;
t293 = t256 * t257;
t194 = t220 * t276 - t252 * t293;
t261 = sin(qJ(3));
t264 = cos(qJ(3));
t290 = t256 * t264;
t202 = -t220 * t261 - t257 * t290;
t292 = t256 * t261;
t277 = t257 * t292;
t203 = t220 * t264 - t277;
t286 = t258 * t265;
t219 = t255 * t262 - t257 * t286;
t304 = Icges(4,5) * t203 + Icges(5,5) * t194 + Icges(4,6) * t202 - Icges(5,6) * t193 + t307 * t219;
t222 = -t255 * t287 + t257 * t265;
t195 = t222 * t252 - t255 * t275;
t294 = t255 * t256;
t196 = t222 * t276 + t252 * t294;
t204 = -t222 * t261 + t255 * t290;
t278 = t255 * t292;
t205 = t222 * t264 + t278;
t221 = t255 * t286 + t257 * t262;
t303 = Icges(4,5) * t205 + Icges(5,5) * t196 + Icges(4,6) * t204 - Icges(5,6) * t195 + t307 * t221;
t291 = t256 * t262;
t212 = t252 * t291 - t258 * t276;
t213 = t258 * t252 + t262 * t275;
t226 = t258 * t264 - t261 * t291;
t288 = t258 * t261;
t227 = t262 * t290 + t288;
t289 = t256 * t265;
t302 = Icges(4,5) * t227 + Icges(5,5) * t213 + Icges(4,6) * t226 - Icges(5,6) * t212 - t307 * t289;
t298 = pkin(6) * t258;
t297 = pkin(3) * t264;
t295 = Icges(2,4) * t255;
t285 = qJD(2) * t256;
t283 = V_base(5) * qJ(1) + V_base(1);
t279 = qJD(1) + V_base(3);
t236 = t255 * t285 + V_base(4);
t247 = qJD(2) * t258 + V_base(6);
t201 = qJD(3) * t221 + t236;
t235 = -t257 * t285 + V_base(5);
t200 = qJD(3) * t219 + t235;
t223 = -qJD(3) * t289 + t247;
t229 = pkin(1) * t255 - pkin(6) * t293;
t274 = -t229 * V_base(6) + V_base(5) * t298 + t283;
t230 = pkin(1) * t257 + pkin(6) * t294;
t273 = V_base(4) * t229 - V_base(5) * t230 + t279;
t272 = V_base(6) * t230 + V_base(2) + (-qJ(1) - t298) * V_base(4);
t191 = pkin(2) * t220 + pkin(7) * t219;
t228 = (pkin(2) * t262 - pkin(7) * t265) * t256;
t271 = -t191 * t247 + t235 * t228 + t274;
t192 = pkin(2) * t222 + pkin(7) * t221;
t270 = t236 * t191 - t235 * t192 + t273;
t269 = t247 * t192 - t228 * t236 + t272;
t188 = pkin(3) * t288 + (-qJ(4) * t265 + t262 * t297) * t256;
t268 = qJD(4) * t221 + t200 * t188 + t271;
t152 = pkin(3) * t278 + qJ(4) * t221 + t222 * t297;
t267 = qJD(4) * t219 + t223 * t152 + t269;
t151 = -pkin(3) * t277 + qJ(4) * t219 + t220 * t297;
t266 = -qJD(4) * t289 + t201 * t151 + t270;
t263 = cos(qJ(5));
t260 = sin(qJ(5));
t253 = Icges(2,4) * t257;
t244 = rSges(2,1) * t257 - rSges(2,2) * t255;
t243 = rSges(2,1) * t255 + rSges(2,2) * t257;
t242 = Icges(2,1) * t257 - t295;
t241 = Icges(2,1) * t255 + t253;
t240 = -Icges(2,2) * t255 + t253;
t239 = Icges(2,2) * t257 + t295;
t234 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t233 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t232 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = t258 * rSges(3,3) + (rSges(3,1) * t262 + rSges(3,2) * t265) * t256;
t211 = Icges(3,5) * t258 + (Icges(3,1) * t262 + Icges(3,4) * t265) * t256;
t210 = Icges(3,6) * t258 + (Icges(3,4) * t262 + Icges(3,2) * t265) * t256;
t209 = Icges(3,3) * t258 + (Icges(3,5) * t262 + Icges(3,6) * t265) * t256;
t208 = V_base(5) * rSges(2,3) - t243 * V_base(6) + t283;
t207 = t244 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t199 = t243 * V_base(4) - t244 * V_base(5) + t279;
t198 = t213 * t263 - t260 * t289;
t197 = -t213 * t260 - t263 * t289;
t190 = qJD(5) * t212 + t223;
t189 = t227 * rSges(4,1) + t226 * rSges(4,2) - rSges(4,3) * t289;
t187 = Icges(4,1) * t227 + Icges(4,4) * t226 - Icges(4,5) * t289;
t186 = Icges(4,4) * t227 + Icges(4,2) * t226 - Icges(4,6) * t289;
t184 = rSges(3,1) * t222 - rSges(3,2) * t221 + rSges(3,3) * t294;
t183 = rSges(3,1) * t220 - rSges(3,2) * t219 - rSges(3,3) * t293;
t182 = pkin(4) * t213 + pkin(8) * t212;
t181 = Icges(3,1) * t222 - Icges(3,4) * t221 + Icges(3,5) * t294;
t180 = Icges(3,1) * t220 - Icges(3,4) * t219 - Icges(3,5) * t293;
t179 = Icges(3,4) * t222 - Icges(3,2) * t221 + Icges(3,6) * t294;
t178 = Icges(3,4) * t220 - Icges(3,2) * t219 - Icges(3,6) * t293;
t177 = Icges(3,5) * t222 - Icges(3,6) * t221 + Icges(3,3) * t294;
t176 = Icges(3,5) * t220 - Icges(3,6) * t219 - Icges(3,3) * t293;
t173 = t213 * rSges(5,1) - t212 * rSges(5,2) - rSges(5,3) * t289;
t172 = Icges(5,1) * t213 - Icges(5,4) * t212 - Icges(5,5) * t289;
t171 = Icges(5,4) * t213 - Icges(5,2) * t212 - Icges(5,6) * t289;
t169 = t196 * t263 + t221 * t260;
t168 = -t196 * t260 + t221 * t263;
t167 = t194 * t263 + t219 * t260;
t166 = -t194 * t260 + t219 * t263;
t165 = qJD(5) * t195 + t201;
t164 = qJD(5) * t193 + t200;
t163 = pkin(4) * t196 + pkin(8) * t195;
t162 = pkin(4) * t194 + pkin(8) * t193;
t160 = rSges(4,1) * t205 + rSges(4,2) * t204 + rSges(4,3) * t221;
t159 = rSges(4,1) * t203 + rSges(4,2) * t202 + rSges(4,3) * t219;
t158 = Icges(4,1) * t205 + Icges(4,4) * t204 + Icges(4,5) * t221;
t157 = Icges(4,1) * t203 + Icges(4,4) * t202 + Icges(4,5) * t219;
t156 = Icges(4,4) * t205 + Icges(4,2) * t204 + Icges(4,6) * t221;
t155 = Icges(4,4) * t203 + Icges(4,2) * t202 + Icges(4,6) * t219;
t150 = rSges(5,1) * t196 - rSges(5,2) * t195 + rSges(5,3) * t221;
t149 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t219;
t148 = Icges(5,1) * t196 - Icges(5,4) * t195 + Icges(5,5) * t221;
t147 = Icges(5,1) * t194 - Icges(5,4) * t193 + Icges(5,5) * t219;
t146 = Icges(5,4) * t196 - Icges(5,2) * t195 + Icges(5,6) * t221;
t145 = Icges(5,4) * t194 - Icges(5,2) * t193 + Icges(5,6) * t219;
t142 = rSges(6,1) * t198 + rSges(6,2) * t197 + rSges(6,3) * t212;
t141 = Icges(6,1) * t198 + Icges(6,4) * t197 + Icges(6,5) * t212;
t140 = Icges(6,4) * t198 + Icges(6,2) * t197 + Icges(6,6) * t212;
t139 = Icges(6,5) * t198 + Icges(6,6) * t197 + Icges(6,3) * t212;
t137 = -t183 * t247 + t214 * t235 + t274;
t136 = t184 * t247 - t214 * t236 + t272;
t134 = t183 * t236 - t184 * t235 + t273;
t133 = rSges(6,1) * t169 + rSges(6,2) * t168 + rSges(6,3) * t195;
t132 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t193;
t131 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t195;
t130 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t193;
t129 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t195;
t128 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t193;
t127 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t195;
t126 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t193;
t125 = -t159 * t223 + t189 * t200 + t271;
t124 = t160 * t223 - t189 * t201 + t269;
t123 = t159 * t201 - t160 * t200 + t270;
t122 = t173 * t200 + (-t149 - t151) * t223 + t268;
t121 = t150 * t223 + (-t173 - t188) * t201 + t267;
t120 = t201 * t149 + (-t150 - t152) * t200 + t266;
t119 = -t132 * t190 + t142 * t164 + t182 * t200 + (-t151 - t162) * t223 + t268;
t118 = t133 * t190 - t142 * t165 + t163 * t223 + (-t182 - t188) * t201 + t267;
t117 = t165 * t132 - t164 * t133 + t201 * t162 + (-t152 - t163) * t200 + t266;
t1 = m(1) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + m(2) * (t199 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t236 * ((t177 * t294 - t179 * t221 + t181 * t222) * t236 + (t176 * t294 - t178 * t221 + t180 * t222) * t235 + (t209 * t294 - t210 * t221 + t211 * t222) * t247) / 0.2e1 + t235 * ((-t177 * t293 - t179 * t219 + t181 * t220) * t236 + (-t176 * t293 - t178 * t219 + t180 * t220) * t235 + (-t209 * t293 - t210 * t219 + t211 * t220) * t247) / 0.2e1 + t247 * ((t176 * t235 + t177 * t236 + t209 * t247) * t258 + ((t179 * t265 + t181 * t262) * t236 + (t178 * t265 + t180 * t262) * t235 + (t210 * t265 + t211 * t262) * t247) * t256) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t165 * ((t127 * t195 + t129 * t168 + t131 * t169) * t165 + (t126 * t195 + t128 * t168 + t130 * t169) * t164 + (t139 * t195 + t140 * t168 + t141 * t169) * t190) / 0.2e1 + t164 * ((t127 * t193 + t129 * t166 + t131 * t167) * t165 + (t126 * t193 + t128 * t166 + t130 * t167) * t164 + (t139 * t193 + t140 * t166 + t141 * t167) * t190) / 0.2e1 + t190 * ((t127 * t212 + t129 * t197 + t131 * t198) * t165 + (t126 * t212 + t128 * t197 + t130 * t198) * t164 + (t139 * t212 + t140 * t197 + t141 * t198) * t190) / 0.2e1 + ((-t171 * t193 + t172 * t194 + t186 * t202 + t187 * t203 + t219 * t302) * t223 + (-t146 * t193 + t148 * t194 + t156 * t202 + t158 * t203 + t219 * t303) * t201 + (-t145 * t193 + t147 * t194 + t155 * t202 + t157 * t203 + t219 * t304) * t200) * t200 / 0.2e1 + ((-t171 * t195 + t172 * t196 + t186 * t204 + t187 * t205 + t221 * t302) * t223 + (-t146 * t195 + t148 * t196 + t156 * t204 + t158 * t205 + t221 * t303) * t201 + (-t145 * t195 + t147 * t196 + t155 * t204 + t157 * t205 + t221 * t304) * t200) * t201 / 0.2e1 + ((-t212 * t171 + t213 * t172 + t226 * t186 + t227 * t187 - t289 * t302) * t223 + (-t212 * t146 + t213 * t148 + t226 * t156 + t227 * t158 - t289 * t303) * t201 + (-t212 * t145 + t213 * t147 + t226 * t155 + t227 * t157 - t289 * t304) * t200) * t223 / 0.2e1 + ((-t239 * t255 + t241 * t257 + Icges(1,4)) * V_base(5) + (-t240 * t255 + t242 * t257 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t239 * t257 + t241 * t255 + Icges(1,2)) * V_base(5) + (t240 * t257 + t242 * t255 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t255 + Icges(2,6) * t257 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t257 - Icges(2,6) * t255 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
