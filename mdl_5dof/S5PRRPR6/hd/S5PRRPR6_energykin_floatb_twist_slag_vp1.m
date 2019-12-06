% Calculate kinetic energy for
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:00
% EndTime: 2019-12-05 16:30:04
% DurationCPUTime: 3.90s
% Computational Cost: add. (2112->355), mult. (4560->507), div. (0->0), fcn. (5454->12), ass. (0->152)
t296 = Icges(4,2) + Icges(5,3);
t249 = sin(pkin(9));
t252 = cos(pkin(9));
t257 = cos(qJ(2));
t253 = cos(pkin(5));
t256 = sin(qJ(2));
t277 = t253 * t256;
t211 = t249 * t257 + t252 * t277;
t250 = sin(pkin(5));
t255 = sin(qJ(3));
t279 = t250 * t255;
t287 = cos(qJ(3));
t197 = t211 * t287 - t252 * t279;
t276 = t253 * t257;
t210 = t249 * t256 - t252 * t276;
t248 = sin(pkin(10));
t251 = cos(pkin(10));
t165 = -t197 * t248 + t210 * t251;
t283 = t210 * t248;
t166 = t197 * t251 + t283;
t267 = t250 * t287;
t196 = t211 * t255 + t252 * t267;
t293 = -Icges(4,4) * t197 + Icges(5,5) * t166 - Icges(4,6) * t210 + Icges(5,6) * t165 + t296 * t196;
t213 = -t249 * t277 + t252 * t257;
t199 = t213 * t287 + t249 * t279;
t212 = t249 * t276 + t252 * t256;
t167 = -t199 * t248 + t212 * t251;
t282 = t212 * t248;
t168 = t199 * t251 + t282;
t198 = t213 * t255 - t249 * t267;
t292 = -Icges(4,4) * t199 + Icges(5,5) * t168 - Icges(4,6) * t212 + Icges(5,6) * t167 + t296 * t198;
t218 = t253 * t255 + t256 * t267;
t278 = t250 * t257;
t194 = -t218 * t248 - t251 * t278;
t268 = t248 * t278;
t195 = t218 * t251 - t268;
t217 = -t253 * t287 + t256 * t279;
t291 = -Icges(4,4) * t218 + Icges(5,5) * t195 + Icges(4,6) * t278 + Icges(5,6) * t194 + t296 * t217;
t286 = pkin(6) * t253;
t285 = pkin(4) * t251;
t284 = Icges(2,4) * t249;
t281 = t249 * t250;
t280 = t250 * t252;
t274 = qJD(2) * t250;
t273 = V_base(5) * qJ(1) + V_base(1);
t269 = qJD(1) + V_base(3);
t226 = t249 * t274 + V_base(4);
t237 = qJD(2) * t253 + V_base(6);
t193 = qJD(3) * t212 + t226;
t225 = -t252 * t274 + V_base(5);
t192 = qJD(3) * t210 + t225;
t214 = -qJD(3) * t278 + t237;
t220 = pkin(1) * t249 - pkin(6) * t280;
t266 = -t220 * V_base(6) + V_base(5) * t286 + t273;
t221 = pkin(1) * t252 + pkin(6) * t281;
t265 = V_base(4) * t220 - t221 * V_base(5) + t269;
t264 = V_base(6) * t221 + V_base(2) + (-qJ(1) - t286) * V_base(4);
t183 = pkin(2) * t211 + pkin(7) * t210;
t219 = (pkin(2) * t256 - pkin(7) * t257) * t250;
t263 = -t183 * t237 + t225 * t219 + t266;
t184 = pkin(2) * t213 + pkin(7) * t212;
t262 = t226 * t183 - t184 * t225 + t265;
t261 = t237 * t184 - t219 * t226 + t264;
t185 = t218 * pkin(3) + t217 * qJ(4);
t260 = qJD(4) * t198 + t192 * t185 + t263;
t157 = pkin(3) * t197 + qJ(4) * t196;
t259 = qJD(4) * t217 + t193 * t157 + t262;
t158 = pkin(3) * t199 + qJ(4) * t198;
t258 = qJD(4) * t196 + t214 * t158 + t261;
t247 = pkin(10) + qJ(5);
t245 = Icges(2,4) * t252;
t244 = cos(t247);
t243 = sin(t247);
t234 = rSges(2,1) * t252 - rSges(2,2) * t249;
t233 = rSges(2,1) * t249 + rSges(2,2) * t252;
t232 = Icges(2,1) * t252 - t284;
t231 = Icges(2,1) * t249 + t245;
t230 = -Icges(2,2) * t249 + t245;
t229 = Icges(2,2) * t252 + t284;
t224 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t223 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t222 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t206 = t253 * rSges(3,3) + (rSges(3,1) * t256 + rSges(3,2) * t257) * t250;
t205 = Icges(3,5) * t253 + (Icges(3,1) * t256 + Icges(3,4) * t257) * t250;
t204 = Icges(3,6) * t253 + (Icges(3,4) * t256 + Icges(3,2) * t257) * t250;
t203 = Icges(3,3) * t253 + (Icges(3,5) * t256 + Icges(3,6) * t257) * t250;
t202 = V_base(5) * rSges(2,3) - t233 * V_base(6) + t273;
t201 = t234 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t191 = t233 * V_base(4) - t234 * V_base(5) + t269;
t188 = t218 * t244 - t243 * t278;
t187 = -t218 * t243 - t244 * t278;
t186 = qJD(5) * t217 + t214;
t182 = t218 * rSges(4,1) - t217 * rSges(4,2) - rSges(4,3) * t278;
t181 = Icges(4,1) * t218 - Icges(4,4) * t217 - Icges(4,5) * t278;
t179 = Icges(4,5) * t218 - Icges(4,6) * t217 - Icges(4,3) * t278;
t178 = rSges(3,1) * t213 - rSges(3,2) * t212 + rSges(3,3) * t281;
t177 = rSges(3,1) * t211 - rSges(3,2) * t210 - rSges(3,3) * t280;
t176 = Icges(3,1) * t213 - Icges(3,4) * t212 + Icges(3,5) * t281;
t175 = Icges(3,1) * t211 - Icges(3,4) * t210 - Icges(3,5) * t280;
t174 = Icges(3,4) * t213 - Icges(3,2) * t212 + Icges(3,6) * t281;
t173 = Icges(3,4) * t211 - Icges(3,2) * t210 - Icges(3,6) * t280;
t172 = Icges(3,5) * t213 - Icges(3,6) * t212 + Icges(3,3) * t281;
t171 = Icges(3,5) * t211 - Icges(3,6) * t210 - Icges(3,3) * t280;
t164 = t199 * t244 + t212 * t243;
t163 = -t199 * t243 + t212 * t244;
t162 = t197 * t244 + t210 * t243;
t161 = -t197 * t243 + t210 * t244;
t160 = qJD(5) * t198 + t193;
t159 = qJD(5) * t196 + t192;
t154 = rSges(5,1) * t195 + rSges(5,2) * t194 + rSges(5,3) * t217;
t153 = rSges(4,1) * t199 - rSges(4,2) * t198 + rSges(4,3) * t212;
t152 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t210;
t151 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t217;
t150 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t217;
t148 = Icges(4,1) * t199 - Icges(4,4) * t198 + Icges(4,5) * t212;
t147 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t210;
t144 = Icges(4,5) * t199 - Icges(4,6) * t198 + Icges(4,3) * t212;
t143 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t210;
t142 = -pkin(4) * t268 + pkin(8) * t217 + t218 * t285;
t141 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t217;
t140 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t217;
t139 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t217;
t138 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t217;
t136 = -t177 * t237 + t206 * t225 + t266;
t135 = t178 * t237 - t206 * t226 + t264;
t134 = rSges(5,1) * t168 + rSges(5,2) * t167 + rSges(5,3) * t198;
t133 = rSges(5,1) * t166 + rSges(5,2) * t165 + rSges(5,3) * t196;
t132 = Icges(5,1) * t168 + Icges(5,4) * t167 + Icges(5,5) * t198;
t131 = Icges(5,1) * t166 + Icges(5,4) * t165 + Icges(5,5) * t196;
t130 = Icges(5,4) * t168 + Icges(5,2) * t167 + Icges(5,6) * t198;
t129 = Icges(5,4) * t166 + Icges(5,2) * t165 + Icges(5,6) * t196;
t126 = t177 * t226 - t178 * t225 + t265;
t125 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t198;
t124 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t196;
t123 = Icges(6,1) * t164 + Icges(6,4) * t163 + Icges(6,5) * t198;
t122 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t196;
t121 = Icges(6,4) * t164 + Icges(6,2) * t163 + Icges(6,6) * t198;
t120 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t196;
t119 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t198;
t118 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t196;
t117 = pkin(4) * t282 + pkin(8) * t198 + t199 * t285;
t116 = pkin(4) * t283 + pkin(8) * t196 + t197 * t285;
t115 = -t152 * t214 + t182 * t192 + t263;
t114 = t153 * t214 - t182 * t193 + t261;
t113 = t152 * t193 - t153 * t192 + t262;
t112 = t154 * t192 + (-t133 - t157) * t214 + t260;
t111 = t134 * t214 + (-t154 - t185) * t193 + t258;
t110 = t133 * t193 + (-t134 - t158) * t192 + t259;
t109 = -t124 * t186 + t141 * t159 + t142 * t192 + (-t116 - t157) * t214 + t260;
t108 = t117 * t214 + t125 * t186 - t141 * t160 + (-t142 - t185) * t193 + t258;
t107 = t116 * t193 + t124 * t160 - t125 * t159 + (-t117 - t158) * t192 + t259;
t1 = m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t191 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(3) * (t126 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t226 * ((t172 * t281 - t174 * t212 + t176 * t213) * t226 + (t171 * t281 - t173 * t212 + t175 * t213) * t225 + (t203 * t281 - t204 * t212 + t205 * t213) * t237) / 0.2e1 + t225 * ((-t172 * t280 - t174 * t210 + t176 * t211) * t226 + (-t171 * t280 - t173 * t210 + t175 * t211) * t225 + (-t203 * t280 - t204 * t210 + t205 * t211) * t237) / 0.2e1 + t237 * ((t171 * t225 + t172 * t226 + t203 * t237) * t253 + ((t174 * t257 + t176 * t256) * t226 + (t173 * t257 + t175 * t256) * t225 + (t204 * t257 + t205 * t256) * t237) * t250) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t160 * ((t198 * t119 + t163 * t121 + t164 * t123) * t160 + (t118 * t198 + t120 * t163 + t122 * t164) * t159 + (t138 * t198 + t139 * t163 + t140 * t164) * t186) / 0.2e1 + t159 * ((t119 * t196 + t121 * t161 + t123 * t162) * t160 + (t196 * t118 + t161 * t120 + t162 * t122) * t159 + (t138 * t196 + t139 * t161 + t140 * t162) * t186) / 0.2e1 + t186 * ((t119 * t217 + t121 * t187 + t123 * t188) * t160 + (t118 * t217 + t120 * t187 + t122 * t188) * t159 + (t138 * t217 + t139 * t187 + t140 * t188) * t186) / 0.2e1 + ((t150 * t165 + t151 * t166 + t179 * t210 + t181 * t197 + t196 * t291) * t214 + (t130 * t165 + t132 * t166 + t144 * t210 + t148 * t197 + t196 * t292) * t193 + (t129 * t165 + t131 * t166 + t143 * t210 + t147 * t197 + t293 * t196) * t192) * t192 / 0.2e1 + ((t150 * t167 + t151 * t168 + t179 * t212 + t181 * t199 + t198 * t291) * t214 + (t130 * t167 + t132 * t168 + t144 * t212 + t148 * t199 + t292 * t198) * t193 + (t129 * t167 + t131 * t168 + t143 * t212 + t147 * t199 + t198 * t293) * t192) * t193 / 0.2e1 + ((t150 * t194 + t151 * t195 - t179 * t278 + t218 * t181 + t291 * t217) * t214 + (t130 * t194 + t132 * t195 - t144 * t278 + t218 * t148 + t217 * t292) * t193 + (t129 * t194 + t131 * t195 - t143 * t278 + t218 * t147 + t217 * t293) * t192) * t214 / 0.2e1 + ((-t229 * t249 + t231 * t252 + Icges(1,4)) * V_base(5) + (-t230 * t249 + t232 * t252 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t229 * t252 + t231 * t249 + Icges(1,2)) * V_base(5) + (t230 * t252 + t232 * t249 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t249 + Icges(2,6) * t252 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t252 - Icges(2,6) * t249 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
