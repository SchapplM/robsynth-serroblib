% Calculate kinetic energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:38
% EndTime: 2019-12-05 16:08:41
% DurationCPUTime: 2.52s
% Computational Cost: add. (1255->271), mult. (1916->378), div. (0->0), fcn. (1876->8), ass. (0->129)
t276 = Icges(5,1) + Icges(6,1);
t275 = -Icges(5,4) + Icges(6,5);
t274 = Icges(6,4) + Icges(5,5);
t273 = Icges(5,2) + Icges(6,3);
t272 = -Icges(6,6) + Icges(5,6);
t271 = -Icges(5,3) - Icges(6,2) - Icges(4,3);
t270 = rSges(6,1) + pkin(4);
t269 = rSges(6,3) + qJ(5);
t208 = qJ(3) + pkin(8);
t204 = sin(t208);
t205 = cos(t208);
t210 = cos(pkin(7));
t209 = sin(pkin(7));
t215 = cos(qJ(2));
t247 = t209 * t215;
t155 = t204 * t247 + t210 * t205;
t156 = -t210 * t204 + t205 * t247;
t213 = sin(qJ(2));
t248 = t209 * t213;
t266 = t273 * t155 + t275 * t156 - t272 * t248;
t244 = t210 * t215;
t157 = t204 * t244 - t209 * t205;
t158 = t209 * t204 + t205 * t244;
t245 = t210 * t213;
t265 = t273 * t157 + t275 * t158 - t272 * t245;
t264 = t275 * t155 + t276 * t156 + t274 * t248;
t263 = t275 * t157 + t276 * t158 + t274 * t245;
t262 = t272 * t215 + (t273 * t204 + t275 * t205) * t213;
t261 = -t274 * t215 + (t275 * t204 + t276 * t205) * t213;
t214 = cos(qJ(3));
t212 = sin(qJ(3));
t243 = t212 * t215;
t170 = -t209 * t243 - t210 * t214;
t242 = t214 * t215;
t246 = t210 * t212;
t171 = t209 * t242 - t246;
t260 = Icges(4,5) * t171 + Icges(4,6) * t170 - t272 * t155 + t274 * t156 - t271 * t248;
t172 = t209 * t214 - t210 * t243;
t249 = t209 * t212;
t173 = t210 * t242 + t249;
t259 = Icges(4,5) * t173 + Icges(4,6) * t172 - t272 * t157 + t274 * t158 - t271 * t245;
t258 = t271 * t215 + (Icges(4,5) * t214 - Icges(4,6) * t212 - t272 * t204 + t274 * t205) * t213;
t254 = pkin(3) * t214;
t252 = Icges(2,4) * t209;
t251 = Icges(3,4) * t213;
t250 = Icges(3,4) * t215;
t241 = rSges(6,2) * t248 + t269 * t155 + t270 * t156;
t240 = rSges(6,2) * t245 + t269 * t157 + t270 * t158;
t239 = -rSges(6,2) * t215 + (t269 * t204 + t270 * t205) * t213;
t238 = qJD(3) * t213;
t237 = qJD(4) * t213;
t236 = V_base(5) * qJ(1) + V_base(1);
t232 = qJD(1) + V_base(3);
t198 = qJD(2) * t209 + V_base(4);
t231 = pkin(2) * t215 + pkin(6) * t213;
t197 = -qJD(2) * t210 + V_base(5);
t230 = rSges(3,1) * t215 - rSges(3,2) * t213;
t229 = Icges(3,1) * t215 - t251;
t228 = -Icges(3,2) * t213 + t250;
t227 = Icges(3,5) * t215 - Icges(3,6) * t213;
t190 = pkin(1) * t210 + pkin(5) * t209;
t226 = -V_base(4) * qJ(1) + V_base(6) * t190 + V_base(2);
t189 = pkin(1) * t209 - pkin(5) * t210;
t225 = V_base(4) * t189 - V_base(5) * t190 + t232;
t224 = qJ(4) * t213 + t215 * t254;
t174 = t231 * t209;
t195 = t213 * pkin(2) - pkin(6) * t215;
t223 = t197 * t195 + (-t174 - t189) * V_base(6) + t236;
t222 = (-Icges(3,3) * t210 + t209 * t227) * t197 + (Icges(3,3) * t209 + t210 * t227) * t198 + (Icges(3,5) * t213 + Icges(3,6) * t215) * V_base(6);
t175 = t231 * t210;
t221 = V_base(6) * t175 - t195 * t198 + t226;
t139 = -qJ(4) * t215 + t213 * t254;
t168 = t209 * t238 + t197;
t220 = t168 * t139 + t210 * t237 + t223;
t219 = t198 * t174 - t197 * t175 + t225;
t132 = pkin(3) * t249 + t210 * t224;
t199 = -qJD(3) * t215 + V_base(6);
t218 = t199 * t132 + t209 * t237 + t221;
t131 = -pkin(3) * t246 + t209 * t224;
t169 = t210 * t238 + t198;
t217 = -qJD(4) * t215 + t169 * t131 + t219;
t151 = -Icges(3,6) * t210 + t209 * t228;
t152 = Icges(3,6) * t209 + t210 * t228;
t153 = -Icges(3,5) * t210 + t209 * t229;
t154 = Icges(3,5) * t209 + t210 * t229;
t192 = Icges(3,2) * t215 + t251;
t193 = Icges(3,1) * t213 + t250;
t216 = (-t152 * t213 + t154 * t215) * t198 + (-t151 * t213 + t153 * t215) * t197 + (-t192 * t213 + t193 * t215) * V_base(6);
t206 = Icges(2,4) * t210;
t194 = t213 * rSges(3,1) + rSges(3,2) * t215;
t188 = rSges(2,1) * t210 - rSges(2,2) * t209;
t187 = rSges(2,1) * t209 + rSges(2,2) * t210;
t186 = Icges(2,1) * t210 - t252;
t185 = Icges(2,1) * t209 + t206;
t184 = -Icges(2,2) * t209 + t206;
t183 = Icges(2,2) * t210 + t252;
t180 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t179 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t178 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = -rSges(4,3) * t215 + (rSges(4,1) * t214 - rSges(4,2) * t212) * t213;
t163 = -Icges(4,5) * t215 + (Icges(4,1) * t214 - Icges(4,4) * t212) * t213;
t162 = -Icges(4,6) * t215 + (Icges(4,4) * t214 - Icges(4,2) * t212) * t213;
t160 = t209 * rSges(3,3) + t210 * t230;
t159 = -t210 * rSges(3,3) + t209 * t230;
t148 = -rSges(5,3) * t215 + (rSges(5,1) * t205 - rSges(5,2) * t204) * t213;
t138 = V_base(5) * rSges(2,3) - t187 * V_base(6) + t236;
t137 = t188 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t136 = t187 * V_base(4) - t188 * V_base(5) + t232;
t134 = rSges(4,1) * t173 + rSges(4,2) * t172 + rSges(4,3) * t245;
t133 = rSges(4,1) * t171 + rSges(4,2) * t170 + rSges(4,3) * t248;
t130 = Icges(4,1) * t173 + Icges(4,4) * t172 + Icges(4,5) * t245;
t129 = Icges(4,1) * t171 + Icges(4,4) * t170 + Icges(4,5) * t248;
t128 = Icges(4,4) * t173 + Icges(4,2) * t172 + Icges(4,6) * t245;
t127 = Icges(4,4) * t171 + Icges(4,2) * t170 + Icges(4,6) * t248;
t121 = rSges(5,1) * t158 - rSges(5,2) * t157 + rSges(5,3) * t245;
t119 = rSges(5,1) * t156 - rSges(5,2) * t155 + rSges(5,3) * t248;
t104 = t194 * t197 + (-t159 - t189) * V_base(6) + t236;
t103 = t160 * V_base(6) - t194 * t198 + t226;
t102 = t159 * t198 - t160 * t197 + t225;
t101 = -t133 * t199 + t165 * t168 + t223;
t100 = t134 * t199 - t165 * t169 + t221;
t99 = t133 * t169 - t134 * t168 + t219;
t98 = t148 * t168 + (-t119 - t131) * t199 + t220;
t97 = t121 * t199 + (-t139 - t148) * t169 + t218;
t96 = t169 * t119 + (-t121 - t132) * t168 + t217;
t95 = qJD(5) * t157 + t239 * t168 + (-t131 - t241) * t199 + t220;
t94 = qJD(5) * t155 + t240 * t199 + (-t139 - t239) * t169 + t218;
t93 = qJD(5) * t213 * t204 + t241 * t169 + (-t132 - t240) * t168 + t217;
t1 = m(1) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(2) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t198 * (t209 * t222 + t210 * t216) / 0.2e1 + t197 * (t209 * t216 - t210 * t222) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + ((-t183 * t209 + t185 * t210 + Icges(1,4)) * V_base(5) + (-t184 * t209 + t186 * t210 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t183 * t210 + t185 * t209 + Icges(1,2)) * V_base(5) + (t184 * t210 + t186 * t209 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t155 * t262 + t156 * t261 + t162 * t170 + t163 * t171 + t248 * t258) * t199 + (t128 * t170 + t130 * t171 + t155 * t265 + t156 * t263 + t248 * t259) * t169 + (t127 * t170 + t129 * t171 + t266 * t155 + t264 * t156 + t260 * t248) * t168) * t168 / 0.2e1 + ((t157 * t262 + t158 * t261 + t162 * t172 + t163 * t173 + t245 * t258) * t199 + (t128 * t172 + t130 * t173 + t265 * t157 + t263 * t158 + t259 * t245) * t169 + (t127 * t172 + t129 * t173 + t157 * t266 + t158 * t264 + t245 * t260) * t168) * t169 / 0.2e1 + ((-t168 * t260 - t169 * t259 - t258 * t199) * t215 + ((-t162 * t212 + t163 * t214 + t204 * t262 + t205 * t261) * t199 + (-t128 * t212 + t130 * t214 + t204 * t265 + t205 * t263) * t169 + (-t127 * t212 + t129 * t214 + t204 * t266 + t264 * t205) * t168) * t213) * t199 / 0.2e1 + ((t152 * t215 + t213 * t154) * t198 + (t151 * t215 + t213 * t153) * t197 + (t192 * t215 + t213 * t193 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t210 - Icges(2,6) * t209 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t209 + Icges(2,6) * t210 + Icges(1,6));
T = t1;
