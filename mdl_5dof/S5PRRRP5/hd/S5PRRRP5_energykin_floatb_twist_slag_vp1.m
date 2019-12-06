% Calculate kinetic energy for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:40
% EndTime: 2019-12-05 16:47:43
% DurationCPUTime: 2.87s
% Computational Cost: add. (1326->283), mult. (2009->419), div. (0->0), fcn. (1959->8), ass. (0->137)
t276 = Icges(5,1) + Icges(6,1);
t275 = Icges(5,4) + Icges(6,4);
t274 = -Icges(6,5) - Icges(5,5);
t273 = Icges(5,2) + Icges(6,2);
t272 = -Icges(6,6) - Icges(5,6);
t271 = -Icges(6,3) - Icges(5,3);
t204 = qJ(3) + qJ(4);
t199 = sin(t204);
t200 = cos(t204);
t206 = cos(pkin(8));
t205 = sin(pkin(8));
t210 = cos(qJ(2));
t246 = t205 * t210;
t150 = -t199 * t246 - t200 * t206;
t151 = -t199 * t206 + t200 * t246;
t208 = sin(qJ(2));
t247 = t205 * t208;
t270 = -t272 * t150 - t274 * t151 - t271 * t247;
t243 = t206 * t210;
t152 = -t199 * t243 + t200 * t205;
t153 = t199 * t205 + t200 * t243;
t244 = t206 * t208;
t267 = -t272 * t152 - t274 * t153 - t271 * t244;
t266 = t273 * t150 + t275 * t151 - t272 * t247;
t265 = t273 * t152 + t275 * t153 - t272 * t244;
t264 = t275 * t150 + t276 * t151 - t274 * t247;
t263 = t275 * t152 + t276 * t153 - t274 * t244;
t262 = t271 * t210 + (t272 * t199 - t274 * t200) * t208;
t261 = t272 * t210 + (-t273 * t199 + t275 * t200) * t208;
t260 = t274 * t210 + (-t275 * t199 + t276 * t200) * t208;
t209 = cos(qJ(3));
t255 = t209 * pkin(3);
t239 = pkin(4) * t200;
t220 = qJ(5) * t208 + t210 * t239;
t229 = pkin(4) * t199;
t253 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t247 + t205 * t220 - t206 * t229;
t252 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t244 + t205 * t229 + t206 * t220;
t251 = Icges(2,4) * t205;
t250 = Icges(3,4) * t208;
t249 = Icges(3,4) * t210;
t207 = sin(qJ(3));
t248 = t205 * t207;
t245 = t206 * t207;
t242 = t207 * t210;
t241 = t209 * t210;
t240 = (-qJ(5) - rSges(6,3)) * t210 + (rSges(6,1) * t200 - rSges(6,2) * t199 + t239) * t208;
t237 = qJD(3) * t208;
t236 = qJD(4) * t208;
t235 = qJD(5) * t208;
t234 = V_base(5) * qJ(1) + V_base(1);
t230 = qJD(1) + V_base(3);
t192 = qJD(2) * t205 + V_base(4);
t161 = t206 * t237 + t192;
t228 = pkin(2) * t210 + pkin(6) * t208;
t191 = -qJD(2) * t206 + V_base(5);
t227 = rSges(3,1) * t210 - rSges(3,2) * t208;
t226 = Icges(3,1) * t210 - t250;
t225 = -Icges(3,2) * t208 + t249;
t224 = Icges(3,5) * t210 - Icges(3,6) * t208;
t160 = t205 * t237 + t191;
t185 = pkin(1) * t206 + pkin(5) * t205;
t223 = -V_base(4) * qJ(1) + V_base(6) * t185 + V_base(2);
t184 = pkin(1) * t205 - pkin(5) * t206;
t222 = V_base(4) * t184 - t185 * V_base(5) + t230;
t221 = pkin(7) * t208 + t210 * t255;
t166 = t228 * t205;
t190 = t208 * pkin(2) - t210 * pkin(6);
t219 = t191 * t190 + (-t166 - t184) * V_base(6) + t234;
t218 = (-Icges(3,3) * t206 + t205 * t224) * t191 + (Icges(3,3) * t205 + t206 * t224) * t192 + (Icges(3,5) * t208 + Icges(3,6) * t210) * V_base(6);
t167 = t228 * t206;
t217 = V_base(6) * t167 - t190 * t192 + t223;
t216 = t192 * t166 - t167 * t191 + t222;
t123 = -pkin(3) * t245 + t205 * t221;
t131 = -pkin(7) * t210 + t208 * t255;
t193 = -qJD(3) * t210 + V_base(6);
t215 = -t123 * t193 + t160 * t131 + t219;
t124 = pkin(3) * t248 + t206 * t221;
t214 = t193 * t124 - t131 * t161 + t217;
t213 = t161 * t123 - t124 * t160 + t216;
t144 = -Icges(3,6) * t206 + t205 * t225;
t145 = Icges(3,6) * t205 + t206 * t225;
t146 = -Icges(3,5) * t206 + t205 * t226;
t147 = Icges(3,5) * t205 + t206 * t226;
t187 = Icges(3,2) * t210 + t250;
t188 = Icges(3,1) * t208 + t249;
t212 = (-t145 * t208 + t147 * t210) * t192 + (-t144 * t208 + t146 * t210) * t191 + (-t187 * t208 + t188 * t210) * V_base(6);
t198 = Icges(2,4) * t206;
t189 = rSges(3,1) * t208 + rSges(3,2) * t210;
t183 = rSges(2,1) * t206 - rSges(2,2) * t205;
t182 = rSges(2,1) * t205 + rSges(2,2) * t206;
t181 = Icges(2,1) * t206 - t251;
t180 = Icges(2,1) * t205 + t198;
t179 = -Icges(2,2) * t205 + t198;
t178 = Icges(2,2) * t206 + t251;
t174 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t173 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t172 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t170 = V_base(6) + (-qJD(3) - qJD(4)) * t210;
t165 = t206 * t241 + t248;
t164 = t205 * t209 - t206 * t242;
t163 = t205 * t241 - t245;
t162 = -t205 * t242 - t206 * t209;
t157 = -rSges(4,3) * t210 + (rSges(4,1) * t209 - rSges(4,2) * t207) * t208;
t156 = -Icges(4,5) * t210 + (Icges(4,1) * t209 - Icges(4,4) * t207) * t208;
t155 = -Icges(4,6) * t210 + (Icges(4,4) * t209 - Icges(4,2) * t207) * t208;
t154 = -Icges(4,3) * t210 + (Icges(4,5) * t209 - Icges(4,6) * t207) * t208;
t149 = rSges(3,3) * t205 + t206 * t227;
t148 = -rSges(3,3) * t206 + t205 * t227;
t141 = -rSges(5,3) * t210 + (rSges(5,1) * t200 - rSges(5,2) * t199) * t208;
t133 = t206 * t236 + t161;
t132 = t205 * t236 + t160;
t129 = V_base(5) * rSges(2,3) - t182 * V_base(6) + t234;
t128 = t183 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t127 = t182 * V_base(4) - t183 * V_base(5) + t230;
t122 = rSges(4,1) * t165 + rSges(4,2) * t164 + rSges(4,3) * t244;
t121 = rSges(4,1) * t163 + rSges(4,2) * t162 + rSges(4,3) * t247;
t120 = Icges(4,1) * t165 + Icges(4,4) * t164 + Icges(4,5) * t244;
t119 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t247;
t118 = Icges(4,4) * t165 + Icges(4,2) * t164 + Icges(4,6) * t244;
t117 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t247;
t116 = Icges(4,5) * t165 + Icges(4,6) * t164 + Icges(4,3) * t244;
t115 = Icges(4,5) * t163 + Icges(4,6) * t162 + Icges(4,3) * t247;
t114 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t244;
t112 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t247;
t96 = t189 * t191 + (-t148 - t184) * V_base(6) + t234;
t95 = t149 * V_base(6) - t189 * t192 + t223;
t92 = t148 * t192 - t149 * t191 + t222;
t91 = -t121 * t193 + t157 * t160 + t219;
t90 = t122 * t193 - t157 * t161 + t217;
t89 = t121 * t161 - t122 * t160 + t216;
t88 = -t112 * t170 + t132 * t141 + t215;
t87 = t114 * t170 - t133 * t141 + t214;
t86 = t112 * t133 - t114 * t132 + t213;
t85 = t132 * t240 - t170 * t253 + t206 * t235 + t215;
t84 = -t133 * t240 + t170 * t252 + t205 * t235 + t214;
t83 = -qJD(5) * t210 - t132 * t252 + t133 * t253 + t213;
t1 = m(1) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(2) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t192 * (t205 * t218 + t206 * t212) / 0.2e1 + t191 * (t205 * t212 - t206 * t218) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t161 * ((t116 * t244 + t118 * t164 + t120 * t165) * t161 + (t115 * t244 + t117 * t164 + t119 * t165) * t160 + (t154 * t244 + t155 * t164 + t156 * t165) * t193) / 0.2e1 + t160 * ((t116 * t247 + t118 * t162 + t120 * t163) * t161 + (t115 * t247 + t117 * t162 + t119 * t163) * t160 + (t154 * t247 + t155 * t162 + t156 * t163) * t193) / 0.2e1 + t193 * ((-t115 * t160 - t116 * t161 - t154 * t193) * t210 + ((-t118 * t207 + t120 * t209) * t161 + (-t117 * t207 + t119 * t209) * t160 + (-t155 * t207 + t156 * t209) * t193) * t208) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t261 * t150 + t260 * t151 + t262 * t247) * t170 + (t265 * t150 + t263 * t151 + t267 * t247) * t133 + (t266 * t150 + t264 * t151 + t270 * t247) * t132) * t132 / 0.2e1 + ((t261 * t152 + t260 * t153 + t262 * t244) * t170 + (t265 * t152 + t263 * t153 + t267 * t244) * t133 + (t266 * t152 + t264 * t153 + t270 * t244) * t132) * t133 / 0.2e1 + ((-t270 * t132 - t267 * t133 - t262 * t170) * t210 + ((-t261 * t199 + t260 * t200) * t170 + (-t265 * t199 + t263 * t200) * t133 + (-t266 * t199 + t264 * t200) * t132) * t208) * t170 / 0.2e1 + ((-t178 * t205 + t180 * t206 + Icges(1,4)) * V_base(5) + (-t179 * t205 + t181 * t206 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t178 * t206 + t180 * t205 + Icges(1,2)) * V_base(5) + (t179 * t206 + t181 * t205 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t210 + t147 * t208) * t192 + (t144 * t210 + t146 * t208) * t191 + (t187 * t210 + t188 * t208 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t206 - Icges(2,6) * t205 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t205 + Icges(2,6) * t206 + Icges(1,6));
T = t1;
