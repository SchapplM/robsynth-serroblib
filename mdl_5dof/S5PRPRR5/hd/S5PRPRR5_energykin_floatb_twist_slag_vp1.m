% Calculate kinetic energy for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:21
% EndTime: 2019-12-05 15:53:24
% DurationCPUTime: 2.89s
% Computational Cost: add. (1418->333), mult. (1929->496), div. (0->0), fcn. (1879->10), ass. (0->149)
t212 = cos(pkin(9));
t256 = t212 * pkin(3);
t211 = sin(pkin(8));
t255 = Icges(2,4) * t211;
t215 = sin(qJ(2));
t254 = Icges(3,4) * t215;
t216 = cos(qJ(2));
t253 = Icges(3,4) * t216;
t210 = sin(pkin(9));
t252 = t211 * t210;
t251 = t211 * t215;
t250 = t211 * t216;
t213 = cos(pkin(8));
t249 = t213 * t210;
t248 = t213 * t215;
t247 = t213 * t216;
t231 = pkin(2) * t216 + qJ(3) * t215;
t166 = t231 * t211;
t184 = pkin(1) * t211 - pkin(5) * t213;
t245 = -t166 - t184;
t209 = pkin(9) + qJ(4);
t203 = cos(t209);
t244 = pkin(4) * t203;
t242 = qJD(3) * t215;
t241 = qJD(4) * t215;
t240 = qJD(5) * t215;
t239 = V_base(5) * qJ(1) + V_base(1);
t235 = qJD(1) + V_base(3);
t192 = qJD(2) * t211 + V_base(4);
t202 = sin(t209);
t234 = pkin(4) * t202;
t165 = t213 * t241 + t192;
t189 = t215 * pkin(2) - qJ(3) * t216;
t191 = -qJD(2) * t213 + V_base(5);
t233 = t189 * t191 + t213 * t242 + t239;
t232 = rSges(3,1) * t216 - rSges(3,2) * t215;
t230 = Icges(3,1) * t216 - t254;
t229 = -Icges(3,2) * t215 + t253;
t228 = Icges(3,5) * t216 - Icges(3,6) * t215;
t164 = t211 * t241 + t191;
t185 = pkin(1) * t213 + pkin(5) * t211;
t227 = -V_base(4) * qJ(1) + t185 * V_base(6) + V_base(2);
t226 = t184 * V_base(4) - t185 * V_base(5) + t235;
t225 = pkin(6) * t215 + t216 * t256;
t224 = pkin(7) * t215 + t216 * t244;
t167 = t231 * t213;
t223 = t167 * V_base(6) + t211 * t242 + t227;
t222 = (-Icges(3,3) * t213 + t211 * t228) * t191 + (Icges(3,3) * t211 + t213 * t228) * t192 + (Icges(3,5) * t215 + Icges(3,6) * t216) * V_base(6);
t221 = -qJD(3) * t216 + t166 * t192 + t226;
t119 = -pkin(3) * t249 + t211 * t225;
t128 = -pkin(6) * t216 + t215 * t256;
t220 = t191 * t128 + (-t119 + t245) * V_base(6) + t233;
t120 = pkin(3) * t252 + t213 * t225;
t219 = V_base(6) * t120 + (-t128 - t189) * t192 + t223;
t218 = t192 * t119 + (-t120 - t167) * t191 + t221;
t144 = -Icges(3,6) * t213 + t211 * t229;
t145 = Icges(3,6) * t211 + t213 * t229;
t146 = -Icges(3,5) * t213 + t211 * t230;
t147 = Icges(3,5) * t211 + t213 * t230;
t187 = Icges(3,2) * t216 + t254;
t188 = Icges(3,1) * t215 + t253;
t217 = (-t145 * t215 + t147 * t216) * t192 + (-t144 * t215 + t146 * t216) * t191 + (-t187 * t215 + t188 * t216) * V_base(6);
t205 = qJ(5) + t209;
t204 = Icges(2,4) * t213;
t199 = cos(t205);
t198 = sin(t205);
t193 = -qJD(4) * t216 + V_base(6);
t190 = rSges(3,1) * t215 + rSges(3,2) * t216;
t183 = rSges(2,1) * t213 - rSges(2,2) * t211;
t182 = rSges(2,1) * t211 + rSges(2,2) * t213;
t181 = Icges(2,1) * t213 - t255;
t180 = Icges(2,1) * t211 + t204;
t179 = -Icges(2,2) * t211 + t204;
t178 = Icges(2,2) * t213 + t255;
t175 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t174 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t173 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t170 = V_base(6) + (-qJD(4) - qJD(5)) * t216;
t163 = t212 * t247 + t252;
t162 = -t210 * t247 + t211 * t212;
t161 = t212 * t250 - t249;
t160 = -t210 * t250 - t212 * t213;
t157 = -rSges(4,3) * t216 + (rSges(4,1) * t212 - rSges(4,2) * t210) * t215;
t156 = t211 * rSges(3,3) + t213 * t232;
t155 = -rSges(3,3) * t213 + t211 * t232;
t154 = -Icges(4,5) * t216 + (Icges(4,1) * t212 - Icges(4,4) * t210) * t215;
t153 = -Icges(4,6) * t216 + (Icges(4,4) * t212 - Icges(4,2) * t210) * t215;
t152 = -Icges(4,3) * t216 + (Icges(4,5) * t212 - Icges(4,6) * t210) * t215;
t151 = t202 * t211 + t203 * t247;
t150 = -t202 * t247 + t203 * t211;
t149 = -t202 * t213 + t203 * t250;
t148 = -t202 * t250 - t203 * t213;
t141 = t198 * t211 + t199 * t247;
t140 = -t198 * t247 + t199 * t211;
t139 = -t198 * t213 + t199 * t250;
t138 = -t198 * t250 - t199 * t213;
t137 = -rSges(5,3) * t216 + (rSges(5,1) * t203 - rSges(5,2) * t202) * t215;
t136 = t213 * t240 + t165;
t135 = t211 * t240 + t164;
t134 = -Icges(5,5) * t216 + (Icges(5,1) * t203 - Icges(5,4) * t202) * t215;
t133 = -Icges(5,6) * t216 + (Icges(5,4) * t203 - Icges(5,2) * t202) * t215;
t132 = -Icges(5,3) * t216 + (Icges(5,5) * t203 - Icges(5,6) * t202) * t215;
t130 = V_base(5) * rSges(2,3) - t182 * V_base(6) + t239;
t129 = t183 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t127 = -rSges(6,3) * t216 + (rSges(6,1) * t199 - rSges(6,2) * t198) * t215;
t126 = -Icges(6,5) * t216 + (Icges(6,1) * t199 - Icges(6,4) * t198) * t215;
t125 = -Icges(6,6) * t216 + (Icges(6,4) * t199 - Icges(6,2) * t198) * t215;
t124 = -Icges(6,3) * t216 + (Icges(6,5) * t199 - Icges(6,6) * t198) * t215;
t123 = t182 * V_base(4) - t183 * V_base(5) + t235;
t121 = -pkin(7) * t216 + t215 * t244;
t118 = rSges(4,1) * t163 + rSges(4,2) * t162 + rSges(4,3) * t248;
t117 = rSges(4,1) * t161 + rSges(4,2) * t160 + rSges(4,3) * t251;
t116 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t248;
t115 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t251;
t114 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t248;
t113 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t251;
t112 = Icges(4,5) * t163 + Icges(4,6) * t162 + Icges(4,3) * t248;
t111 = Icges(4,5) * t161 + Icges(4,6) * t160 + Icges(4,3) * t251;
t109 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t248;
t108 = rSges(5,1) * t149 + rSges(5,2) * t148 + rSges(5,3) * t251;
t107 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t248;
t106 = Icges(5,1) * t149 + Icges(5,4) * t148 + Icges(5,5) * t251;
t105 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t248;
t104 = Icges(5,4) * t149 + Icges(5,2) * t148 + Icges(5,6) * t251;
t103 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t248;
t102 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t251;
t100 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t248;
t99 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t251;
t98 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t248;
t97 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t251;
t96 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t248;
t95 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t251;
t94 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t248;
t93 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t251;
t92 = t190 * t191 + (-t155 - t184) * V_base(6) + t239;
t91 = t156 * V_base(6) - t190 * t192 + t227;
t90 = t211 * t234 + t213 * t224;
t89 = t211 * t224 - t213 * t234;
t88 = t155 * t192 - t156 * t191 + t226;
t87 = t157 * t191 + (-t117 + t245) * V_base(6) + t233;
t86 = t118 * V_base(6) + (-t157 - t189) * t192 + t223;
t85 = t192 * t117 + (-t118 - t167) * t191 + t221;
t84 = -t108 * t193 + t137 * t164 + t220;
t83 = t109 * t193 - t137 * t165 + t219;
t82 = t108 * t165 - t164 * t109 + t218;
t81 = t121 * t164 + t127 * t135 - t170 * t99 - t193 * t89 + t220;
t80 = t100 * t170 - t121 * t165 - t127 * t136 + t193 * t90 + t219;
t79 = -t100 * t135 + t136 * t99 - t164 * t90 + t165 * t89 + t218;
t1 = m(1) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(3) * (t88 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + t165 * ((t103 * t248 + t105 * t150 + t107 * t151) * t165 + (t102 * t248 + t104 * t150 + t106 * t151) * t164 + (t132 * t248 + t133 * t150 + t134 * t151) * t193) / 0.2e1 + t164 * ((t103 * t251 + t105 * t148 + t107 * t149) * t165 + (t102 * t251 + t104 * t148 + t106 * t149) * t164 + (t132 * t251 + t133 * t148 + t134 * t149) * t193) / 0.2e1 + t193 * ((-t102 * t164 - t103 * t165 - t132 * t193) * t216 + ((-t105 * t202 + t107 * t203) * t165 + (-t104 * t202 + t106 * t203) * t164 + (-t133 * t202 + t134 * t203) * t193) * t215) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + t136 * ((t140 * t96 + t141 * t98 + t248 * t94) * t136 + (t140 * t95 + t141 * t97 + t248 * t93) * t135 + (t124 * t248 + t125 * t140 + t126 * t141) * t170) / 0.2e1 + t135 * ((t138 * t96 + t139 * t98 + t251 * t94) * t136 + (t138 * t95 + t139 * t97 + t251 * t93) * t135 + (t124 * t251 + t125 * t138 + t126 * t139) * t170) / 0.2e1 + t170 * ((-t124 * t170 - t135 * t93 - t136 * t94) * t216 + ((-t198 * t96 + t199 * t98) * t136 + (-t198 * t95 + t199 * t97) * t135 + (-t125 * t198 + t126 * t199) * t170) * t215) / 0.2e1 + (t211 * t217 - t213 * t222 + (t112 * t251 + t114 * t160 + t116 * t161) * t192 + (t111 * t251 + t113 * t160 + t115 * t161) * t191 + (t152 * t251 + t153 * t160 + t154 * t161) * V_base(6)) * t191 / 0.2e1 + (t211 * t222 + t213 * t217 + (t112 * t248 + t114 * t162 + t116 * t163) * t192 + (t111 * t248 + t113 * t162 + t115 * t163) * t191 + (t152 * t248 + t153 * t162 + t154 * t163) * V_base(6)) * t192 / 0.2e1 + ((-t178 * t211 + t180 * t213 + Icges(1,4)) * V_base(5) + (-t179 * t211 + t181 * t213 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t178 * t213 + t180 * t211 + Icges(1,2)) * V_base(5) + (t179 * t213 + t181 * t211 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t145 * t216 + t215 * t147) * t192 + (t144 * t216 + t215 * t146) * t191 + (-t111 * t191 - t112 * t192) * t216 + ((-t114 * t210 + t116 * t212) * t192 + (-t113 * t210 + t115 * t212) * t191) * t215 + (Icges(1,3) + Icges(2,3) + (t187 - t152) * t216 + (-t153 * t210 + t154 * t212 + t188) * t215) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t213 - Icges(2,6) * t211 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t211 + Icges(2,6) * t213 + Icges(1,6));
T = t1;
