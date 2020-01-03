% Calculate kinetic energy for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:12:57
% EndTime: 2019-12-31 21:12:59
% DurationCPUTime: 2.17s
% Computational Cost: add. (1428->280), mult. (1408->407), div. (0->0), fcn. (1246->10), ass. (0->144)
t261 = Icges(4,3) + Icges(5,3);
t195 = qJ(2) + qJ(3);
t184 = pkin(9) + t195;
t181 = sin(t184);
t182 = cos(t184);
t189 = sin(t195);
t190 = cos(t195);
t260 = Icges(4,5) * t190 + Icges(5,5) * t182 - Icges(4,6) * t189 - Icges(5,6) * t181;
t198 = sin(qJ(1));
t201 = cos(qJ(1));
t244 = Icges(5,4) * t182;
t219 = -Icges(5,2) * t181 + t244;
t111 = -Icges(5,6) * t201 + t219 * t198;
t112 = Icges(5,6) * t198 + t219 * t201;
t245 = Icges(5,4) * t181;
t222 = Icges(5,1) * t182 - t245;
t113 = -Icges(5,5) * t201 + t222 * t198;
t114 = Icges(5,5) * t198 + t222 * t201;
t246 = Icges(4,4) * t190;
t220 = -Icges(4,2) * t189 + t246;
t125 = -Icges(4,6) * t201 + t220 * t198;
t126 = Icges(4,6) * t198 + t220 * t201;
t247 = Icges(4,4) * t189;
t223 = Icges(4,1) * t190 - t247;
t127 = -Icges(4,5) * t201 + t223 * t198;
t128 = Icges(4,5) * t198 + t223 * t201;
t147 = Icges(5,2) * t182 + t245;
t148 = Icges(5,1) * t181 + t244;
t153 = Icges(4,2) * t190 + t247;
t154 = Icges(4,1) * t189 + t246;
t157 = V_base(5) + (-qJD(2) - qJD(3)) * t201;
t180 = qJD(2) * t198 + V_base(4);
t158 = qJD(3) * t198 + t180;
t185 = V_base(6) + qJD(1);
t259 = (-t147 * t181 + t148 * t182 - t153 * t189 + t154 * t190) * t185 + (-t112 * t181 + t114 * t182 - t126 * t189 + t128 * t190) * t158 + (-t111 * t181 + t113 * t182 - t125 * t189 + t127 * t190) * t157;
t258 = (Icges(4,5) * t189 + Icges(5,5) * t181 + Icges(4,6) * t190 + Icges(5,6) * t182) * t185 + (t261 * t198 + t260 * t201) * t158 + (t260 * t198 - t261 * t201) * t157;
t197 = sin(qJ(2));
t254 = pkin(2) * t197;
t253 = pkin(3) * t189;
t200 = cos(qJ(2));
t252 = t200 * pkin(2);
t250 = Icges(2,4) * t198;
t249 = Icges(3,4) * t197;
t248 = Icges(3,4) * t200;
t243 = t181 * t198;
t242 = t181 * t201;
t196 = sin(qJ(5));
t241 = t196 * t198;
t240 = t196 * t201;
t199 = cos(qJ(5));
t239 = t198 * t199;
t238 = t199 * t201;
t119 = -pkin(7) * t201 + t252 * t198;
t177 = t198 * pkin(1) - t201 * pkin(6);
t237 = -t119 - t177;
t236 = pkin(3) * t190;
t234 = qJD(5) * t181;
t233 = V_base(5) * pkin(5) + V_base(1);
t100 = -qJ(4) * t201 + t236 * t198;
t230 = -t100 + t237;
t179 = -qJD(2) * t201 + V_base(5);
t229 = t179 * t254 + t233;
t228 = pkin(4) * t182 + pkin(8) * t181;
t227 = rSges(3,1) * t200 - rSges(3,2) * t197;
t226 = rSges(4,1) * t190 - rSges(4,2) * t189;
t225 = rSges(5,1) * t182 - rSges(5,2) * t181;
t224 = Icges(3,1) * t200 - t249;
t221 = -Icges(3,2) * t197 + t248;
t218 = Icges(3,5) * t200 - Icges(3,6) * t197;
t215 = qJD(4) * t198 + t157 * t253 + t229;
t178 = t201 * pkin(1) + t198 * pkin(6);
t214 = -V_base(4) * pkin(5) + t185 * t178 + V_base(2);
t213 = V_base(4) * t177 - t178 * V_base(5) + V_base(3);
t210 = (-Icges(3,3) * t201 + t218 * t198) * t179 + (Icges(3,3) * t198 + t218 * t201) * t180 + (Icges(3,5) * t197 + Icges(3,6) * t200) * t185;
t120 = pkin(7) * t198 + t252 * t201;
t209 = t180 * t119 - t120 * t179 + t213;
t208 = t185 * t120 - t180 * t254 + t214;
t207 = t158 * t100 + t209;
t101 = qJ(4) * t198 + t236 * t201;
t206 = -qJD(4) * t201 + t185 * t101 + t208;
t135 = -Icges(3,6) * t201 + t221 * t198;
t136 = Icges(3,6) * t198 + t221 * t201;
t137 = -Icges(3,5) * t201 + t224 * t198;
t138 = Icges(3,5) * t198 + t224 * t201;
t168 = Icges(3,2) * t200 + t249;
t171 = Icges(3,1) * t197 + t248;
t203 = (-t136 * t197 + t138 * t200) * t180 + (-t135 * t197 + t137 * t200) * t179 + (-t168 * t197 + t171 * t200) * t185;
t191 = Icges(2,4) * t201;
t176 = rSges(2,1) * t201 - rSges(2,2) * t198;
t175 = rSges(2,1) * t198 + rSges(2,2) * t201;
t174 = rSges(3,1) * t197 + rSges(3,2) * t200;
t173 = Icges(2,1) * t201 - t250;
t172 = Icges(2,1) * t198 + t191;
t170 = -Icges(2,2) * t198 + t191;
t169 = Icges(2,2) * t201 + t250;
t164 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t163 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t162 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t156 = -qJD(5) * t182 + t185;
t155 = rSges(4,1) * t189 + rSges(4,2) * t190;
t150 = pkin(4) * t181 - pkin(8) * t182;
t149 = rSges(5,1) * t181 + rSges(5,2) * t182;
t145 = t182 * t238 + t241;
t144 = -t182 * t240 + t239;
t143 = t182 * t239 - t240;
t142 = -t182 * t241 - t238;
t140 = rSges(3,3) * t198 + t227 * t201;
t139 = -rSges(3,3) * t201 + t227 * t198;
t132 = t228 * t201;
t131 = t228 * t198;
t130 = rSges(4,3) * t198 + t226 * t201;
t129 = -rSges(4,3) * t201 + t226 * t198;
t122 = t201 * t234 + t158;
t121 = t198 * t234 + t157;
t118 = rSges(5,3) * t198 + t225 * t201;
t117 = -rSges(5,3) * t201 + t225 * t198;
t116 = V_base(5) * rSges(2,3) - t175 * t185 + t233;
t115 = t176 * t185 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t108 = t175 * V_base(4) - t176 * V_base(5) + V_base(3);
t106 = -rSges(6,3) * t182 + (rSges(6,1) * t199 - rSges(6,2) * t196) * t181;
t105 = -Icges(6,5) * t182 + (Icges(6,1) * t199 - Icges(6,4) * t196) * t181;
t104 = -Icges(6,6) * t182 + (Icges(6,4) * t199 - Icges(6,2) * t196) * t181;
t103 = -Icges(6,3) * t182 + (Icges(6,5) * t199 - Icges(6,6) * t196) * t181;
t97 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t242;
t96 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t243;
t95 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t242;
t94 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t243;
t93 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t242;
t92 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t243;
t91 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t242;
t90 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t243;
t89 = t174 * t179 + (-t139 - t177) * t185 + t233;
t88 = t140 * t185 - t174 * t180 + t214;
t87 = t139 * t180 - t140 * t179 + t213;
t86 = t155 * t157 + (-t129 + t237) * t185 + t229;
t85 = t130 * t185 - t155 * t158 + t208;
t84 = t129 * t158 - t130 * t157 + t209;
t83 = t149 * t157 + (-t117 + t230) * t185 + t215;
t82 = t118 * t185 + (-t149 - t253) * t158 + t206;
t81 = t117 * t158 + (-t101 - t118) * t157 + t207;
t80 = t106 * t121 + t150 * t157 - t156 * t96 + (-t131 + t230) * t185 + t215;
t79 = -t106 * t122 + t132 * t185 + t156 * t97 + (-t150 - t253) * t158 + t206;
t78 = -t121 * t97 + t122 * t96 + t131 * t158 + (-t101 - t132) * t157 + t207;
t1 = m(1) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t180 * (t210 * t198 + t203 * t201) / 0.2e1 + t179 * (t203 * t198 - t210 * t201) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t122 * ((t144 * t93 + t145 * t95 + t91 * t242) * t122 + (t144 * t92 + t145 * t94 + t90 * t242) * t121 + (t103 * t242 + t104 * t144 + t105 * t145) * t156) / 0.2e1 + t121 * ((t142 * t93 + t143 * t95 + t91 * t243) * t122 + (t142 * t92 + t143 * t94 + t90 * t243) * t121 + (t103 * t243 + t104 * t142 + t105 * t143) * t156) / 0.2e1 + t156 * ((-t103 * t156 - t90 * t121 - t91 * t122) * t182 + ((-t196 * t93 + t199 * t95) * t122 + (-t196 * t92 + t199 * t94) * t121 + (-t104 * t196 + t105 * t199) * t156) * t181) / 0.2e1 + (t259 * t198 - t258 * t201) * t157 / 0.2e1 + (t258 * t198 + t259 * t201) * t158 / 0.2e1 + ((-t169 * t198 + t172 * t201 + Icges(1,4)) * V_base(5) + (-t170 * t198 + t173 * t201 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t169 * t201 + t172 * t198 + Icges(1,2)) * V_base(5) + (t170 * t201 + t173 * t198 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t200 + t138 * t197) * t180 + (t135 * t200 + t137 * t197) * t179 + (t112 * t182 + t114 * t181 + t126 * t190 + t128 * t189) * t158 + (t111 * t182 + t113 * t181 + t125 * t190 + t127 * t189) * t157 + (t147 * t182 + t148 * t181 + t153 * t190 + t154 * t189 + t168 * t200 + t171 * t197 + Icges(2,3)) * t185) * t185 / 0.2e1 + t185 * V_base(4) * (Icges(2,5) * t201 - Icges(2,6) * t198) + V_base(5) * t185 * (Icges(2,5) * t198 + Icges(2,6) * t201) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
