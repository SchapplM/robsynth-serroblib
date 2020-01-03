% Calculate kinetic energy for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:22
% EndTime: 2019-12-31 18:29:24
% DurationCPUTime: 2.45s
% Computational Cost: add. (1393->312), mult. (1558->451), div. (0->0), fcn. (1452->10), ass. (0->151)
t201 = sin(pkin(8));
t256 = pkin(2) * t201;
t203 = cos(pkin(8));
t255 = pkin(2) * t203;
t202 = cos(pkin(9));
t254 = pkin(4) * t202;
t206 = sin(qJ(1));
t253 = Icges(2,4) * t206;
t252 = Icges(3,4) * t201;
t251 = Icges(3,4) * t203;
t199 = pkin(8) + qJ(3);
t190 = sin(t199);
t250 = Icges(4,4) * t190;
t192 = cos(t199);
t249 = Icges(4,4) * t192;
t248 = t190 * t206;
t207 = cos(qJ(1));
t247 = t190 * t207;
t246 = t192 * t207;
t200 = sin(pkin(9));
t245 = t200 * t207;
t244 = t202 * t207;
t198 = pkin(9) + qJ(5);
t189 = sin(t198);
t243 = t206 * t189;
t191 = cos(t198);
t242 = t206 * t191;
t241 = t206 * t200;
t240 = t206 * t202;
t125 = -pkin(6) * t207 + t206 * t255;
t180 = t206 * pkin(1) - qJ(2) * t207;
t237 = -t125 - t180;
t236 = qJD(4) * t190;
t235 = qJD(5) * t190;
t234 = V_base(4) * t180 + V_base(3);
t233 = V_base(5) * pkin(5) + V_base(1);
t225 = pkin(3) * t192 + qJ(4) * t190;
t149 = t225 * t206;
t230 = -t149 + t237;
t185 = qJD(3) * t206 + V_base(4);
t193 = V_base(6) + qJD(1);
t229 = qJD(2) * t206 + t233;
t228 = V_base(5) * t256 + t229;
t184 = -qJD(3) * t207 + V_base(5);
t227 = rSges(3,1) * t203 - rSges(3,2) * t201;
t226 = rSges(4,1) * t192 - rSges(4,2) * t190;
t224 = Icges(3,1) * t203 - t252;
t223 = Icges(4,1) * t192 - t250;
t222 = -Icges(3,2) * t201 + t251;
t221 = -Icges(4,2) * t190 + t249;
t220 = Icges(3,5) * t203 - Icges(3,6) * t201;
t219 = Icges(4,5) * t192 - Icges(4,6) * t190;
t182 = pkin(1) * t207 + t206 * qJ(2);
t218 = -qJD(2) * t207 + t193 * t182 + V_base(2);
t161 = pkin(3) * t190 - qJ(4) * t192;
t217 = t184 * t161 + t207 * t236 + t228;
t216 = (-Icges(4,3) * t207 + t206 * t219) * t184 + (Icges(4,3) * t206 + t207 * t219) * t185 + (Icges(4,5) * t190 + Icges(4,6) * t192) * t193;
t215 = pkin(7) * t190 + t192 * t254;
t126 = pkin(6) * t206 + t207 * t255;
t214 = V_base(4) * t125 + (-t126 - t182) * V_base(5) + t234;
t213 = (-Icges(3,3) * t207 + t206 * t220) * V_base(5) + (Icges(3,3) * t206 + t207 * t220) * V_base(4) + (Icges(3,5) * t201 + Icges(3,6) * t203) * t193;
t212 = t193 * t126 + (-pkin(5) - t256) * V_base(4) + t218;
t211 = -qJD(4) * t192 + t185 * t149 + t214;
t150 = t225 * t207;
t210 = t193 * t150 + t206 * t236 + t212;
t130 = -Icges(4,6) * t207 + t206 * t221;
t131 = Icges(4,6) * t206 + t207 * t221;
t132 = -Icges(4,5) * t207 + t206 * t223;
t133 = Icges(4,5) * t206 + t207 * t223;
t159 = Icges(4,2) * t192 + t250;
t160 = Icges(4,1) * t190 + t249;
t209 = (-t131 * t190 + t133 * t192) * t185 + (-t130 * t190 + t132 * t192) * t184 + (-t159 * t190 + t160 * t192) * t193;
t143 = -Icges(3,6) * t207 + t206 * t222;
t144 = Icges(3,6) * t206 + t207 * t222;
t145 = -Icges(3,5) * t207 + t206 * t224;
t146 = Icges(3,5) * t206 + t207 * t224;
t169 = Icges(3,2) * t203 + t252;
t170 = Icges(3,1) * t201 + t251;
t208 = (-t144 * t201 + t146 * t203) * V_base(4) + (-t143 * t201 + t145 * t203) * V_base(5) + (-t169 * t201 + t170 * t203) * t193;
t196 = Icges(2,4) * t207;
t183 = rSges(2,1) * t207 - t206 * rSges(2,2);
t181 = t206 * rSges(2,1) + rSges(2,2) * t207;
t177 = Icges(2,1) * t207 - t253;
t176 = Icges(2,1) * t206 + t196;
t175 = -Icges(2,2) * t206 + t196;
t174 = Icges(2,2) * t207 + t253;
t173 = Icges(2,5) * t207 - Icges(2,6) * t206;
t172 = Icges(2,5) * t206 + Icges(2,6) * t207;
t171 = rSges(3,1) * t201 + rSges(3,2) * t203;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t164 = -qJD(5) * t192 + t193;
t162 = rSges(4,1) * t190 + rSges(4,2) * t192;
t156 = t207 * t235 + t185;
t155 = t206 * t235 + t184;
t154 = t192 * t244 + t241;
t153 = -t192 * t245 + t240;
t152 = t192 * t240 - t245;
t151 = -t192 * t241 - t244;
t148 = t206 * rSges(3,3) + t207 * t227;
t147 = -rSges(3,3) * t207 + t206 * t227;
t140 = t191 * t246 + t243;
t139 = -t189 * t246 + t242;
t138 = -t189 * t207 + t192 * t242;
t137 = -t191 * t207 - t192 * t243;
t135 = t206 * rSges(4,3) + t207 * t226;
t134 = -rSges(4,3) * t207 + t206 * t226;
t124 = V_base(5) * rSges(2,3) - t181 * t193 + t233;
t123 = t183 * t193 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t122 = -rSges(5,3) * t192 + (rSges(5,1) * t202 - rSges(5,2) * t200) * t190;
t121 = -Icges(5,5) * t192 + (Icges(5,1) * t202 - Icges(5,4) * t200) * t190;
t120 = -Icges(5,6) * t192 + (Icges(5,4) * t202 - Icges(5,2) * t200) * t190;
t119 = -Icges(5,3) * t192 + (Icges(5,5) * t202 - Icges(5,6) * t200) * t190;
t117 = t181 * V_base(4) - t183 * V_base(5) + V_base(3);
t114 = -rSges(6,3) * t192 + (rSges(6,1) * t191 - rSges(6,2) * t189) * t190;
t113 = -Icges(6,5) * t192 + (Icges(6,1) * t191 - Icges(6,4) * t189) * t190;
t112 = -Icges(6,6) * t192 + (Icges(6,4) * t191 - Icges(6,2) * t189) * t190;
t111 = -Icges(6,3) * t192 + (Icges(6,5) * t191 - Icges(6,6) * t189) * t190;
t110 = -pkin(7) * t192 + t190 * t254;
t109 = t154 * rSges(5,1) + t153 * rSges(5,2) + rSges(5,3) * t247;
t108 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t248;
t107 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t247;
t106 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t248;
t105 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t247;
t104 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t248;
t103 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t247;
t102 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t248;
t101 = pkin(4) * t241 + t207 * t215;
t100 = -pkin(4) * t245 + t206 * t215;
t99 = t140 * rSges(6,1) + t139 * rSges(6,2) + rSges(6,3) * t247;
t98 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t248;
t97 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t247;
t96 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t248;
t95 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t247;
t94 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t248;
t93 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t247;
t92 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t248;
t91 = t171 * V_base(5) + (-t147 - t180) * t193 + t229;
t90 = t193 * t148 + (-pkin(5) - t171) * V_base(4) + t218;
t89 = t147 * V_base(4) + (-t148 - t182) * V_base(5) + t234;
t88 = t162 * t184 + (-t134 + t237) * t193 + t228;
t87 = t193 * t135 - t185 * t162 + t212;
t86 = t134 * t185 - t135 * t184 + t214;
t85 = t122 * t184 + (-t108 + t230) * t193 + t217;
t84 = t193 * t109 + (-t122 - t161) * t185 + t210;
t83 = t108 * t185 + (-t109 - t150) * t184 + t211;
t82 = t110 * t184 + t114 * t155 - t164 * t98 + (-t100 + t230) * t193 + t217;
t81 = t193 * t101 - t156 * t114 + t164 * t99 + (-t110 - t161) * t185 + t210;
t80 = t100 * t185 - t155 * t99 + t156 * t98 + (-t101 - t150) * t184 + t211;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t156 * ((t139 * t95 + t140 * t97 + t93 * t247) * t156 + (t139 * t94 + t140 * t96 + t247 * t92) * t155 + (t111 * t247 + t139 * t112 + t140 * t113) * t164) / 0.2e1 + t155 * ((t137 * t95 + t138 * t97 + t248 * t93) * t156 + (t137 * t94 + t138 * t96 + t92 * t248) * t155 + (t111 * t248 + t112 * t137 + t113 * t138) * t164) / 0.2e1 + t164 * ((-t111 * t164 - t92 * t155 - t93 * t156) * t192 + ((-t189 * t95 + t191 * t97) * t156 + (-t189 * t94 + t191 * t96) * t155 + (-t112 * t189 + t113 * t191) * t164) * t190) / 0.2e1 + (t209 * t206 - t216 * t207 + (t103 * t248 + t105 * t151 + t107 * t152) * t185 + (t102 * t248 + t104 * t151 + t106 * t152) * t184 + (t119 * t248 + t120 * t151 + t121 * t152) * t193) * t184 / 0.2e1 + (t216 * t206 + t209 * t207 + (t103 * t247 + t153 * t105 + t154 * t107) * t185 + (t102 * t247 + t153 * t104 + t154 * t106) * t184 + (t119 * t247 + t153 * t120 + t154 * t121) * t193) * t185 / 0.2e1 + (t173 * t193 + t213 * t206 + t208 * t207 + (-t206 * t174 + t176 * t207 + Icges(1,4)) * V_base(5) + (-t206 * t175 + t177 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t172 * t193 + t208 * t206 - t213 * t207 + (t174 * t207 + t206 * t176 + Icges(1,2)) * V_base(5) + (t175 * t207 + t206 * t177 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t131 * t192 + t133 * t190) * t185 + (t130 * t192 + t132 * t190) * t184 + (-t102 * t184 - t103 * t185) * t192 + ((-t105 * t200 + t107 * t202) * t185 + (-t104 * t200 + t106 * t202) * t184) * t190 + (t143 * t203 + t145 * t201 + t172) * V_base(5) + (t144 * t203 + t146 * t201 + t173) * V_base(4) + (t169 * t203 + t170 * t201 + Icges(2,3) + (t159 - t119) * t192 + (-t120 * t200 + t121 * t202 + t160) * t190) * t193) * t193 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
