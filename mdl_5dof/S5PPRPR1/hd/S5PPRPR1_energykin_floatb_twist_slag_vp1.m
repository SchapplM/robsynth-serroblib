% Calculate kinetic energy for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:49
% EndTime: 2019-12-05 15:00:51
% DurationCPUTime: 2.39s
% Computational Cost: add. (1348->312), mult. (1558->447), div. (0->0), fcn. (1452->10), ass. (0->150)
t201 = sin(pkin(7));
t204 = cos(pkin(7));
t260 = Icges(2,5) * t204 - Icges(2,6) * t201 + Icges(1,5);
t259 = Icges(2,5) * t201 + Icges(2,6) * t204 + Icges(1,6);
t200 = sin(pkin(8));
t256 = pkin(2) * t200;
t203 = cos(pkin(8));
t255 = pkin(2) * t203;
t202 = cos(pkin(9));
t254 = pkin(4) * t202;
t253 = Icges(2,4) * t201;
t252 = Icges(3,4) * t200;
t251 = Icges(3,4) * t203;
t198 = pkin(8) + qJ(3);
t192 = sin(t198);
t250 = Icges(4,4) * t192;
t194 = cos(t198);
t249 = Icges(4,4) * t194;
t248 = t192 * t201;
t247 = t192 * t204;
t246 = t194 * t201;
t245 = t194 * t204;
t199 = sin(pkin(9));
t244 = t199 * t201;
t243 = t199 * t204;
t242 = t201 * t202;
t241 = t202 * t204;
t123 = -pkin(5) * t204 + t201 * t255;
t178 = t201 * pkin(1) - t204 * qJ(2);
t238 = -t123 - t178;
t237 = qJD(4) * t192;
t236 = qJD(5) * t192;
t235 = V_base(5) * qJ(1) + V_base(1);
t231 = qJD(1) + V_base(3);
t224 = pkin(3) * t194 + qJ(4) * t192;
t149 = t224 * t201;
t230 = -t149 + t238;
t185 = qJD(3) * t201 + V_base(4);
t229 = qJD(2) * t201 + t235;
t228 = V_base(4) * t178 + t231;
t227 = V_base(5) * t256 + t229;
t184 = -qJD(3) * t204 + V_base(5);
t226 = rSges(3,1) * t203 - rSges(3,2) * t200;
t225 = rSges(4,1) * t194 - rSges(4,2) * t192;
t223 = Icges(3,1) * t203 - t252;
t222 = Icges(4,1) * t194 - t250;
t221 = -Icges(3,2) * t200 + t251;
t220 = -Icges(4,2) * t192 + t249;
t219 = Icges(3,5) * t203 - Icges(3,6) * t200;
t218 = Icges(4,5) * t194 - Icges(4,6) * t192;
t180 = t204 * pkin(1) + t201 * qJ(2);
t217 = -qJD(2) * t204 + V_base(6) * t180 + V_base(2);
t160 = pkin(3) * t192 - qJ(4) * t194;
t216 = t184 * t160 + t204 * t237 + t227;
t215 = pkin(6) * t192 + t194 * t254;
t214 = (-Icges(4,3) * t204 + t201 * t218) * t184 + (Icges(4,3) * t201 + t204 * t218) * t185 + (Icges(4,5) * t192 + Icges(4,6) * t194) * V_base(6);
t124 = pkin(5) * t201 + t204 * t255;
t213 = V_base(4) * t123 + (-t124 - t180) * V_base(5) + t228;
t212 = (-Icges(3,3) * t204 + t201 * t219) * V_base(5) + (Icges(3,3) * t201 + t204 * t219) * V_base(4) + (Icges(3,5) * t200 + Icges(3,6) * t203) * V_base(6);
t211 = V_base(6) * t124 + (-qJ(1) - t256) * V_base(4) + t217;
t210 = -qJD(4) * t194 + t185 * t149 + t213;
t150 = t224 * t204;
t209 = V_base(6) * t150 + t201 * t237 + t211;
t127 = -Icges(4,6) * t204 + t201 * t220;
t128 = Icges(4,6) * t201 + t204 * t220;
t129 = -Icges(4,5) * t204 + t201 * t222;
t130 = Icges(4,5) * t201 + t204 * t222;
t158 = Icges(4,2) * t194 + t250;
t159 = Icges(4,1) * t192 + t249;
t208 = (-t128 * t192 + t130 * t194) * t185 + (-t127 * t192 + t129 * t194) * t184 + (-t158 * t192 + t159 * t194) * V_base(6);
t143 = -Icges(3,6) * t204 + t201 * t221;
t144 = Icges(3,6) * t201 + t204 * t221;
t145 = -Icges(3,5) * t204 + t201 * t223;
t146 = Icges(3,5) * t201 + t204 * t223;
t170 = Icges(3,2) * t203 + t252;
t173 = Icges(3,1) * t200 + t251;
t207 = (-t144 * t200 + t146 * t203) * V_base(4) + (-t143 * t200 + t145 * t203) * V_base(5) + (-t170 * t200 + t173 * t203) * V_base(6);
t197 = pkin(9) + qJ(5);
t195 = Icges(2,4) * t204;
t193 = cos(t197);
t191 = sin(t197);
t181 = rSges(2,1) * t204 - rSges(2,2) * t201;
t179 = rSges(2,1) * t201 + rSges(2,2) * t204;
t177 = rSges(3,1) * t200 + rSges(3,2) * t203;
t176 = -qJD(5) * t194 + V_base(6);
t175 = Icges(2,1) * t204 - t253;
t174 = Icges(2,1) * t201 + t195;
t172 = -Icges(2,2) * t201 + t195;
t171 = Icges(2,2) * t204 + t253;
t166 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t165 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t164 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = rSges(4,1) * t192 + rSges(4,2) * t194;
t156 = t204 * t236 + t185;
t155 = t201 * t236 + t184;
t154 = t194 * t241 + t244;
t153 = -t194 * t243 + t242;
t152 = t194 * t242 - t243;
t151 = -t194 * t244 - t241;
t148 = rSges(3,3) * t201 + t204 * t226;
t147 = -rSges(3,3) * t204 + t201 * t226;
t140 = t191 * t201 + t193 * t245;
t139 = -t191 * t245 + t193 * t201;
t138 = -t191 * t204 + t193 * t246;
t137 = -t191 * t246 - t193 * t204;
t134 = rSges(4,3) * t201 + t204 * t225;
t133 = -rSges(4,3) * t204 + t201 * t225;
t132 = V_base(5) * rSges(2,3) - t179 * V_base(6) + t235;
t131 = t181 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t122 = -rSges(5,3) * t194 + (rSges(5,1) * t202 - rSges(5,2) * t199) * t192;
t121 = -Icges(5,5) * t194 + (Icges(5,1) * t202 - Icges(5,4) * t199) * t192;
t120 = -Icges(5,6) * t194 + (Icges(5,4) * t202 - Icges(5,2) * t199) * t192;
t119 = -Icges(5,3) * t194 + (Icges(5,5) * t202 - Icges(5,6) * t199) * t192;
t115 = -rSges(6,3) * t194 + (rSges(6,1) * t193 - rSges(6,2) * t191) * t192;
t114 = -Icges(6,5) * t194 + (Icges(6,1) * t193 - Icges(6,4) * t191) * t192;
t113 = -Icges(6,6) * t194 + (Icges(6,4) * t193 - Icges(6,2) * t191) * t192;
t112 = -Icges(6,3) * t194 + (Icges(6,5) * t193 - Icges(6,6) * t191) * t192;
t111 = t179 * V_base(4) - t181 * V_base(5) + t231;
t110 = -pkin(6) * t194 + t192 * t254;
t109 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t247;
t108 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t248;
t107 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t247;
t106 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t248;
t105 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t247;
t104 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t248;
t103 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t247;
t102 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t248;
t101 = pkin(4) * t244 + t204 * t215;
t100 = -pkin(4) * t243 + t201 * t215;
t99 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t247;
t98 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t248;
t97 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t247;
t96 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t248;
t95 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t247;
t94 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t248;
t93 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t247;
t92 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t248;
t91 = t177 * V_base(5) + (-t147 - t178) * V_base(6) + t229;
t90 = t148 * V_base(6) + (-qJ(1) - t177) * V_base(4) + t217;
t89 = t147 * V_base(4) + (-t148 - t180) * V_base(5) + t228;
t88 = t161 * t184 + (-t133 + t238) * V_base(6) + t227;
t87 = t134 * V_base(6) - t161 * t185 + t211;
t86 = t133 * t185 - t134 * t184 + t213;
t85 = t122 * t184 + (-t108 + t230) * V_base(6) + t216;
t84 = t109 * V_base(6) + (-t122 - t160) * t185 + t209;
t83 = t108 * t185 + (-t109 - t150) * t184 + t210;
t82 = t110 * t184 + t115 * t155 - t176 * t98 + (-t100 + t230) * V_base(6) + t216;
t81 = t101 * V_base(6) - t115 * t156 + t176 * t99 + (-t110 - t160) * t185 + t209;
t80 = t100 * t185 - t155 * t99 + t156 * t98 + (-t101 - t150) * t184 + t210;
t1 = m(1) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t156 * ((t139 * t95 + t140 * t97 + t93 * t247) * t156 + (t139 * t94 + t140 * t96 + t247 * t92) * t155 + (t112 * t247 + t113 * t139 + t114 * t140) * t176) / 0.2e1 + t155 * ((t137 * t95 + t138 * t97 + t248 * t93) * t156 + (t137 * t94 + t138 * t96 + t92 * t248) * t155 + (t112 * t248 + t113 * t137 + t114 * t138) * t176) / 0.2e1 + t176 * ((-t112 * t176 - t155 * t92 - t156 * t93) * t194 + ((-t191 * t95 + t193 * t97) * t156 + (-t191 * t94 + t193 * t96) * t155 + (-t113 * t191 + t114 * t193) * t176) * t192) / 0.2e1 + (t208 * t201 - t214 * t204 + (t103 * t248 + t105 * t151 + t107 * t152) * t185 + (t102 * t248 + t104 * t151 + t106 * t152) * t184 + (t119 * t248 + t120 * t151 + t121 * t152) * V_base(6)) * t184 / 0.2e1 + (t214 * t201 + t208 * t204 + (t103 * t247 + t105 * t153 + t107 * t154) * t185 + (t102 * t247 + t104 * t153 + t106 * t154) * t184 + (t119 * t247 + t120 * t153 + t121 * t154) * V_base(6)) * t185 / 0.2e1 + (t212 * t201 + t207 * t204 + t260 * V_base(6) + (-t171 * t201 + t174 * t204 + Icges(1,4)) * V_base(5) + (-t172 * t201 + t175 * t204 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t207 * t201 - t212 * t204 + t259 * V_base(6) + (t171 * t204 + t174 * t201 + Icges(1,2)) * V_base(5) + (t172 * t204 + t175 * t201 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t128 * t194 + t130 * t192) * t185 + (t127 * t194 + t129 * t192) * t184 + (-t102 * t184 - t103 * t185) * t194 + ((-t105 * t199 + t107 * t202) * t185 + (-t104 * t199 + t106 * t202) * t184) * t192 + (t143 * t203 + t145 * t200 + t259) * V_base(5) + (t144 * t203 + t146 * t200 + t260) * V_base(4) + (t170 * t203 + t173 * t200 + Icges(1,3) + Icges(2,3) + (t158 - t119) * t194 + (-t120 * t199 + t121 * t202 + t159) * t192) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
