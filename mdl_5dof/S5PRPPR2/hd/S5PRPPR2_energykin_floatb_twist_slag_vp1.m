% Calculate kinetic energy for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:53
% EndTime: 2019-12-05 15:23:56
% DurationCPUTime: 2.62s
% Computational Cost: add. (1370->306), mult. (1580->427), div. (0->0), fcn. (1474->10), ass. (0->147)
t265 = Icges(3,3) + Icges(4,3);
t198 = qJ(2) + pkin(8);
t192 = sin(t198);
t194 = cos(t198);
t205 = sin(qJ(2));
t206 = cos(qJ(2));
t264 = Icges(3,5) * t206 + Icges(4,5) * t194 - Icges(3,6) * t205 - Icges(4,6) * t192;
t200 = sin(pkin(7));
t202 = cos(pkin(7));
t249 = Icges(4,4) * t194;
t221 = -Icges(4,2) * t192 + t249;
t127 = -Icges(4,6) * t202 + t200 * t221;
t128 = Icges(4,6) * t200 + t202 * t221;
t250 = Icges(4,4) * t192;
t223 = Icges(4,1) * t194 - t250;
t129 = -Icges(4,5) * t202 + t200 * t223;
t130 = Icges(4,5) * t200 + t202 * t223;
t251 = Icges(3,4) * t206;
t222 = -Icges(3,2) * t205 + t251;
t143 = -Icges(3,6) * t202 + t200 * t222;
t144 = Icges(3,6) * t200 + t202 * t222;
t252 = Icges(3,4) * t205;
t224 = Icges(3,1) * t206 - t252;
t145 = -Icges(3,5) * t202 + t200 * t224;
t146 = Icges(3,5) * t200 + t202 * t224;
t158 = Icges(4,2) * t194 + t250;
t159 = Icges(4,1) * t192 + t249;
t182 = Icges(3,2) * t206 + t252;
t183 = Icges(3,1) * t205 + t251;
t185 = -qJD(2) * t202 + V_base(5);
t186 = qJD(2) * t200 + V_base(4);
t261 = (-t158 * t192 + t159 * t194 - t182 * t205 + t183 * t206) * V_base(6) + (-t128 * t192 + t130 * t194 - t144 * t205 + t146 * t206) * t186 + (-t127 * t192 + t129 * t194 - t143 * t205 + t145 * t206) * t185;
t260 = (Icges(3,5) * t205 + Icges(4,5) * t192 + Icges(3,6) * t206 + Icges(4,6) * t194) * V_base(6) + (t265 * t200 + t264 * t202) * t186 + (t264 * t200 - t265 * t202) * t185;
t257 = pkin(2) * t205;
t256 = pkin(2) * t206;
t201 = cos(pkin(9));
t255 = pkin(4) * t201;
t253 = Icges(2,4) * t200;
t248 = t192 * t200;
t247 = t192 * t202;
t246 = t194 * t200;
t245 = t194 * t202;
t199 = sin(pkin(9));
t244 = t199 * t200;
t243 = t199 * t202;
t242 = t200 * t201;
t241 = t201 * t202;
t123 = -qJ(3) * t202 + t200 * t256;
t179 = pkin(1) * t200 - pkin(5) * t202;
t239 = -t123 - t179;
t124 = qJ(3) * t200 + t202 * t256;
t225 = pkin(3) * t194 + qJ(4) * t192;
t150 = t225 * t202;
t238 = -t124 - t150;
t237 = qJD(4) * t192;
t236 = qJD(5) * t192;
t235 = V_base(5) * qJ(1) + V_base(1);
t231 = qJD(1) + V_base(3);
t149 = t225 * t200;
t230 = -t149 + t239;
t160 = pkin(3) * t192 - qJ(4) * t194;
t229 = -t160 - t257;
t228 = qJD(3) * t200 + t185 * t257 + t235;
t227 = rSges(3,1) * t206 - rSges(3,2) * t205;
t226 = rSges(4,1) * t194 - rSges(4,2) * t192;
t180 = pkin(1) * t202 + pkin(5) * t200;
t218 = -V_base(4) * qJ(1) + V_base(6) * t180 + V_base(2);
t217 = t185 * t160 + t202 * t237 + t228;
t216 = V_base(4) * t179 - t180 * V_base(5) + t231;
t215 = pkin(6) * t192 + t194 * t255;
t214 = t186 * t123 + t216;
t211 = -qJD(3) * t202 + V_base(6) * t124 + t218;
t210 = V_base(6) * t150 + t200 * t237 + t211;
t209 = -qJD(4) * t194 + t186 * t149 + t214;
t197 = pkin(9) + qJ(5);
t195 = Icges(2,4) * t202;
t193 = cos(t197);
t191 = sin(t197);
t184 = t205 * rSges(3,1) + rSges(3,2) * t206;
t176 = rSges(2,1) * t202 - rSges(2,2) * t200;
t175 = rSges(2,1) * t200 + rSges(2,2) * t202;
t174 = -qJD(5) * t194 + V_base(6);
t173 = Icges(2,1) * t202 - t253;
t172 = Icges(2,1) * t200 + t195;
t171 = -Icges(2,2) * t200 + t195;
t170 = Icges(2,2) * t202 + t253;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = rSges(4,1) * t192 + rSges(4,2) * t194;
t156 = t202 * t236 + t186;
t155 = t200 * t236 + t185;
t154 = t194 * t241 + t244;
t153 = -t194 * t243 + t242;
t152 = t194 * t242 - t243;
t151 = -t194 * t244 - t241;
t148 = t200 * rSges(3,3) + t202 * t227;
t147 = -t202 * rSges(3,3) + t200 * t227;
t140 = t191 * t200 + t193 * t245;
t139 = -t191 * t245 + t193 * t200;
t138 = -t191 * t202 + t193 * t246;
t137 = -t191 * t246 - t193 * t202;
t134 = rSges(4,3) * t200 + t202 * t226;
t133 = -rSges(4,3) * t202 + t200 * t226;
t132 = V_base(5) * rSges(2,3) - t175 * V_base(6) + t235;
t131 = t176 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t122 = -rSges(5,3) * t194 + (rSges(5,1) * t201 - rSges(5,2) * t199) * t192;
t121 = -Icges(5,5) * t194 + (Icges(5,1) * t201 - Icges(5,4) * t199) * t192;
t120 = -Icges(5,6) * t194 + (Icges(5,4) * t201 - Icges(5,2) * t199) * t192;
t119 = -Icges(5,3) * t194 + (Icges(5,5) * t201 - Icges(5,6) * t199) * t192;
t116 = -rSges(6,3) * t194 + (rSges(6,1) * t193 - rSges(6,2) * t191) * t192;
t115 = -Icges(6,5) * t194 + (Icges(6,1) * t193 - Icges(6,4) * t191) * t192;
t114 = -Icges(6,6) * t194 + (Icges(6,4) * t193 - Icges(6,2) * t191) * t192;
t113 = -Icges(6,3) * t194 + (Icges(6,5) * t193 - Icges(6,6) * t191) * t192;
t112 = t175 * V_base(4) - t176 * V_base(5) + t231;
t111 = -pkin(6) * t194 + t192 * t255;
t109 = rSges(5,1) * t154 + rSges(5,2) * t153 + rSges(5,3) * t247;
t108 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t248;
t107 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t247;
t106 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t248;
t105 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t247;
t104 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t248;
t103 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t247;
t102 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t248;
t101 = pkin(4) * t244 + t202 * t215;
t100 = -pkin(4) * t243 + t200 * t215;
t99 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t247;
t98 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t248;
t97 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t247;
t96 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t248;
t95 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t247;
t94 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t248;
t93 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t247;
t92 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t248;
t91 = t184 * t185 + (-t147 - t179) * V_base(6) + t235;
t90 = t148 * V_base(6) - t184 * t186 + t218;
t89 = t147 * t186 - t148 * t185 + t216;
t88 = t161 * t185 + (-t133 + t239) * V_base(6) + t228;
t87 = t134 * V_base(6) + (-t161 - t257) * t186 + t211;
t86 = t133 * t186 + (-t124 - t134) * t185 + t214;
t85 = t122 * t185 + (-t108 + t230) * V_base(6) + t217;
t84 = t109 * V_base(6) + (-t122 + t229) * t186 + t210;
t83 = t108 * t186 + (-t109 + t238) * t185 + t209;
t82 = t111 * t185 + t116 * t155 - t174 * t98 + (-t100 + t230) * V_base(6) + t217;
t81 = t101 * V_base(6) - t116 * t156 + t174 * t99 + (-t111 + t229) * t186 + t210;
t80 = t100 * t186 - t155 * t99 + t156 * t98 + (-t101 + t238) * t185 + t209;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t112 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t156 * ((t139 * t95 + t140 * t97 + t93 * t247) * t156 + (t139 * t94 + t140 * t96 + t247 * t92) * t155 + (t113 * t247 + t114 * t139 + t115 * t140) * t174) / 0.2e1 + t155 * ((t137 * t95 + t138 * t97 + t248 * t93) * t156 + (t137 * t94 + t138 * t96 + t92 * t248) * t155 + (t113 * t248 + t114 * t137 + t115 * t138) * t174) / 0.2e1 + t174 * ((-t113 * t174 - t92 * t155 - t93 * t156) * t194 + ((-t191 * t95 + t193 * t97) * t156 + (-t191 * t94 + t193 * t96) * t155 + (-t114 * t191 + t115 * t193) * t174) * t192) / 0.2e1 + ((-t170 * t200 + t172 * t202 + Icges(1,4)) * V_base(5) + (-t171 * t200 + t173 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t170 * t202 + t172 * t200 + Icges(1,2)) * V_base(5) + (t171 * t202 + t173 * t200 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t103 * t248 + t105 * t151 + t107 * t152) * t186 + (t102 * t248 + t104 * t151 + t106 * t152) * t185 + (t119 * t248 + t120 * t151 + t121 * t152) * V_base(6) - t260 * t202 + t261 * t200) * t185 / 0.2e1 + ((t103 * t247 + t105 * t153 + t107 * t154) * t186 + (t102 * t247 + t104 * t153 + t106 * t154) * t185 + (t119 * t247 + t120 * t153 + t121 * t154) * V_base(6) + t261 * t202 + t260 * t200) * t186 / 0.2e1 + ((t144 * t206 + t205 * t146 + (t128 - t103) * t194 + (-t105 * t199 + t107 * t201 + t130) * t192) * t186 + (t143 * t206 + t205 * t145 + (t127 - t102) * t194 + (-t104 * t199 + t106 * t201 + t129) * t192) * t185 + (t182 * t206 + t205 * t183 + Icges(1,3) + Icges(2,3) + (t158 - t119) * t194 + (-t120 * t199 + t121 * t201 + t159) * t192) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t202 - Icges(2,6) * t200 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t200 + Icges(2,6) * t202 + Icges(1,6));
T = t1;
