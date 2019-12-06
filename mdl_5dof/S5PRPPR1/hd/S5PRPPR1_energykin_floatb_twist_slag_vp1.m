% Calculate kinetic energy for
% S5PRPPR1
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:43
% EndTime: 2019-12-05 15:21:45
% DurationCPUTime: 2.28s
% Computational Cost: add. (1305->296), mult. (1336->401), div. (0->0), fcn. (1228->10), ass. (0->141)
t197 = cos(pkin(7));
t243 = pkin(1) * t197;
t195 = cos(pkin(9));
t242 = pkin(4) * t195;
t241 = -pkin(5) - qJ(1);
t194 = sin(pkin(7));
t240 = Icges(2,4) * t194;
t191 = pkin(7) + qJ(2);
t183 = sin(t191);
t239 = Icges(3,4) * t183;
t193 = sin(pkin(8));
t238 = Icges(4,4) * t193;
t196 = cos(pkin(8));
t237 = Icges(4,4) * t196;
t192 = sin(pkin(9));
t236 = t183 * t192;
t235 = t183 * t193;
t234 = t183 * t196;
t185 = cos(t191);
t233 = t185 * t192;
t232 = t185 * t193;
t231 = t185 * t196;
t230 = t192 * t196;
t229 = t195 * t196;
t211 = pkin(3) * t196 + qJ(4) * t193;
t141 = t211 * t183;
t151 = pkin(2) * t183 - qJ(3) * t185;
t227 = -t141 - t151;
t226 = qJD(4) * t193;
t225 = qJD(5) * t193;
t217 = pkin(1) * V_base(6);
t224 = t197 * t217 + V_base(2);
t223 = V_base(5) * qJ(1) + V_base(1);
t219 = qJD(1) + V_base(3);
t170 = t193 * pkin(3) - t196 * qJ(4);
t218 = -t170 + t241;
t187 = V_base(6) + qJD(2);
t153 = pkin(2) * t185 + qJ(3) * t183;
t216 = -t153 - t243;
t215 = V_base(4) * t194 * pkin(1) + t219;
t142 = t211 * t185;
t214 = -t142 + t216;
t213 = V_base(4) * t151 + t215;
t212 = rSges(4,1) * t196 - rSges(4,2) * t193;
t210 = Icges(4,1) * t196 - t238;
t209 = -Icges(4,2) * t193 + t237;
t208 = Icges(4,5) * t196 - Icges(4,6) * t193;
t207 = -qJD(3) * t185 + t187 * t153 + t224;
t206 = pkin(6) * t193 + t196 * t242;
t205 = V_base(5) * pkin(5) - t194 * t217 + t223;
t204 = t187 * t142 + t183 * t226 + t207;
t203 = -qJD(4) * t196 + V_base(4) * t141 + t213;
t202 = qJD(3) * t183 + t205;
t201 = (-Icges(4,3) * t185 + t183 * t208) * V_base(5) + (Icges(4,3) * t183 + t185 * t208) * V_base(4) + (Icges(4,5) * t193 + Icges(4,6) * t196) * t187;
t200 = V_base(5) * t170 + t185 * t226 + t202;
t114 = -Icges(4,6) * t185 + t183 * t209;
t115 = Icges(4,6) * t183 + t185 * t209;
t116 = -Icges(4,5) * t185 + t183 * t210;
t117 = Icges(4,5) * t183 + t185 * t210;
t164 = Icges(4,2) * t196 + t238;
t167 = Icges(4,1) * t193 + t237;
t199 = (-t115 * t193 + t117 * t196) * V_base(4) + (-t114 * t193 + t116 * t196) * V_base(5) + (-t164 * t193 + t167 * t196) * t187;
t190 = pkin(9) + qJ(5);
t186 = Icges(2,4) * t197;
t184 = cos(t190);
t182 = sin(t190);
t180 = Icges(3,4) * t185;
t176 = -qJD(5) * t196 + t187;
t173 = rSges(2,1) * t197 - rSges(2,2) * t194;
t172 = rSges(2,1) * t194 + rSges(2,2) * t197;
t171 = rSges(4,1) * t193 + rSges(4,2) * t196;
t169 = Icges(2,1) * t197 - t240;
t168 = Icges(2,1) * t194 + t186;
t166 = -Icges(2,2) * t194 + t186;
t165 = Icges(2,2) * t197 + t240;
t160 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t159 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t158 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t157 = t185 * t225 + V_base(4);
t156 = t183 * t225 + V_base(5);
t154 = rSges(3,1) * t185 - rSges(3,2) * t183;
t152 = rSges(3,1) * t183 + rSges(3,2) * t185;
t150 = Icges(3,1) * t185 - t239;
t149 = Icges(3,1) * t183 + t180;
t148 = -Icges(3,2) * t183 + t180;
t147 = Icges(3,2) * t185 + t239;
t146 = Icges(3,5) * t185 - Icges(3,6) * t183;
t145 = Icges(3,5) * t183 + Icges(3,6) * t185;
t140 = t185 * t229 + t236;
t139 = t183 * t195 - t185 * t230;
t138 = t183 * t229 - t233;
t137 = -t183 * t230 - t185 * t195;
t136 = -rSges(5,3) * t196 + (rSges(5,1) * t195 - rSges(5,2) * t192) * t193;
t135 = -Icges(5,5) * t196 + (Icges(5,1) * t195 - Icges(5,4) * t192) * t193;
t134 = -Icges(5,6) * t196 + (Icges(5,4) * t195 - Icges(5,2) * t192) * t193;
t133 = -Icges(5,3) * t196 + (Icges(5,5) * t195 - Icges(5,6) * t192) * t193;
t131 = t182 * t183 + t184 * t231;
t130 = -t182 * t231 + t183 * t184;
t129 = -t182 * t185 + t184 * t234;
t128 = -t182 * t234 - t184 * t185;
t127 = -rSges(6,3) * t196 + (rSges(6,1) * t184 - rSges(6,2) * t182) * t193;
t126 = V_base(5) * rSges(2,3) - t172 * V_base(6) + t223;
t125 = t173 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t124 = -Icges(6,5) * t196 + (Icges(6,1) * t184 - Icges(6,4) * t182) * t193;
t123 = -Icges(6,6) * t196 + (Icges(6,4) * t184 - Icges(6,2) * t182) * t193;
t122 = -Icges(6,3) * t196 + (Icges(6,5) * t184 - Icges(6,6) * t182) * t193;
t120 = -pkin(6) * t196 + t193 * t242;
t119 = rSges(4,3) * t183 + t185 * t212;
t118 = -rSges(4,3) * t185 + t183 * t212;
t111 = t172 * V_base(4) - t173 * V_base(5) + t219;
t110 = V_base(5) * rSges(3,3) - t152 * t187 + t205;
t109 = t154 * t187 + (-rSges(3,3) + t241) * V_base(4) + t224;
t108 = t152 * V_base(4) + (-t154 - t243) * V_base(5) + t215;
t107 = rSges(5,1) * t140 + rSges(5,2) * t139 + rSges(5,3) * t232;
t106 = rSges(5,1) * t138 + rSges(5,2) * t137 + rSges(5,3) * t235;
t105 = pkin(4) * t236 + t185 * t206;
t104 = -pkin(4) * t233 + t183 * t206;
t103 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t232;
t102 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t235;
t101 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t232;
t100 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t235;
t99 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t232;
t98 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t235;
t97 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t232;
t96 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t235;
t95 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t232;
t94 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t235;
t93 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t232;
t92 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t235;
t91 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t232;
t90 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t235;
t89 = t171 * V_base(5) + (-t118 - t151) * t187 + t202;
t88 = t119 * t187 + (-t171 + t241) * V_base(4) + t207;
t87 = t118 * V_base(4) + (-t119 + t216) * V_base(5) + t213;
t86 = t136 * V_base(5) + (-t106 + t227) * t187 + t200;
t85 = t107 * t187 + (-t136 + t218) * V_base(4) + t204;
t84 = t106 * V_base(4) + (-t107 + t214) * V_base(5) + t203;
t83 = t120 * V_base(5) + t127 * t156 - t176 * t96 + (-t104 + t227) * t187 + t200;
t82 = t105 * t187 - t127 * t157 + t176 * t97 + (-t120 + t218) * V_base(4) + t204;
t81 = t104 * V_base(4) - t156 * t97 + t157 * t96 + (-t105 + t214) * V_base(5) + t203;
t1 = m(1) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t157 * ((t130 * t93 + t131 * t95 + t232 * t91) * t157 + (t130 * t92 + t131 * t94 + t232 * t90) * t156 + (t122 * t232 + t123 * t130 + t124 * t131) * t176) / 0.2e1 + t156 * ((t128 * t93 + t129 * t95 + t235 * t91) * t157 + (t128 * t92 + t129 * t94 + t90 * t235) * t156 + (t122 * t235 + t123 * t128 + t124 * t129) * t176) / 0.2e1 + t176 * ((-t122 * t176 - t90 * t156 - t91 * t157) * t196 + ((-t182 * t93 + t184 * t95) * t157 + (-t182 * t92 + t184 * t94) * t156 + (-t123 * t182 + t124 * t184) * t176) * t193) / 0.2e1 + ((t145 + (t114 - t98) * t196 + (-t100 * t192 + t102 * t195 + t116) * t193) * V_base(5) + (t146 + (t115 - t99) * t196 + (-t101 * t192 + t103 * t195 + t117) * t193) * V_base(4) + (Icges(3,3) + (t164 - t133) * t196 + (-t134 * t192 + t135 * t195 + t167) * t193) * t187) * t187 / 0.2e1 + (t183 * t201 + t185 * t199 + (t133 * t232 + t134 * t139 + t135 * t140 + t146) * t187 + (t100 * t139 + t102 * t140 - t147 * t183 + t149 * t185 - t165 * t194 + t168 * t197 + t232 * t98 + Icges(1,4)) * V_base(5) + (t139 * t101 + t140 * t103 - t183 * t148 + t185 * t150 - t194 * t166 + t197 * t169 + t232 * t99 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t183 * t199 - t201 * t185 + (t133 * t235 + t134 * t137 + t135 * t138 + t145) * t187 + (t137 * t100 + t138 * t102 + t185 * t147 + t183 * t149 + t197 * t165 + t194 * t168 + t235 * t98 + Icges(1,2)) * V_base(5) + (t101 * t137 + t103 * t138 + t148 * t185 + t150 * t183 + t166 * t197 + t169 * t194 + t235 * t99 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t194 + Icges(2,6) * t197 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t197 - Icges(2,6) * t194 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
