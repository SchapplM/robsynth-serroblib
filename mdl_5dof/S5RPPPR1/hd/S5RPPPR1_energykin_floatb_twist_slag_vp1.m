% Calculate kinetic energy for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:09
% EndTime: 2022-01-20 09:12:11
% DurationCPUTime: 2.03s
% Computational Cost: add. (1318->297), mult. (1336->398), div. (0->0), fcn. (1228->10), ass. (0->141)
t191 = qJ(1) + pkin(7);
t183 = sin(t191);
t185 = cos(t191);
t197 = sin(qJ(1));
t198 = cos(qJ(1));
t245 = Icges(2,5) * t197 + Icges(3,5) * t183 + Icges(2,6) * t198 + Icges(3,6) * t185;
t244 = Icges(2,5) * t198 + Icges(3,5) * t185 - Icges(2,6) * t197 - Icges(3,6) * t183;
t242 = pkin(1) * t197;
t241 = pkin(1) * t198;
t194 = cos(pkin(9));
t240 = pkin(4) * t194;
t239 = -pkin(5) - qJ(2);
t238 = Icges(2,4) * t197;
t237 = Icges(3,4) * t183;
t193 = sin(pkin(8));
t236 = Icges(4,4) * t193;
t195 = cos(pkin(8));
t235 = Icges(4,4) * t195;
t192 = sin(pkin(9));
t234 = t183 * t192;
t233 = t183 * t193;
t232 = t183 * t195;
t231 = t185 * t192;
t230 = t185 * t193;
t229 = t185 * t195;
t228 = t192 * t195;
t227 = t194 * t195;
t225 = qJD(4) * t193;
t224 = qJD(5) * t193;
t186 = V_base(6) + qJD(1);
t223 = t186 * t241 + V_base(2);
t222 = V_base(5) * pkin(5) + V_base(1);
t164 = pkin(3) * t193 - qJ(4) * t195;
t219 = -t164 + t239;
t151 = pkin(2) * t183 - qJ(3) * t185;
t218 = -t151 - t242;
t153 = pkin(2) * t185 + qJ(3) * t183;
t217 = -t153 - t241;
t216 = V_base(5) * qJ(2) + t222;
t215 = V_base(4) * t242 + qJD(2) + V_base(3);
t209 = pkin(3) * t195 + qJ(4) * t193;
t141 = t209 * t183;
t214 = -t141 + t218;
t142 = t209 * t185;
t213 = -t142 + t217;
t212 = qJD(3) * t183 + t216;
t211 = V_base(4) * t151 + t215;
t210 = rSges(4,1) * t195 - rSges(4,2) * t193;
t208 = Icges(4,1) * t195 - t236;
t207 = -Icges(4,2) * t193 + t235;
t206 = Icges(4,5) * t195 - Icges(4,6) * t193;
t205 = -qJD(3) * t185 + t186 * t153 + t223;
t204 = V_base(5) * t164 + t185 * t225 + t212;
t203 = pkin(6) * t193 + t195 * t240;
t202 = t186 * t142 + t183 * t225 + t205;
t201 = -qJD(4) * t195 + V_base(4) * t141 + t211;
t200 = (-Icges(4,3) * t185 + t183 * t206) * V_base(5) + (Icges(4,3) * t183 + t185 * t206) * V_base(4) + (Icges(4,5) * t193 + Icges(4,6) * t195) * t186;
t114 = -Icges(4,6) * t185 + t183 * t207;
t115 = Icges(4,6) * t183 + t185 * t207;
t116 = -Icges(4,5) * t185 + t183 * t208;
t117 = Icges(4,5) * t183 + t185 * t208;
t162 = Icges(4,2) * t195 + t236;
t163 = Icges(4,1) * t193 + t235;
t199 = (-t115 * t193 + t117 * t195) * V_base(4) + (-t114 * t193 + t116 * t195) * V_base(5) + (-t162 * t193 + t163 * t195) * t186;
t190 = pkin(9) + qJ(5);
t188 = Icges(2,4) * t198;
t184 = cos(t190);
t182 = sin(t190);
t180 = Icges(3,4) * t185;
t176 = rSges(2,1) * t198 - t197 * rSges(2,2);
t175 = t197 * rSges(2,1) + rSges(2,2) * t198;
t174 = Icges(2,1) * t198 - t238;
t173 = Icges(2,1) * t197 + t188;
t172 = -Icges(2,2) * t197 + t188;
t171 = Icges(2,2) * t198 + t238;
t168 = -qJD(5) * t195 + t186;
t165 = rSges(4,1) * t193 + rSges(4,2) * t195;
t160 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t159 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t158 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t157 = t185 * t224 + V_base(4);
t156 = t183 * t224 + V_base(5);
t154 = rSges(3,1) * t185 - rSges(3,2) * t183;
t152 = rSges(3,1) * t183 + rSges(3,2) * t185;
t150 = Icges(3,1) * t185 - t237;
t149 = Icges(3,1) * t183 + t180;
t148 = -Icges(3,2) * t183 + t180;
t147 = Icges(3,2) * t185 + t237;
t140 = t185 * t227 + t234;
t139 = t183 * t194 - t185 * t228;
t138 = t183 * t227 - t231;
t137 = -t183 * t228 - t185 * t194;
t136 = -rSges(5,3) * t195 + (rSges(5,1) * t194 - rSges(5,2) * t192) * t193;
t135 = -Icges(5,5) * t195 + (Icges(5,1) * t194 - Icges(5,4) * t192) * t193;
t134 = -Icges(5,6) * t195 + (Icges(5,4) * t194 - Icges(5,2) * t192) * t193;
t133 = -Icges(5,3) * t195 + (Icges(5,5) * t194 - Icges(5,6) * t192) * t193;
t131 = t182 * t183 + t184 * t229;
t130 = -t182 * t229 + t183 * t184;
t129 = -t182 * t185 + t184 * t232;
t128 = -t182 * t232 - t184 * t185;
t127 = -rSges(6,3) * t195 + (rSges(6,1) * t184 - rSges(6,2) * t182) * t193;
t126 = -Icges(6,5) * t195 + (Icges(6,1) * t184 - Icges(6,4) * t182) * t193;
t125 = -Icges(6,6) * t195 + (Icges(6,4) * t184 - Icges(6,2) * t182) * t193;
t124 = -Icges(6,3) * t195 + (Icges(6,5) * t184 - Icges(6,6) * t182) * t193;
t122 = V_base(5) * rSges(2,3) - t175 * t186 + t222;
t121 = t176 * t186 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t120 = -pkin(6) * t195 + t193 * t240;
t119 = rSges(4,3) * t183 + t185 * t210;
t118 = -rSges(4,3) * t185 + t183 * t210;
t111 = t175 * V_base(4) - t176 * V_base(5) + V_base(3);
t110 = V_base(5) * rSges(3,3) + (-t152 - t242) * t186 + t216;
t109 = t154 * t186 + (-rSges(3,3) + t239) * V_base(4) + t223;
t108 = V_base(4) * t152 + (-t154 - t241) * V_base(5) + t215;
t107 = rSges(5,1) * t140 + rSges(5,2) * t139 + rSges(5,3) * t230;
t106 = rSges(5,1) * t138 + rSges(5,2) * t137 + rSges(5,3) * t233;
t105 = pkin(4) * t234 + t185 * t203;
t104 = -pkin(4) * t231 + t183 * t203;
t103 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t230;
t102 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t233;
t101 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t230;
t100 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t233;
t99 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t230;
t98 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t233;
t97 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t230;
t96 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t233;
t95 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t230;
t94 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t233;
t93 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t230;
t92 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t233;
t91 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t230;
t90 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t233;
t89 = t165 * V_base(5) + (-t118 + t218) * t186 + t212;
t88 = t119 * t186 + (-t165 + t239) * V_base(4) + t205;
t87 = V_base(4) * t118 + (-t119 + t217) * V_base(5) + t211;
t86 = t136 * V_base(5) + (-t106 + t214) * t186 + t204;
t85 = t107 * t186 + (-t136 + t219) * V_base(4) + t202;
t84 = V_base(4) * t106 + (-t107 + t213) * V_base(5) + t201;
t83 = t120 * V_base(5) + t127 * t156 - t168 * t96 + (-t104 + t214) * t186 + t204;
t82 = t105 * t186 - t127 * t157 + t168 * t97 + (-t120 + t219) * V_base(4) + t202;
t81 = V_base(4) * t104 - t156 * t97 + t157 * t96 + (-t105 + t213) * V_base(5) + t201;
t1 = m(1) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t157 * ((t130 * t93 + t131 * t95 + t91 * t230) * t157 + (t130 * t92 + t131 * t94 + t230 * t90) * t156 + (t124 * t230 + t125 * t130 + t126 * t131) * t168) / 0.2e1 + t156 * ((t128 * t93 + t129 * t95 + t233 * t91) * t157 + (t128 * t92 + t129 * t94 + t90 * t233) * t156 + (t124 * t233 + t125 * t128 + t126 * t129) * t168) / 0.2e1 + t168 * ((-t124 * t168 - t90 * t156 - t91 * t157) * t195 + ((-t182 * t93 + t184 * t95) * t157 + (-t182 * t92 + t184 * t94) * t156 + (-t125 * t182 + t126 * t184) * t168) * t193) / 0.2e1 + (((t114 - t98) * t195 + (-t100 * t192 + t102 * t194 + t116) * t193 + t245) * V_base(5) + ((t115 - t99) * t195 + (-t101 * t192 + t103 * t194 + t117) * t193 + t244) * V_base(4) + (Icges(2,3) + Icges(3,3) + (t162 - t133) * t195 + (-t134 * t192 + t135 * t194 + t163) * t193) * t186) * t186 / 0.2e1 + (t200 * t183 + t199 * t185 + (t133 * t230 + t134 * t139 + t135 * t140 + t244) * t186 + (t100 * t139 + t102 * t140 - t147 * t183 + t149 * t185 - t197 * t171 + t173 * t198 + t230 * t98 + Icges(1,4)) * V_base(5) + (t139 * t101 + t140 * t103 - t183 * t148 + t185 * t150 - t197 * t172 + t198 * t174 + t230 * t99 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t199 * t183 - t200 * t185 + (t133 * t233 + t134 * t137 + t135 * t138 + t245) * t186 + (t137 * t100 + t138 * t102 + t185 * t147 + t183 * t149 + t198 * t171 + t197 * t173 + t233 * t98 + Icges(1,2)) * V_base(5) + (t101 * t137 + t103 * t138 + t148 * t185 + t150 * t183 + t172 * t198 + t197 * t174 + t233 * t99 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
