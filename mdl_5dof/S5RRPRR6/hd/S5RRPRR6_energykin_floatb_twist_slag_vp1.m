% Calculate kinetic energy for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:37
% EndTime: 2019-12-05 18:35:40
% DurationCPUTime: 2.60s
% Computational Cost: add. (1410->293), mult. (1396->424), div. (0->0), fcn. (1288->10), ass. (0->140)
t192 = qJ(1) + qJ(2);
t186 = sin(t192);
t188 = cos(t192);
t193 = sin(pkin(9));
t194 = cos(pkin(9));
t231 = Icges(4,4) * t194;
t211 = -Icges(4,2) * t193 + t231;
t115 = Icges(4,6) * t188 - t186 * t211;
t232 = Icges(4,4) * t193;
t212 = Icges(4,1) * t194 - t232;
t117 = Icges(4,5) * t188 - t186 * t212;
t233 = Icges(3,4) * t188;
t245 = -Icges(3,1) * t186 - t115 * t193 + t117 * t194 - t233;
t116 = Icges(4,6) * t186 + t188 * t211;
t118 = Icges(4,5) * t186 + t188 * t212;
t234 = Icges(3,4) * t186;
t244 = Icges(3,1) * t188 - t116 * t193 + t118 * t194 - t234;
t197 = cos(qJ(4));
t238 = pkin(4) * t197;
t243 = pkin(8) * t193 + t194 * t238;
t242 = -pkin(5) - pkin(6);
t196 = sin(qJ(1));
t240 = pkin(1) * t196;
t198 = cos(qJ(1));
t239 = pkin(1) * t198;
t236 = Icges(2,4) * t196;
t235 = Icges(2,4) * t198;
t230 = t186 * t193;
t229 = t186 * t194;
t195 = sin(qJ(4));
t228 = t186 * t195;
t227 = t188 * t193;
t226 = t188 * t194;
t225 = t188 * t195;
t224 = t194 * t195;
t223 = t194 * t197;
t222 = qJD(4) * t193;
t221 = -qJD(4) - qJD(5);
t220 = V_base(6) * pkin(5) + V_base(2);
t158 = t188 * t222 + V_base(6);
t184 = V_base(4) + qJD(1);
t217 = V_base(5) * t239 + V_base(6) * t240 + V_base(1);
t182 = qJD(2) + t184;
t155 = pkin(2) * t188 + qJ(3) * t186;
t216 = V_base(5) * t155 + t217;
t215 = pkin(3) * t194 + pkin(7) * t193;
t214 = rSges(4,1) * t194 - rSges(4,2) * t193;
t213 = -t184 * t240 + V_base(3);
t210 = Icges(4,5) * t194 - Icges(4,6) * t193;
t165 = Icges(4,2) * t194 + t232;
t166 = Icges(4,1) * t193 + t231;
t207 = t165 * t193 - t166 * t194;
t153 = -pkin(2) * t186 + qJ(3) * t188;
t206 = qJD(3) * t186 + t182 * t153 + t213;
t205 = V_base(6) * pkin(6) - t184 * t239 + t220;
t204 = qJD(3) * t188 + t205;
t203 = (Icges(4,3) * t188 - t186 * t210) * V_base(5) + (Icges(4,3) * t186 + t188 * t210) * V_base(6) + (Icges(4,5) * t193 + Icges(4,6) * t194) * t182;
t143 = t215 * t186;
t144 = t215 * t188;
t202 = V_base(5) * t144 + (t143 - t153) * V_base(6) + t216;
t168 = t193 * pkin(3) - t194 * pkin(7);
t201 = V_base(6) * t168 + (-t144 - t155) * t182 + t204;
t200 = -t182 * t143 + (-t168 + t242) * V_base(5) + t206;
t191 = qJ(4) + qJ(5);
t187 = cos(t191);
t185 = sin(t191);
t177 = rSges(2,1) * t198 - rSges(2,2) * t196;
t176 = -rSges(2,1) * t196 - rSges(2,2) * t198;
t174 = Icges(2,1) * t198 - t236;
t173 = -Icges(2,1) * t196 - t235;
t172 = -Icges(2,2) * t196 + t235;
t171 = -Icges(2,2) * t198 - t236;
t167 = rSges(4,1) * t193 + rSges(4,2) * t194;
t163 = -qJD(4) * t194 + t182;
t162 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t161 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t160 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t159 = -t186 * t222 + V_base(5);
t156 = rSges(3,1) * t188 - rSges(3,2) * t186;
t154 = -rSges(3,1) * t186 - rSges(3,2) * t188;
t150 = -Icges(3,2) * t186 + t233;
t149 = -Icges(3,2) * t188 - t234;
t148 = Icges(3,5) * t188 - Icges(3,6) * t186;
t147 = -Icges(3,5) * t186 - Icges(3,6) * t188;
t146 = t194 * t221 + t182;
t142 = t188 * t223 + t228;
t141 = t186 * t197 - t188 * t224;
t140 = -t186 * t223 + t225;
t139 = t186 * t224 + t188 * t197;
t137 = -rSges(5,3) * t194 + (rSges(5,1) * t197 - rSges(5,2) * t195) * t193;
t136 = -Icges(5,5) * t194 + (Icges(5,1) * t197 - Icges(5,4) * t195) * t193;
t135 = -Icges(5,6) * t194 + (Icges(5,4) * t197 - Icges(5,2) * t195) * t193;
t134 = -Icges(5,3) * t194 + (Icges(5,5) * t197 - Icges(5,6) * t195) * t193;
t132 = t221 * t230 + V_base(5);
t131 = qJD(5) * t227 + t158;
t130 = t185 * t186 + t187 * t226;
t129 = -t185 * t226 + t186 * t187;
t128 = t185 * t188 - t187 * t229;
t127 = t185 * t229 + t187 * t188;
t126 = -rSges(6,3) * t194 + (rSges(6,1) * t187 - rSges(6,2) * t185) * t193;
t125 = -Icges(6,5) * t194 + (Icges(6,1) * t187 - Icges(6,4) * t185) * t193;
t124 = -Icges(6,6) * t194 + (Icges(6,4) * t187 - Icges(6,2) * t185) * t193;
t123 = -Icges(6,3) * t194 + (Icges(6,5) * t187 - Icges(6,6) * t185) * t193;
t122 = -pkin(8) * t194 + t193 * t238;
t121 = rSges(4,3) * t186 + t188 * t214;
t120 = rSges(4,3) * t188 - t186 * t214;
t112 = V_base(6) * rSges(2,3) - t177 * t184 + t220;
t111 = t176 * t184 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t110 = -t176 * V_base(6) + t177 * V_base(5) + V_base(1);
t109 = V_base(6) * rSges(3,3) - t156 * t182 + t205;
t108 = t154 * t182 + (-rSges(3,3) + t242) * V_base(5) + t213;
t107 = -t154 * V_base(6) + t156 * V_base(5) + t217;
t106 = rSges(5,1) * t142 + rSges(5,2) * t141 + rSges(5,3) * t227;
t105 = rSges(5,1) * t140 + rSges(5,2) * t139 - rSges(5,3) * t230;
t104 = pkin(4) * t228 + t188 * t243;
t103 = pkin(4) * t225 - t186 * t243;
t102 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t227;
t101 = Icges(5,1) * t140 + Icges(5,4) * t139 - Icges(5,5) * t230;
t100 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t227;
t99 = Icges(5,4) * t140 + Icges(5,2) * t139 - Icges(5,6) * t230;
t98 = Icges(5,5) * t142 + Icges(5,6) * t141 + Icges(5,3) * t227;
t97 = Icges(5,5) * t140 + Icges(5,6) * t139 - Icges(5,3) * t230;
t96 = rSges(6,1) * t130 + rSges(6,2) * t129 + rSges(6,3) * t227;
t95 = rSges(6,1) * t128 + rSges(6,2) * t127 - rSges(6,3) * t230;
t94 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t227;
t93 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t230;
t92 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t227;
t91 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t230;
t90 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t227;
t89 = Icges(6,5) * t128 + Icges(6,6) * t127 - Icges(6,3) * t230;
t88 = t167 * V_base(6) + (-t121 - t155) * t182 + t204;
t87 = t120 * t182 + (-t167 + t242) * V_base(5) + t206;
t86 = t121 * V_base(5) + (-t120 - t153) * V_base(6) + t216;
t85 = -t106 * t163 + t137 * t158 + t201;
t84 = t105 * t163 - t137 * t159 + t200;
t83 = -t105 * t158 + t106 * t159 + t202;
t82 = -t104 * t163 + t122 * t158 + t126 * t131 - t146 * t96 + t201;
t81 = t103 * t163 - t122 * t159 - t126 * t132 + t146 * t95 + t200;
t80 = -t103 * t158 + t104 * t159 - t131 * t95 + t132 * t96 + t202;
t1 = m(1) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(2) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(3) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t163 * ((-t134 * t163 - t98 * t158 - t97 * t159) * t194 + ((-t135 * t195 + t136 * t197) * t163 + (t101 * t197 - t195 * t99) * t159 + (-t100 * t195 + t102 * t197) * t158) * t193) / 0.2e1 + t159 * ((-t134 * t230 + t135 * t139 + t136 * t140) * t163 + (t101 * t140 + t139 * t99 - t97 * t230) * t159 + (t100 * t139 + t102 * t140 - t230 * t98) * t158) / 0.2e1 + t158 * ((t134 * t227 + t135 * t141 + t136 * t142) * t163 + (t101 * t142 + t141 * t99 + t227 * t97) * t159 + (t100 * t141 + t102 * t142 + t98 * t227) * t158) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t146 * ((-t123 * t146 - t90 * t131 - t89 * t132) * t194 + ((-t124 * t185 + t125 * t187) * t146 + (-t185 * t91 + t187 * t93) * t132 + (-t185 * t92 + t187 * t94) * t131) * t193) / 0.2e1 + t132 * ((-t123 * t230 + t124 * t127 + t125 * t128) * t146 + (t127 * t91 + t128 * t93 - t89 * t230) * t132 + (t127 * t92 + t128 * t94 - t230 * t90) * t131) / 0.2e1 + t131 * ((t123 * t227 + t124 * t129 + t125 * t130) * t146 + (t129 * t91 + t130 * t93 + t227 * t89) * t132 + (t129 * t92 + t130 * t94 + t90 * t227) * t131) / 0.2e1 + ((t116 * t194 + t118 * t193 + t148) * V_base(6) + (t115 * t194 + t117 * t193 + t147) * V_base(5) + (t165 * t194 + t166 * t193 + Icges(3,3)) * t182) * t182 / 0.2e1 + (t203 * t188 + (t207 * t186 + t147) * t182 + (-t150 * t188 - t172 * t198 - t174 * t196 - t186 * t244 + Icges(1,6)) * V_base(6) + (-t149 * t188 - t171 * t198 - t173 * t196 - t245 * t186 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t203 * t186 + (-t207 * t188 + t148) * t182 + (-t150 * t186 - t172 * t196 + t174 * t198 + t244 * t188 + Icges(1,3)) * V_base(6) + (-t149 * t186 - t171 * t196 + t173 * t198 + t245 * t188 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t196 - Icges(2,6) * t198) * V_base(5) + (Icges(2,5) * t198 - Icges(2,6) * t196) * V_base(6) + Icges(2,3) * t184 / 0.2e1) * t184;
T = t1;
