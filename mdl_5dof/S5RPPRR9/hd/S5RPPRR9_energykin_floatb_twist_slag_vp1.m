% Calculate kinetic energy for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:15
% EndTime: 2019-12-31 18:02:18
% DurationCPUTime: 2.18s
% Computational Cost: add. (889->242), mult. (1630->333), div. (0->0), fcn. (1794->8), ass. (0->118)
t226 = Icges(2,4) - Icges(3,5);
t225 = Icges(2,1) + Icges(3,1);
t224 = Icges(3,4) + Icges(2,5);
t223 = Icges(2,2) + Icges(3,3);
t222 = Icges(2,6) - Icges(3,6);
t211 = sin(qJ(1));
t221 = t226 * t211;
t212 = cos(qJ(1));
t220 = t226 * t212;
t219 = -t223 * t212 - t221;
t218 = t223 * t211 - t220;
t217 = t225 * t211 + t220;
t216 = t225 * t212 - t221;
t210 = cos(pkin(8));
t209 = sin(pkin(8));
t131 = -t209 * t211 - t210 * t212;
t208 = Icges(4,4) * t131;
t173 = sin(qJ(4));
t207 = Icges(5,4) * t173;
t175 = cos(qJ(4));
t206 = Icges(5,4) * t175;
t205 = t131 * t173;
t132 = t209 * t212 - t210 * t211;
t204 = t132 * t173;
t172 = sin(qJ(5));
t203 = t172 * t175;
t174 = cos(qJ(5));
t202 = t174 * t175;
t201 = qJD(5) * t173;
t154 = pkin(1) * t211 - qJ(2) * t212;
t200 = V_base(4) * t154 + V_base(3);
t199 = V_base(5) * pkin(5) + V_base(1);
t196 = t212 * pkin(2);
t195 = t211 * pkin(2);
t126 = qJD(4) * t132 + V_base(4);
t125 = -qJD(4) * t131 + V_base(5);
t166 = V_base(6) + qJD(1);
t192 = qJD(2) * t211 + t199;
t191 = -t154 - t195;
t157 = pkin(1) * t212 + qJ(2) * t211;
t190 = -t157 - t196;
t189 = -pkin(4) * t175 - pkin(7) * t173;
t188 = V_base(4) * t195 - qJD(3) + t200;
t187 = -rSges(5,1) * t175 + rSges(5,2) * t173;
t186 = -Icges(5,1) * t175 + t207;
t185 = Icges(5,2) * t173 - t206;
t184 = -Icges(5,5) * t175 + Icges(5,6) * t173;
t116 = -pkin(3) * t132 - pkin(6) * t131;
t183 = -t116 + t191;
t182 = -qJD(2) * t212 + t166 * t157 + V_base(2);
t181 = -V_base(5) * qJ(3) + t192;
t180 = t125 * (-Icges(5,3) * t131 + t132 * t184) + t126 * (Icges(5,3) * t132 + t131 * t184) + (-Icges(5,5) * t173 - Icges(5,6) * t175) * t166;
t179 = V_base(4) * qJ(3) + t166 * t196 + t182;
t117 = -pkin(3) * t131 + pkin(6) * t132;
t178 = -V_base(4) * pkin(5) + t166 * t117 + t179;
t177 = V_base(4) * t116 + (-t117 + t190) * V_base(5) + t188;
t142 = -Icges(5,2) * t175 - t207;
t147 = -Icges(5,1) * t173 - t206;
t90 = -Icges(5,6) * t131 + t132 * t185;
t91 = Icges(5,6) * t132 + t131 * t185;
t92 = -Icges(5,5) * t131 + t132 * t186;
t93 = Icges(5,5) * t132 + t131 * t186;
t176 = (t173 * t91 - t175 * t93) * t126 + (t173 * t90 - t175 * t92) * t125 + (t142 * t173 - t147 * t175) * t166;
t160 = -t173 * pkin(4) + pkin(7) * t175;
t159 = rSges(2,1) * t212 - rSges(2,2) * t211;
t158 = rSges(3,1) * t212 + rSges(3,3) * t211;
t156 = rSges(2,1) * t211 + rSges(2,2) * t212;
t155 = rSges(3,1) * t211 - rSges(3,3) * t212;
t153 = -t173 * rSges(5,1) - rSges(5,2) * t175;
t152 = qJD(5) * t175 + t166;
t136 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t135 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t134 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = Icges(4,4) * t132;
t124 = rSges(6,3) * t175 + (-rSges(6,1) * t174 + rSges(6,2) * t172) * t173;
t123 = Icges(6,5) * t175 + (-Icges(6,1) * t174 + Icges(6,4) * t172) * t173;
t122 = Icges(6,6) * t175 + (-Icges(6,4) * t174 + Icges(6,2) * t172) * t173;
t121 = Icges(6,3) * t175 + (-Icges(6,5) * t174 + Icges(6,6) * t172) * t173;
t120 = V_base(5) * rSges(2,3) - t156 * t166 + t199;
t119 = t159 * t166 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t118 = t156 * V_base(4) - t159 * V_base(5) + V_base(3);
t115 = -rSges(4,1) * t131 - rSges(4,2) * t132;
t114 = -rSges(4,1) * t132 + rSges(4,2) * t131;
t113 = -Icges(4,1) * t131 - t129;
t112 = -Icges(4,1) * t132 + t208;
t111 = -Icges(4,2) * t132 - t208;
t110 = Icges(4,2) * t131 - t129;
t106 = t189 * t131;
t105 = t189 * t132;
t103 = -t131 * t202 + t132 * t172;
t102 = t131 * t203 + t132 * t174;
t101 = -t131 * t172 - t132 * t202;
t100 = -t131 * t174 + t132 * t203;
t99 = -t131 * t201 + t126;
t98 = -t132 * t201 + t125;
t97 = V_base(5) * rSges(3,2) + (-t154 - t155) * t166 + t192;
t96 = t166 * t158 + (-pkin(5) - rSges(3,2)) * V_base(4) + t182;
t95 = t132 * rSges(5,3) + t131 * t187;
t94 = -t131 * rSges(5,3) + t132 * t187;
t87 = t155 * V_base(4) + (-t157 - t158) * V_base(5) + t200;
t86 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t114 + t191) * t166 + t192;
t85 = t166 * t115 + (rSges(4,3) - pkin(5)) * V_base(4) + t179;
t84 = rSges(6,1) * t103 + rSges(6,2) * t102 - rSges(6,3) * t205;
t83 = rSges(6,1) * t101 + rSges(6,2) * t100 - rSges(6,3) * t204;
t82 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t205;
t81 = Icges(6,1) * t101 + Icges(6,4) * t100 - Icges(6,5) * t204;
t80 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t205;
t79 = Icges(6,4) * t101 + Icges(6,2) * t100 - Icges(6,6) * t204;
t78 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t205;
t77 = Icges(6,5) * t101 + Icges(6,6) * t100 - Icges(6,3) * t204;
t76 = V_base(4) * t114 + (-t115 + t190) * V_base(5) + t188;
t75 = t125 * t153 + (-t94 + t183) * t166 + t181;
t74 = -t126 * t153 + t166 * t95 + t178;
t73 = -t125 * t95 + t126 * t94 + t177;
t72 = t98 * t124 + t125 * t160 - t152 * t83 + (-t105 + t183) * t166 + t181;
t71 = t166 * t106 - t99 * t124 - t126 * t160 + t152 * t84 + t178;
t70 = t126 * t105 - t125 * t106 + t99 * t83 - t98 * t84 + t177;
t1 = m(1) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(2) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t126 * (t176 * t131 + t180 * t132) / 0.2e1 + t125 * (-t180 * t131 + t176 * t132) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t99 * ((t102 * t80 + t103 * t82 - t205 * t78) * t99 + (t102 * t79 + t103 * t81 - t205 * t77) * t98 + (t102 * t122 + t103 * t123 - t121 * t205) * t152) / 0.2e1 + t98 * ((t100 * t80 + t101 * t82 - t204 * t78) * t99 + (t100 * t79 + t101 * t81 - t204 * t77) * t98 + (t100 * t122 + t101 * t123 - t121 * t204) * t152) / 0.2e1 + t152 * ((t121 * t152 + t77 * t98 + t78 * t99) * t175 + ((t172 * t80 - t174 * t82) * t99 + (t172 * t79 - t174 * t81) * t98 + (t122 * t172 - t123 * t174) * t152) * t173) / 0.2e1 + ((-t173 * t93 - t175 * t91) * t126 + (-t173 * t92 - t175 * t90) * t125 + (-t142 * t175 - t173 * t147 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t166) * t166 / 0.2e1 + ((-t110 * t132 - t112 * t131 + t211 * t219 + t217 * t212 + Icges(1,4)) * V_base(5) + (-t111 * t132 - t113 * t131 + t218 * t211 + t216 * t212 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t110 * t131 - t112 * t132 + t217 * t211 - t219 * t212 + Icges(1,2)) * V_base(5) + (t111 * t131 - t113 * t132 + t211 * t216 - t212 * t218 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t166 * (Icges(4,5) * t132 - Icges(4,6) * t131 + t224 * t211 + t222 * t212) + V_base(4) * t166 * (Icges(4,5) * t131 + Icges(4,6) * t132 - t222 * t211 + t224 * t212) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
