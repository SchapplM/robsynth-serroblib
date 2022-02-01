% Calculate kinetic energy for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:16
% EndTime: 2022-01-20 09:51:17
% DurationCPUTime: 1.40s
% Computational Cost: add. (1114->232), mult. (736->299), div. (0->0), fcn. (514->10), ass. (0->121)
t156 = qJ(1) + qJ(2);
t147 = pkin(8) + t156;
t141 = sin(t147);
t142 = cos(t147);
t149 = sin(t156);
t150 = cos(t156);
t207 = Icges(3,5) * t149 + Icges(4,5) * t141 + Icges(3,6) * t150 + Icges(4,6) * t142;
t206 = Icges(3,5) * t150 + Icges(4,5) * t142 - Icges(3,6) * t149 - Icges(4,6) * t141;
t205 = -pkin(5) - pkin(6);
t160 = sin(qJ(1));
t203 = pkin(1) * t160;
t161 = cos(qJ(1));
t202 = pkin(1) * t161;
t201 = pkin(2) * t149;
t200 = pkin(2) * t150;
t157 = sin(pkin(9));
t199 = pkin(4) * t157;
t158 = cos(pkin(9));
t198 = pkin(4) * t158;
t197 = Icges(2,4) * t160;
t196 = Icges(3,4) * t149;
t195 = Icges(4,4) * t141;
t194 = Icges(5,4) * t157;
t193 = Icges(5,4) * t158;
t155 = pkin(9) + qJ(5);
t145 = sin(t155);
t192 = Icges(6,4) * t145;
t146 = cos(t155);
t191 = Icges(6,4) * t146;
t189 = -qJ(3) + t205;
t148 = V_base(6) + qJD(1);
t188 = t148 * t202 + V_base(2);
t187 = V_base(4) * t203 + V_base(3);
t186 = V_base(5) * pkin(5) + V_base(1);
t101 = pkin(3) * t141 - qJ(4) * t142;
t183 = -t101 - t201;
t143 = qJD(2) + t148;
t182 = t143 * t200 + t188;
t181 = -t200 - t202;
t180 = V_base(4) * t201 + qJD(3) + t187;
t179 = rSges(5,1) * t158 - rSges(5,2) * t157;
t178 = rSges(6,1) * t146 - rSges(6,2) * t145;
t177 = Icges(5,1) * t158 - t194;
t176 = Icges(6,1) * t146 - t192;
t175 = -Icges(5,2) * t157 + t193;
t174 = -Icges(6,2) * t145 + t191;
t173 = Icges(5,5) * t158 - Icges(5,6) * t157;
t172 = Icges(6,5) * t146 - Icges(6,6) * t145;
t171 = t101 * V_base(4) + t180;
t103 = pkin(3) * t142 + qJ(4) * t141;
t170 = -t103 + t181;
t169 = V_base(5) * pkin(6) - t148 * t203 + t186;
t168 = -qJD(4) * t142 + t103 * t143 + t182;
t118 = -qJD(5) * t142 + V_base(5);
t119 = qJD(5) * t141 + V_base(4);
t167 = (Icges(6,5) * t145 + Icges(6,6) * t146) * t143 + t118 * (-Icges(6,3) * t142 + t141 * t172) + t119 * (Icges(6,3) * t141 + t142 * t172);
t166 = V_base(5) * qJ(3) + t169;
t165 = qJD(4) * t141 + t166;
t164 = (Icges(5,5) * t157 + Icges(5,6) * t158) * t143 + (-Icges(5,3) * t142 + t141 * t173) * V_base(5) + (Icges(5,3) * t141 + t142 * t173) * V_base(4);
t106 = Icges(6,2) * t146 + t192;
t107 = Icges(6,1) * t145 + t191;
t76 = -Icges(6,6) * t142 + t141 * t174;
t77 = Icges(6,6) * t141 + t142 * t174;
t78 = -Icges(6,5) * t142 + t141 * t176;
t79 = Icges(6,5) * t141 + t142 * t176;
t163 = (-t145 * t77 + t146 * t79) * t119 + (-t145 * t76 + t146 * t78) * t118 + (-t106 * t145 + t107 * t146) * t143;
t124 = Icges(5,2) * t158 + t194;
t125 = Icges(5,1) * t157 + t193;
t84 = -Icges(5,6) * t142 + t141 * t175;
t85 = Icges(5,6) * t141 + t142 * t175;
t86 = -Icges(5,5) * t142 + t141 * t177;
t87 = Icges(5,5) * t141 + t142 * t177;
t162 = (-t157 * t85 + t158 * t87) * V_base(4) + (-t157 * t84 + t158 * t86) * V_base(5) + (-t124 * t157 + t125 * t158) * t143;
t152 = Icges(2,4) * t161;
t140 = Icges(3,4) * t150;
t138 = Icges(4,4) * t142;
t134 = rSges(2,1) * t161 - rSges(2,2) * t160;
t133 = rSges(2,1) * t160 + rSges(2,2) * t161;
t132 = Icges(2,1) * t161 - t197;
t131 = Icges(2,1) * t160 + t152;
t130 = -Icges(2,2) * t160 + t152;
t129 = Icges(2,2) * t161 + t197;
t126 = rSges(5,1) * t157 + rSges(5,2) * t158;
t122 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t121 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t120 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t116 = rSges(3,1) * t150 - rSges(3,2) * t149;
t115 = rSges(3,1) * t149 + rSges(3,2) * t150;
t114 = Icges(3,1) * t150 - t196;
t113 = Icges(3,1) * t149 + t140;
t112 = -Icges(3,2) * t149 + t140;
t111 = Icges(3,2) * t150 + t196;
t108 = rSges(6,1) * t145 + rSges(6,2) * t146;
t104 = rSges(4,1) * t142 - rSges(4,2) * t141;
t102 = rSges(4,1) * t141 + rSges(4,2) * t142;
t100 = Icges(4,1) * t142 - t195;
t99 = Icges(4,1) * t141 + t138;
t98 = -Icges(4,2) * t141 + t138;
t97 = Icges(4,2) * t142 + t195;
t92 = V_base(5) * rSges(2,3) - t133 * t148 + t186;
t91 = t134 * t148 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t133 * V_base(4) - t134 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t141 + t142 * t179;
t88 = -rSges(5,3) * t142 + t141 * t179;
t81 = rSges(6,3) * t141 + t142 * t178;
t80 = -rSges(6,3) * t142 + t141 * t178;
t73 = pkin(7) * t141 + t142 * t198;
t72 = -pkin(7) * t142 + t141 * t198;
t71 = V_base(5) * rSges(3,3) - t115 * t143 + t169;
t70 = t116 * t143 + (-rSges(3,3) + t205) * V_base(4) + t188;
t69 = V_base(4) * t115 + (-t116 - t202) * V_base(5) + t187;
t68 = V_base(5) * rSges(4,3) + (-t102 - t201) * t143 + t166;
t67 = t104 * t143 + (-rSges(4,3) + t189) * V_base(4) + t182;
t66 = V_base(4) * t102 + (-t104 + t181) * V_base(5) + t180;
t65 = t126 * V_base(5) + (t183 - t88) * t143 + t165;
t64 = t143 * t89 + (-t126 + t189) * V_base(4) + t168;
t63 = V_base(4) * t88 + (t170 - t89) * V_base(5) + t171;
t62 = V_base(5) * t199 + t108 * t118 + (t183 - t72 - t80) * t143 + t165;
t61 = -t108 * t119 + (t73 + t81) * t143 + (t189 - t199) * V_base(4) + t168;
t60 = -t118 * t81 + t119 * t80 + V_base(4) * t72 + (t170 - t73) * V_base(5) + t171;
t1 = m(1) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t119 * (t167 * t141 + t163 * t142) / 0.2e1 + t118 * (t163 * t141 - t167 * t142) / 0.2e1 + ((t145 * t79 + t146 * t77) * t119 + (t145 * t78 + t146 * t76) * t118 + (t157 * t86 + t158 * t84 + t207) * V_base(5) + (t157 * t87 + t158 * t85 + t206) * V_base(4) + (t106 * t146 + t107 * t145 + t124 * t158 + t125 * t157 + Icges(3,3) + Icges(4,3)) * t143) * t143 / 0.2e1 + (t164 * t141 + t162 * t142 + t206 * t143 + (-t111 * t149 + t113 * t150 - t129 * t160 + t131 * t161 - t141 * t97 + t142 * t99 + Icges(1,4)) * V_base(5) + (t100 * t142 - t149 * t112 + t114 * t150 - t160 * t130 + t132 * t161 - t141 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t162 * t141 - t164 * t142 + t207 * t143 + (t111 * t150 + t113 * t149 + t129 * t161 + t160 * t131 + t141 * t99 + t142 * t97 + Icges(1,2)) * V_base(5) + (t100 * t141 + t112 * t150 + t114 * t149 + t130 * t161 + t132 * t160 + t142 * t98 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t160 + Icges(2,6) * t161) * V_base(5) + (Icges(2,5) * t161 - Icges(2,6) * t160) * V_base(4) + Icges(2,3) * t148 / 0.2e1) * t148;
T = t1;
