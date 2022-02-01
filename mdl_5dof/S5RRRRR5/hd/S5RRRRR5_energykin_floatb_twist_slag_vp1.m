% Calculate kinetic energy for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:37
% EndTime: 2022-01-20 12:01:38
% DurationCPUTime: 1.43s
% Computational Cost: add. (1207->228), mult. (760->313), div. (0->0), fcn. (538->10), ass. (0->119)
t205 = -pkin(5) - pkin(6);
t161 = sin(qJ(1));
t203 = pkin(1) * t161;
t163 = cos(qJ(1));
t202 = pkin(1) * t163;
t159 = qJ(1) + qJ(2);
t150 = sin(t159);
t201 = pkin(2) * t150;
t152 = cos(t159);
t200 = pkin(2) * t152;
t160 = sin(qJ(4));
t199 = pkin(4) * t160;
t162 = cos(qJ(4));
t198 = pkin(4) * t162;
t196 = Icges(2,4) * t161;
t195 = Icges(3,4) * t150;
t154 = qJ(3) + t159;
t145 = sin(t154);
t194 = Icges(4,4) * t145;
t193 = Icges(5,4) * t160;
t192 = Icges(5,4) * t162;
t158 = qJ(4) + qJ(5);
t149 = sin(t158);
t191 = Icges(6,4) * t149;
t151 = cos(t158);
t190 = Icges(6,4) * t151;
t189 = -pkin(7) + t205;
t148 = V_base(6) + qJD(1);
t188 = t148 * t202 + V_base(2);
t187 = V_base(4) * t203 + V_base(3);
t186 = V_base(5) * pkin(5) + V_base(1);
t124 = qJD(4) * t145 + V_base(4);
t144 = qJD(2) + t148;
t183 = t144 * t200 + t188;
t182 = V_base(4) * t201 + t187;
t181 = -t200 - t202;
t180 = rSges(5,1) * t162 - rSges(5,2) * t160;
t179 = rSges(6,1) * t151 - rSges(6,2) * t149;
t178 = Icges(5,1) * t162 - t193;
t177 = Icges(6,1) * t151 - t191;
t176 = -Icges(5,2) * t160 + t192;
t175 = -Icges(6,2) * t149 + t190;
t174 = Icges(5,5) * t162 - Icges(5,6) * t160;
t173 = Icges(6,5) * t151 - Icges(6,6) * t149;
t172 = V_base(5) * pkin(6) - t148 * t203 + t186;
t141 = qJD(3) + t144;
t146 = cos(t154);
t94 = V_base(5) + (-qJD(4) - qJD(5)) * t146;
t95 = qJD(5) * t145 + t124;
t171 = (Icges(6,5) * t149 + Icges(6,6) * t151) * t141 + (-Icges(6,3) * t146 + t145 * t173) * t94 + (Icges(6,3) * t145 + t146 * t173) * t95;
t123 = -qJD(4) * t146 + V_base(5);
t170 = t123 * (-Icges(5,3) * t146 + t145 * t174) + t124 * (Icges(5,3) * t145 + t146 * t174) + (Icges(5,5) * t160 + Icges(5,6) * t162) * t141;
t106 = t146 * pkin(3) + t145 * pkin(8);
t169 = t141 * t106 + t189 * V_base(4) + t183;
t168 = V_base(5) * pkin(7) - t144 * t201 + t172;
t105 = t145 * pkin(3) - t146 * pkin(8);
t167 = V_base(4) * t105 + (-t106 + t181) * V_base(5) + t182;
t110 = Icges(6,2) * t151 + t191;
t113 = Icges(6,1) * t149 + t190;
t76 = -Icges(6,6) * t146 + t145 * t175;
t77 = Icges(6,6) * t145 + t146 * t175;
t78 = -Icges(6,5) * t146 + t145 * t177;
t79 = Icges(6,5) * t145 + t146 * t177;
t166 = (-t149 * t77 + t151 * t79) * t95 + (-t149 * t76 + t151 * t78) * t94 + (-t110 * t149 + t113 * t151) * t141;
t128 = Icges(5,2) * t162 + t193;
t131 = Icges(5,1) * t160 + t192;
t84 = -Icges(5,6) * t146 + t145 * t176;
t85 = Icges(5,6) * t145 + t146 * t176;
t86 = -Icges(5,5) * t146 + t145 * t178;
t87 = Icges(5,5) * t145 + t146 * t178;
t165 = (-t160 * t85 + t162 * t87) * t124 + (-t160 * t84 + t162 * t86) * t123 + (-t128 * t160 + t131 * t162) * t141;
t153 = Icges(2,4) * t163;
t143 = Icges(3,4) * t152;
t140 = Icges(4,4) * t146;
t136 = rSges(2,1) * t163 - rSges(2,2) * t161;
t135 = rSges(2,1) * t161 + rSges(2,2) * t163;
t134 = rSges(5,1) * t160 + rSges(5,2) * t162;
t133 = Icges(2,1) * t163 - t196;
t132 = Icges(2,1) * t161 + t153;
t130 = -Icges(2,2) * t161 + t153;
t129 = Icges(2,2) * t163 + t196;
t122 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t121 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t120 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t118 = rSges(3,1) * t152 - rSges(3,2) * t150;
t117 = rSges(3,1) * t150 + rSges(3,2) * t152;
t116 = rSges(6,1) * t149 + rSges(6,2) * t151;
t115 = Icges(3,1) * t152 - t195;
t114 = Icges(3,1) * t150 + t143;
t112 = -Icges(3,2) * t150 + t143;
t111 = Icges(3,2) * t152 + t195;
t104 = rSges(4,1) * t146 - rSges(4,2) * t145;
t103 = rSges(4,1) * t145 + rSges(4,2) * t146;
t102 = Icges(4,1) * t146 - t194;
t101 = Icges(4,1) * t145 + t140;
t100 = -Icges(4,2) * t145 + t140;
t99 = Icges(4,2) * t146 + t194;
t92 = V_base(5) * rSges(2,3) - t135 * t148 + t186;
t91 = t136 * t148 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t135 * V_base(4) - t136 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t145 + t146 * t180;
t88 = -rSges(5,3) * t146 + t145 * t180;
t81 = rSges(6,3) * t145 + t146 * t179;
t80 = -rSges(6,3) * t146 + t145 * t179;
t73 = pkin(9) * t145 + t146 * t198;
t72 = -pkin(9) * t146 + t145 * t198;
t71 = V_base(5) * rSges(3,3) - t117 * t144 + t172;
t70 = t118 * t144 + (-rSges(3,3) + t205) * V_base(4) + t188;
t69 = t117 * V_base(4) + (-t118 - t202) * V_base(5) + t187;
t68 = V_base(5) * rSges(4,3) - t103 * t141 + t168;
t67 = t104 * t141 + (-rSges(4,3) + t189) * V_base(4) + t183;
t66 = t103 * V_base(4) + (-t104 + t181) * V_base(5) + t182;
t65 = t123 * t134 + (-t105 - t88) * t141 + t168;
t64 = -t124 * t134 + t141 * t89 + t169;
t63 = -t123 * t89 + t124 * t88 + t167;
t62 = t123 * t199 + t116 * t94 + (-t105 - t72 - t80) * t141 + t168;
t61 = -t124 * t199 - t116 * t95 + (t73 + t81) * t141 + t169;
t60 = -t123 * t73 + t124 * t72 + t80 * t95 - t81 * t94 + t167;
t1 = m(1) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t124 * (t170 * t145 + t165 * t146) / 0.2e1 + t123 * (t165 * t145 - t170 * t146) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t95 * (t171 * t145 + t166 * t146) / 0.2e1 + t94 * (t166 * t145 - t171 * t146) / 0.2e1 + ((t160 * t87 + t162 * t85) * t124 + (t160 * t86 + t162 * t84) * t123 + (t149 * t79 + t151 * t77) * t95 + (t149 * t78 + t151 * t76) * t94 + (t110 * t151 + t113 * t149 + t128 * t162 + t131 * t160 + Icges(4,3)) * t141) * t141 / 0.2e1 + ((t101 * t146 - t111 * t150 + t114 * t152 - t129 * t161 + t132 * t163 - t145 * t99 + Icges(1,4)) * V_base(5) + (-t100 * t145 + t102 * t146 - t150 * t112 + t115 * t152 - t130 * t161 + t133 * t163 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t145 + t111 * t152 + t114 * t150 + t129 * t163 + t132 * t161 + t146 * t99 + Icges(1,2)) * V_base(5) + (t100 * t146 + t102 * t145 + t112 * t152 + t115 * t150 + t130 * t163 + t133 * t161 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t141 * (Icges(4,5) * t146 - Icges(4,6) * t145) + V_base(5) * t141 * (Icges(4,5) * t145 + Icges(4,6) * t146) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t161 + Icges(2,6) * t163) * V_base(5) + (Icges(2,5) * t163 - Icges(2,6) * t161) * V_base(4) + Icges(2,3) * t148 / 0.2e1) * t148 + ((Icges(3,5) * t150 + Icges(3,6) * t152) * V_base(5) + (Icges(3,5) * t152 - Icges(3,6) * t150) * V_base(4) + Icges(3,3) * t144 / 0.2e1) * t144;
T = t1;
