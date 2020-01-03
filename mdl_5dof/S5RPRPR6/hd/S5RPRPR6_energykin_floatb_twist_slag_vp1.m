% Calculate kinetic energy for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:40
% EndTime: 2019-12-31 18:17:41
% DurationCPUTime: 1.37s
% Computational Cost: add. (860->203), mult. (602->250), div. (0->0), fcn. (380->8), ass. (0->102)
t198 = Icges(4,4) + Icges(5,6);
t197 = Icges(4,1) + Icges(5,2);
t196 = -Icges(5,4) + Icges(4,5);
t195 = Icges(5,5) - Icges(4,6);
t194 = Icges(4,2) + Icges(5,3);
t140 = qJ(1) + pkin(8);
t134 = qJ(3) + t140;
t130 = cos(t134);
t193 = t198 * t130;
t129 = sin(t134);
t192 = t198 * t129;
t191 = -t130 * t194 - t192;
t190 = t129 * t194 - t193;
t189 = t129 * t197 + t193;
t188 = t130 * t197 - t192;
t105 = qJD(5) * t129 + V_base(5);
t106 = qJD(5) * t130 + V_base(4);
t141 = sin(qJ(5));
t143 = cos(qJ(5));
t170 = Icges(6,4) * t143;
t114 = -Icges(6,2) * t141 + t170;
t171 = Icges(6,4) * t141;
t117 = Icges(6,1) * t143 - t171;
t135 = V_base(6) + qJD(1);
t131 = qJD(3) + t135;
t152 = Icges(6,2) * t143 + t171;
t68 = Icges(6,6) * t130 + t129 * t152;
t69 = Icges(6,6) * t129 - t130 * t152;
t153 = Icges(6,1) * t141 + t170;
t70 = Icges(6,5) * t130 + t129 * t153;
t71 = Icges(6,5) * t129 - t130 * t153;
t183 = (t141 * t70 + t143 * t68) * t106 + (t141 * t71 + t143 * t69) * t105 + (t114 * t143 + t117 * t141) * t131;
t142 = sin(qJ(1));
t180 = pkin(1) * t142;
t144 = cos(qJ(1));
t179 = pkin(1) * t144;
t132 = sin(t140);
t178 = pkin(2) * t132;
t133 = cos(t140);
t177 = pkin(2) * t133;
t176 = pkin(7) * t129;
t175 = -pkin(5) - qJ(2);
t174 = Icges(2,4) * t142;
t173 = Icges(3,4) * t132;
t167 = -pkin(6) + t175;
t166 = t135 * t179 + V_base(2);
t165 = V_base(5) * pkin(5) + V_base(1);
t162 = t135 * t177 + t166;
t161 = V_base(5) * qJ(2) + t165;
t160 = t180 * V_base(4) + qJD(2) + V_base(3);
t94 = pkin(3) * t130 + qJ(4) * t129;
t159 = t131 * t94 + t162;
t158 = -t177 - t179;
t157 = t178 * V_base(4) + t160;
t156 = rSges(6,1) * t141 + rSges(6,2) * t143;
t151 = Icges(6,5) * t141 + Icges(6,6) * t143;
t149 = t158 - t94;
t91 = pkin(3) * t129 - qJ(4) * t130;
t148 = t91 * V_base(4) + t157;
t147 = t105 * (Icges(6,3) * t129 - t130 * t151) + t106 * (Icges(6,3) * t130 + t129 * t151) + (Icges(6,5) * t143 - Icges(6,6) * t141) * t131;
t146 = V_base(5) * pkin(6) + (-t178 - t180) * t135 + t161;
t145 = qJD(4) * t129 + t146;
t137 = Icges(2,4) * t144;
t128 = Icges(3,4) * t133;
t122 = rSges(2,1) * t144 - rSges(2,2) * t142;
t121 = rSges(6,1) * t143 - rSges(6,2) * t141;
t120 = rSges(2,1) * t142 + rSges(2,2) * t144;
t119 = Icges(2,1) * t144 - t174;
t118 = Icges(2,1) * t142 + t137;
t116 = -Icges(2,2) * t142 + t137;
t115 = Icges(2,2) * t144 + t174;
t109 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t108 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t107 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t104 = rSges(3,1) * t133 - rSges(3,2) * t132;
t103 = rSges(3,1) * t132 + rSges(3,2) * t133;
t102 = Icges(3,1) * t133 - t173;
t101 = Icges(3,1) * t132 + t128;
t100 = -Icges(3,2) * t132 + t128;
t99 = Icges(3,2) * t133 + t173;
t96 = rSges(4,1) * t130 - rSges(4,2) * t129;
t95 = -rSges(5,2) * t130 + rSges(5,3) * t129;
t93 = rSges(4,1) * t129 + rSges(4,2) * t130;
t92 = -rSges(5,2) * t129 - rSges(5,3) * t130;
t76 = V_base(5) * rSges(2,3) - t120 * t135 + t165;
t75 = t122 * t135 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t74 = t120 * V_base(4) - t122 * V_base(5) + V_base(3);
t73 = rSges(6,3) * t129 - t130 * t156;
t72 = rSges(6,3) * t130 + t129 * t156;
t65 = V_base(5) * rSges(3,3) + (-t103 - t180) * t135 + t161;
t64 = t104 * t135 + (-rSges(3,3) + t175) * V_base(4) + t166;
t63 = V_base(4) * t103 + (-t104 - t179) * V_base(5) + t160;
t62 = V_base(5) * rSges(4,3) - t131 * t93 + t146;
t61 = t131 * t96 + (-rSges(4,3) + t167) * V_base(4) + t162;
t60 = V_base(4) * t93 + (t158 - t96) * V_base(5) + t157;
t59 = V_base(5) * rSges(5,1) + (-t91 - t92) * t131 + t145;
t58 = -qJD(4) * t130 + t131 * t95 + (-rSges(5,1) + t167) * V_base(4) + t159;
t57 = V_base(4) * t92 + (t149 - t95) * V_base(5) + t148;
t56 = V_base(5) * pkin(4) + t105 * t121 + (-t73 - t91 - t176) * t131 + t145;
t55 = -t106 * t121 + t131 * t72 + (pkin(7) * t131 - qJD(4)) * t130 + (-pkin(4) + t167) * V_base(4) + t159;
t54 = V_base(4) * t176 - t105 * t72 + t106 * t73 + (-pkin(7) * t130 + t149) * V_base(5) + t148;
t1 = m(1) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + t106 * (t183 * t129 + t147 * t130) / 0.2e1 + t105 * (t147 * t129 - t183 * t130) / 0.2e1 + ((-t141 * t68 + t143 * t70) * t106 + (-t141 * t69 + t143 * t71) * t105 + (-t141 * t114 + t143 * t117 + Icges(5,1) + Icges(4,3)) * t131) * t131 / 0.2e1 + ((t101 * t133 - t142 * t115 + t118 * t144 + t129 * t191 + t130 * t189 - t132 * t99 + Icges(1,4)) * V_base(5) + (-t132 * t100 + t133 * t102 - t142 * t116 + t144 * t119 + t190 * t129 + t188 * t130 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t132 * t101 + t144 * t115 + t142 * t118 + t189 * t129 - t191 * t130 + t133 * t99 + Icges(1,2)) * V_base(5) + (t100 * t133 + t102 * t132 + t116 * t144 + t142 * t119 + t129 * t188 - t130 * t190 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t131 * (t129 * t196 - t130 * t195) + V_base(4) * t131 * (t129 * t195 + t130 * t196) + ((Icges(2,5) * t144 + Icges(3,5) * t133 - Icges(2,6) * t142 - Icges(3,6) * t132) * V_base(4) + (Icges(2,5) * t142 + Icges(3,5) * t132 + Icges(2,6) * t144 + Icges(3,6) * t133) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t135) * t135 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
