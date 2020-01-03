% Calculate kinetic energy for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:26
% EndTime: 2019-12-31 19:26:27
% DurationCPUTime: 1.34s
% Computational Cost: add. (871->203), mult. (602->249), div. (0->0), fcn. (380->8), ass. (0->103)
t196 = Icges(4,4) + Icges(5,6);
t195 = Icges(4,1) + Icges(5,2);
t194 = -Icges(5,4) + Icges(4,5);
t193 = Icges(5,5) - Icges(4,6);
t192 = Icges(4,2) + Icges(5,3);
t140 = qJ(1) + qJ(2);
t132 = pkin(8) + t140;
t130 = cos(t132);
t191 = t196 * t130;
t129 = sin(t132);
t190 = t196 * t129;
t189 = -t130 * t192 - t190;
t188 = t129 * t192 - t191;
t187 = t129 * t195 + t191;
t186 = t130 * t195 - t190;
t106 = qJD(5) * t129 + V_base(5);
t107 = qJD(5) * t130 + V_base(4);
t141 = sin(qJ(5));
t143 = cos(qJ(5));
t171 = Icges(6,4) * t143;
t114 = -Icges(6,2) * t141 + t171;
t172 = Icges(6,4) * t141;
t117 = Icges(6,1) * t143 - t172;
t133 = V_base(6) + qJD(1);
t131 = qJD(2) + t133;
t153 = Icges(6,2) * t143 + t172;
t68 = Icges(6,6) * t130 + t129 * t153;
t69 = Icges(6,6) * t129 - t130 * t153;
t154 = Icges(6,1) * t141 + t171;
t70 = Icges(6,5) * t130 + t129 * t154;
t71 = Icges(6,5) * t129 - t130 * t154;
t183 = (t141 * t70 + t143 * t68) * t107 + (t141 * t71 + t143 * t69) * t106 + (t114 * t143 + t117 * t141) * t131;
t182 = -pkin(5) - pkin(6);
t142 = sin(qJ(1));
t180 = pkin(1) * t142;
t144 = cos(qJ(1));
t179 = pkin(1) * t144;
t134 = sin(t140);
t178 = pkin(2) * t134;
t135 = cos(t140);
t177 = pkin(2) * t135;
t176 = pkin(7) * t129;
t175 = Icges(2,4) * t142;
t174 = Icges(3,4) * t134;
t168 = -qJ(3) + t182;
t167 = t133 * t179 + V_base(2);
t166 = V_base(4) * t180 + V_base(3);
t165 = V_base(5) * pkin(5) + V_base(1);
t91 = pkin(3) * t129 - qJ(4) * t130;
t162 = -t91 - t178;
t161 = t131 * t177 + t167;
t94 = pkin(3) * t130 + qJ(4) * t129;
t160 = t131 * t94 + t161;
t159 = -t177 - t179;
t158 = V_base(4) * t178 + qJD(3) + t166;
t157 = rSges(6,1) * t141 + rSges(6,2) * t143;
t152 = Icges(6,5) * t141 + Icges(6,6) * t143;
t150 = t159 - t94;
t149 = V_base(4) * t91 + t158;
t148 = V_base(5) * pkin(6) - t133 * t180 + t165;
t147 = (Icges(6,3) * t129 - t130 * t152) * t106 + (Icges(6,3) * t130 + t129 * t152) * t107 + (Icges(6,5) * t143 - Icges(6,6) * t141) * t131;
t146 = V_base(5) * qJ(3) + t148;
t145 = qJD(4) * t129 + t146;
t137 = Icges(2,4) * t144;
t128 = Icges(3,4) * t135;
t122 = rSges(2,1) * t144 - t142 * rSges(2,2);
t121 = rSges(6,1) * t143 - rSges(6,2) * t141;
t120 = t142 * rSges(2,1) + rSges(2,2) * t144;
t119 = Icges(2,1) * t144 - t175;
t118 = Icges(2,1) * t142 + t137;
t116 = -Icges(2,2) * t142 + t137;
t115 = Icges(2,2) * t144 + t175;
t110 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t109 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t108 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t104 = rSges(3,1) * t135 - rSges(3,2) * t134;
t103 = rSges(3,1) * t134 + rSges(3,2) * t135;
t102 = Icges(3,1) * t135 - t174;
t101 = Icges(3,1) * t134 + t128;
t100 = -Icges(3,2) * t134 + t128;
t99 = Icges(3,2) * t135 + t174;
t96 = rSges(4,1) * t130 - rSges(4,2) * t129;
t95 = -rSges(5,2) * t130 + rSges(5,3) * t129;
t93 = rSges(4,1) * t129 + rSges(4,2) * t130;
t92 = -rSges(5,2) * t129 - rSges(5,3) * t130;
t76 = V_base(5) * rSges(2,3) - t120 * t133 + t165;
t75 = t122 * t133 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t74 = t120 * V_base(4) - t122 * V_base(5) + V_base(3);
t73 = rSges(6,3) * t129 - t130 * t157;
t72 = rSges(6,3) * t130 + t129 * t157;
t65 = V_base(5) * rSges(3,3) - t103 * t131 + t148;
t64 = t104 * t131 + (-rSges(3,3) + t182) * V_base(4) + t167;
t63 = V_base(4) * t103 + (-t104 - t179) * V_base(5) + t166;
t62 = V_base(5) * rSges(4,3) + (-t93 - t178) * t131 + t146;
t61 = t131 * t96 + (-rSges(4,3) + t168) * V_base(4) + t161;
t60 = V_base(4) * t93 + (t159 - t96) * V_base(5) + t158;
t59 = V_base(5) * rSges(5,1) + (t162 - t92) * t131 + t145;
t58 = -qJD(4) * t130 + t131 * t95 + (-rSges(5,1) + t168) * V_base(4) + t160;
t57 = V_base(4) * t92 + (t150 - t95) * V_base(5) + t149;
t56 = V_base(5) * pkin(4) + t106 * t121 + (t162 - t73 - t176) * t131 + t145;
t55 = -t107 * t121 + t131 * t72 + (pkin(7) * t131 - qJD(4)) * t130 + (-pkin(4) + t168) * V_base(4) + t160;
t54 = V_base(4) * t176 - t106 * t72 + t107 * t73 + (-pkin(7) * t130 + t150) * V_base(5) + t149;
t1 = m(1) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + t107 * (t183 * t129 + t147 * t130) / 0.2e1 + t106 * (t147 * t129 - t183 * t130) / 0.2e1 + ((-t141 * t68 + t143 * t70) * t107 + (-t141 * t69 + t143 * t71) * t106 + (-t141 * t114 + t143 * t117 + Icges(5,1) + Icges(3,3) + Icges(4,3)) * t131) * t131 / 0.2e1 + ((t101 * t135 - t142 * t115 + t118 * t144 + t129 * t189 + t187 * t130 - t134 * t99 + Icges(1,4)) * V_base(5) + (-t134 * t100 + t135 * t102 - t142 * t116 + t144 * t119 + t188 * t129 + t186 * t130 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t134 * t101 + t144 * t115 + t142 * t118 + t187 * t129 - t189 * t130 + t135 * t99 + Icges(1,2)) * V_base(5) + (t100 * t135 + t102 * t134 + t116 * t144 + t142 * t119 + t129 * t186 - t130 * t188 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t131 * (Icges(3,5) * t134 + Icges(3,6) * t135 + t129 * t194 - t130 * t193) + V_base(4) * t131 * (Icges(3,5) * t135 - Icges(3,6) * t134 + t129 * t193 + t130 * t194) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t142 + Icges(2,6) * t144) * V_base(5) + (Icges(2,5) * t144 - Icges(2,6) * t142) * V_base(4) + Icges(2,3) * t133 / 0.2e1) * t133;
T = t1;
