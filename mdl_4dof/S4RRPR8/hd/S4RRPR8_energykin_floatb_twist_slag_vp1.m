% Calculate kinetic energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:45
% EndTime: 2019-12-31 17:07:46
% DurationCPUTime: 1.81s
% Computational Cost: add. (550->197), mult. (1032->275), div. (0->0), fcn. (942->6), ass. (0->98)
t198 = Icges(3,4) - Icges(4,5);
t197 = Icges(3,1) + Icges(4,1);
t196 = Icges(3,2) + Icges(4,3);
t144 = cos(qJ(2));
t195 = t198 * t144;
t141 = sin(qJ(2));
t194 = t198 * t141;
t193 = Icges(4,4) + Icges(3,5);
t192 = Icges(3,6) - Icges(4,6);
t191 = t196 * t141 - t195;
t190 = t197 * t144 - t194;
t142 = sin(qJ(1));
t145 = cos(qJ(1));
t189 = t191 * t142 + t192 * t145;
t188 = -t192 * t142 + t191 * t145;
t187 = t190 * t142 - t193 * t145;
t186 = t193 * t142 + t190 * t145;
t185 = Icges(4,2) + Icges(3,3);
t184 = -t196 * t144 - t194;
t183 = t197 * t141 + t195;
t182 = -t192 * t141 + t193 * t144;
t132 = -qJD(2) * t145 + V_base(5);
t133 = qJD(2) * t142 + V_base(4);
t136 = V_base(6) + qJD(1);
t181 = (t184 * t141 + t183 * t144) * t136 + (t188 * t141 + t186 * t144) * t133 + (t189 * t141 + t187 * t144) * t132;
t180 = (t193 * t141 + t192 * t144) * t136 + (t185 * t142 + t182 * t145) * t133 + (t182 * t142 - t185 * t145) * t132;
t176 = pkin(3) * t141;
t175 = pkin(3) * t144;
t174 = Icges(2,4) * t142;
t161 = pkin(2) * t144 + qJ(3) * t141;
t100 = t161 * t142;
t130 = t142 * pkin(1) - pkin(5) * t145;
t169 = -t100 - t130;
t168 = qJD(3) * t141;
t167 = V_base(5) * pkin(4) + V_base(1);
t125 = pkin(2) * t141 - qJ(3) * t144;
t164 = t132 * t125 + t145 * t168 + t167;
t163 = rSges(3,1) * t144 - rSges(3,2) * t141;
t162 = rSges(4,1) * t144 + rSges(4,3) * t141;
t140 = sin(qJ(4));
t143 = cos(qJ(4));
t106 = -t140 * t144 + t141 * t143;
t154 = t140 * t141 + t143 * t144;
t131 = pkin(1) * t145 + t142 * pkin(5);
t153 = -V_base(4) * pkin(4) + t136 * t131 + V_base(2);
t152 = V_base(4) * t130 - t131 * V_base(5) + V_base(3);
t101 = t161 * t145;
t149 = t136 * t101 + t142 * t168 + t153;
t148 = -qJD(3) * t144 + t133 * t100 + t152;
t138 = Icges(2,4) * t145;
t129 = rSges(2,1) * t145 - t142 * rSges(2,2);
t128 = t142 * rSges(2,1) + rSges(2,2) * t145;
t127 = rSges(3,1) * t141 + rSges(3,2) * t144;
t126 = rSges(4,1) * t141 - rSges(4,3) * t144;
t124 = Icges(2,1) * t145 - t174;
t123 = Icges(2,1) * t142 + t138;
t120 = -Icges(2,2) * t142 + t138;
t119 = Icges(2,2) * t145 + t174;
t112 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t111 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t110 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t109 = -t142 * pkin(6) + t145 * t175;
t108 = pkin(6) * t145 + t142 * t175;
t104 = -qJD(4) * t142 + t133;
t103 = V_base(5) + (-qJD(2) + qJD(4)) * t145;
t99 = t154 * t145;
t98 = t106 * t145;
t97 = t154 * t142;
t96 = t106 * t142;
t94 = t142 * rSges(3,3) + t145 * t163;
t93 = t142 * rSges(4,2) + t145 * t162;
t92 = -rSges(3,3) * t145 + t142 * t163;
t91 = -rSges(4,2) * t145 + t142 * t162;
t76 = V_base(5) * rSges(2,3) - t128 * t136 + t167;
t75 = t129 * t136 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t74 = t128 * V_base(4) - t129 * V_base(5) + V_base(3);
t73 = rSges(5,1) * t106 - rSges(5,2) * t154;
t72 = Icges(5,1) * t106 - Icges(5,4) * t154;
t71 = Icges(5,4) * t106 - Icges(5,2) * t154;
t70 = Icges(5,5) * t106 - Icges(5,6) * t154;
t69 = rSges(5,1) * t99 + rSges(5,2) * t98 - rSges(5,3) * t142;
t68 = t97 * rSges(5,1) + t96 * rSges(5,2) + rSges(5,3) * t145;
t67 = Icges(5,1) * t99 + Icges(5,4) * t98 - Icges(5,5) * t142;
t66 = Icges(5,1) * t97 + Icges(5,4) * t96 + Icges(5,5) * t145;
t65 = Icges(5,4) * t99 + Icges(5,2) * t98 - Icges(5,6) * t142;
t64 = Icges(5,4) * t97 + Icges(5,2) * t96 + Icges(5,6) * t145;
t63 = Icges(5,5) * t99 + Icges(5,6) * t98 - Icges(5,3) * t142;
t62 = Icges(5,5) * t97 + Icges(5,6) * t96 + Icges(5,3) * t145;
t61 = t127 * t132 + (-t130 - t92) * t136 + t167;
t60 = -t127 * t133 + t136 * t94 + t153;
t59 = -t132 * t94 + t133 * t92 + t152;
t58 = t126 * t132 + (-t91 + t169) * t136 + t164;
t57 = t136 * t93 + (-t125 - t126) * t133 + t149;
t56 = t133 * t91 + (-t101 - t93) * t132 + t148;
t55 = t132 * t176 + t103 * t73 + (-t108 - t68 + t169) * t136 + t164;
t54 = -t104 * t73 + (t109 + t69) * t136 + (-t125 - t176) * t133 + t149;
t53 = -t103 * t69 + t104 * t68 + t108 * t133 + (-t101 - t109) * t132 + t148;
t1 = m(1) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + t104 * ((-t142 * t63 + t65 * t98 + t67 * t99) * t104 + (-t142 * t62 + t64 * t98 + t66 * t99) * t103 + (-t142 * t70 + t71 * t98 + t72 * t99) * t136) / 0.2e1 + t103 * ((t145 * t63 + t96 * t65 + t97 * t67) * t104 + (t145 * t62 + t96 * t64 + t97 * t66) * t103 + (t145 * t70 + t96 * t71 + t97 * t72) * t136) / 0.2e1 + (t181 * t142 - t180 * t145) * t132 / 0.2e1 + (t180 * t142 + t181 * t145) * t133 / 0.2e1 + ((-t142 * t119 + t123 * t145 + Icges(1,4)) * V_base(5) + (-t142 * t120 + t124 * t145 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t119 * t145 + t142 * t123 + Icges(1,2)) * V_base(5) + (t120 * t145 + t142 * t124 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t106 * t67 - t154 * t65) * t104 + (t106 * t66 - t154 * t64) * t103 + (t186 * t141 - t188 * t144) * t133 + (t187 * t141 - t189 * t144) * t132 + (t106 * t72 + t183 * t141 - t184 * t144 - t154 * t71 + Icges(2,3)) * t136) * t136 / 0.2e1 + V_base(4) * t136 * (Icges(2,5) * t145 - Icges(2,6) * t142) + V_base(5) * t136 * (Icges(2,5) * t142 + Icges(2,6) * t145) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
