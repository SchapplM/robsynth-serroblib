% Calculate kinetic energy for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:34
% EndTime: 2019-12-31 17:39:36
% DurationCPUTime: 1.61s
% Computational Cost: add. (935->203), mult. (884->257), div. (0->0), fcn. (804->8), ass. (0->104)
t209 = Icges(3,4) - Icges(4,5);
t208 = Icges(3,1) + Icges(4,1);
t207 = Icges(4,4) + Icges(3,5);
t206 = Icges(3,2) + Icges(4,3);
t205 = Icges(3,6) - Icges(4,6);
t182 = pkin(8) + qJ(2);
t173 = sin(t182);
t204 = t209 * t173;
t146 = cos(t182);
t203 = t209 * t146;
t200 = -t206 * t146 - t204;
t199 = t206 * t173 - t203;
t196 = t208 * t173 + t203;
t195 = t208 * t146 - t204;
t192 = cos(qJ(4));
t191 = sin(qJ(4));
t153 = cos(pkin(8));
t190 = pkin(1) * t153;
t189 = pkin(3) * t146;
t188 = -pkin(5) - qJ(1);
t152 = sin(pkin(8));
t187 = Icges(2,4) * t152;
t100 = -t146 * t192 - t173 * t191;
t186 = Icges(5,4) * t100;
t154 = sin(qJ(5));
t185 = Icges(6,4) * t154;
t155 = cos(qJ(5));
t184 = Icges(6,4) * t155;
t175 = pkin(1) * V_base(6);
t181 = t153 * t175 + V_base(2);
t180 = V_base(5) * qJ(1) + V_base(1);
t176 = qJD(1) + V_base(3);
t148 = V_base(6) + qJD(2);
t117 = t146 * pkin(2) + qJ(3) * t173;
t174 = -t117 - t190;
t172 = V_base(4) * t152 * pkin(1) + t176;
t114 = pkin(2) * t173 - t146 * qJ(3);
t171 = V_base(4) * t114 + t172;
t170 = t173 * pkin(3);
t169 = -rSges(6,1) * t155 + rSges(6,2) * t154;
t167 = -Icges(6,1) * t155 + t185;
t166 = Icges(6,2) * t154 - t184;
t165 = -Icges(6,5) * t155 + Icges(6,6) * t154;
t164 = V_base(4) * t170 + t171;
t163 = t174 - t189;
t162 = -qJD(3) * t146 + t148 * t117 + t181;
t101 = t146 * t191 - t173 * t192;
t145 = -qJD(4) + t148;
t95 = -qJD(5) * t100 + V_base(5);
t96 = qJD(5) * t101 + V_base(4);
t161 = (-Icges(6,5) * t154 - Icges(6,6) * t155) * t145 + (-Icges(6,3) * t100 + t101 * t165) * t95 + (Icges(6,3) * t101 + t100 * t165) * t96;
t160 = V_base(5) * pkin(5) - t152 * t175 + t180;
t159 = V_base(4) * pkin(6) + t148 * t189 + t162;
t158 = qJD(3) * t173 + t160;
t157 = (-t170 - t114) * t148 + t158;
t133 = -Icges(6,2) * t155 - t185;
t134 = -Icges(6,1) * t154 - t184;
t73 = -Icges(6,6) * t100 + t101 * t166;
t74 = Icges(6,6) * t101 + t100 * t166;
t75 = -Icges(6,5) * t100 + t101 * t167;
t76 = Icges(6,5) * t101 + t100 * t167;
t156 = (t154 * t74 - t155 * t76) * t96 + (t154 * t73 - t155 * t75) * t95 + (t133 * t154 - t134 * t155) * t145;
t147 = Icges(2,4) * t153;
t135 = -t154 * rSges(6,1) - rSges(6,2) * t155;
t131 = rSges(2,1) * t153 - rSges(2,2) * t152;
t130 = rSges(2,1) * t152 + rSges(2,2) * t153;
t129 = Icges(2,1) * t153 - t187;
t128 = Icges(2,1) * t152 + t147;
t127 = -Icges(2,2) * t152 + t147;
t126 = Icges(2,2) * t153 + t187;
t122 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t121 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t120 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t119 = t146 * rSges(3,1) - rSges(3,2) * t173;
t118 = t146 * rSges(4,1) + rSges(4,3) * t173;
t116 = rSges(3,1) * t173 + t146 * rSges(3,2);
t115 = rSges(4,1) * t173 - t146 * rSges(4,3);
t98 = Icges(5,4) * t101;
t94 = V_base(5) * rSges(2,3) - t130 * V_base(6) + t180;
t93 = t131 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t92 = t130 * V_base(4) - t131 * V_base(5) + t176;
t91 = -pkin(4) * t100 + pkin(7) * t101;
t90 = -pkin(4) * t101 - pkin(7) * t100;
t89 = -rSges(5,1) * t100 - rSges(5,2) * t101;
t88 = -rSges(5,1) * t101 + rSges(5,2) * t100;
t87 = -Icges(5,1) * t100 - t98;
t86 = -Icges(5,1) * t101 + t186;
t85 = -Icges(5,2) * t101 - t186;
t84 = Icges(5,2) * t100 - t98;
t81 = V_base(5) * rSges(3,3) - t116 * t148 + t160;
t80 = t119 * t148 + (-rSges(3,3) + t188) * V_base(4) + t181;
t79 = t116 * V_base(4) + (-t119 - t190) * V_base(5) + t172;
t78 = t101 * rSges(6,3) + t100 * t169;
t77 = -t100 * rSges(6,3) + t101 * t169;
t70 = V_base(5) * rSges(4,2) + (-t114 - t115) * t148 + t158;
t69 = t118 * t148 + (-rSges(4,2) + t188) * V_base(4) + t162;
t68 = t115 * V_base(4) + (-t118 + t174) * V_base(5) + t171;
t67 = -t145 * t88 + (-rSges(5,3) - pkin(6)) * V_base(5) + t157;
t66 = t145 * t89 + (rSges(5,3) + t188) * V_base(4) + t159;
t65 = t88 * V_base(4) + (t163 - t89) * V_base(5) + t164;
t64 = -V_base(5) * pkin(6) + t95 * t135 + (-t77 - t90) * t145 + t157;
t63 = -t135 * t96 + t188 * V_base(4) + (t78 + t91) * t145 + t159;
t62 = t77 * t96 - t78 * t95 + t90 * V_base(4) + (t163 - t91) * V_base(5) + t164;
t1 = m(1) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(2) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t96 * (t156 * t100 + t161 * t101) / 0.2e1 + t95 * (-t161 * t100 + t156 * t101) / 0.2e1 + ((-t154 * t76 - t155 * t74) * t96 + (-t154 * t75 - t155 * t73) * t95 + (-t155 * t133 - t154 * t134 + Icges(5,3)) * t145) * t145 / 0.2e1 + ((-t100 * t86 - t101 * t84 - t126 * t152 + t128 * t153 + t146 * t196 + t173 * t200 + Icges(1,4)) * V_base(5) + (-t100 * t87 - t101 * t85 - t152 * t127 + t153 * t129 + t195 * t146 + t199 * t173 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t100 * t84 - t101 * t86 + t153 * t126 + t152 * t128 - t200 * t146 + t196 * t173 + Icges(1,2)) * V_base(5) + (t100 * t85 - t101 * t87 + t127 * t153 + t129 * t152 - t146 * t199 + t173 * t195 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t145 * (Icges(5,5) * t100 + Icges(5,6) * t101) + V_base(5) * t145 * (Icges(5,5) * t101 - Icges(5,6) * t100) + ((Icges(2,5) * t152 + Icges(2,6) * t153 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t153 - Icges(2,6) * t152 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t148 + (-t205 * V_base(4) + t207 * V_base(5)) * t173 + (t205 * V_base(5) + t207 * V_base(4)) * t146) * t148;
T = t1;
