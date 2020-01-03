% Calculate kinetic energy for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:38
% EndTime: 2019-12-31 19:27:39
% DurationCPUTime: 1.45s
% Computational Cost: add. (968->204), mult. (884->253), div. (0->0), fcn. (804->8), ass. (0->103)
t204 = Icges(3,4) - Icges(4,5);
t203 = Icges(3,1) + Icges(4,1);
t202 = Icges(4,4) + Icges(3,5);
t201 = Icges(3,2) + Icges(4,3);
t200 = Icges(3,6) - Icges(4,6);
t179 = qJ(1) + qJ(2);
t173 = sin(t179);
t199 = t204 * t173;
t147 = cos(t179);
t198 = t204 * t147;
t197 = -t201 * t147 - t199;
t196 = t201 * t173 - t198;
t195 = t203 * t173 + t198;
t194 = t203 * t147 - t199;
t191 = -pkin(5) - pkin(6);
t153 = sin(qJ(1));
t189 = pkin(1) * t153;
t155 = cos(qJ(1));
t188 = pkin(1) * t155;
t187 = pkin(3) * t147;
t184 = sin(pkin(8));
t185 = cos(pkin(8));
t99 = -t147 * t185 - t173 * t184;
t186 = Icges(5,4) * t99;
t183 = Icges(2,4) * t153;
t152 = sin(qJ(5));
t182 = Icges(6,4) * t152;
t154 = cos(qJ(5));
t181 = Icges(6,4) * t154;
t146 = V_base(6) + qJD(1);
t178 = t146 * t188 + V_base(2);
t177 = V_base(4) * t189 + V_base(3);
t176 = V_base(5) * pkin(5) + V_base(1);
t117 = t147 * pkin(2) + qJ(3) * t173;
t172 = -t117 - t188;
t114 = pkin(2) * t173 - t147 * qJ(3);
t171 = V_base(4) * t114 + t177;
t170 = t173 * pkin(3);
t168 = -rSges(6,1) * t154 + rSges(6,2) * t152;
t167 = -Icges(6,1) * t154 + t182;
t166 = Icges(6,2) * t152 - t181;
t165 = -Icges(6,5) * t154 + Icges(6,6) * t152;
t164 = t172 - t187;
t163 = V_base(4) * t170 - qJD(4) + t171;
t162 = -t114 - t170;
t145 = qJD(2) + t146;
t161 = -qJD(3) * t147 + t145 * t117 + t178;
t160 = V_base(5) * pkin(6) - t146 * t189 + t176;
t100 = t147 * t184 - t173 * t185;
t95 = -qJD(5) * t99 + V_base(5);
t96 = qJD(5) * t100 + V_base(4);
t159 = (-Icges(6,5) * t152 - Icges(6,6) * t154) * t145 + (-Icges(6,3) * t99 + t100 * t165) * t95 + (Icges(6,3) * t100 + t165 * t99) * t96;
t158 = qJD(3) * t173 + t160;
t157 = V_base(4) * qJ(4) + t145 * t187 + t161;
t127 = -Icges(6,2) * t154 - t182;
t130 = -Icges(6,1) * t152 - t181;
t73 = -Icges(6,6) * t99 + t100 * t166;
t74 = Icges(6,6) * t100 + t166 * t99;
t75 = -Icges(6,5) * t99 + t100 * t167;
t76 = Icges(6,5) * t100 + t167 * t99;
t156 = (t152 * t74 - t154 * t76) * t96 + (t152 * t73 - t154 * t75) * t95 + (t127 * t152 - t130 * t154) * t145;
t149 = Icges(2,4) * t155;
t135 = rSges(2,1) * t155 - t153 * rSges(2,2);
t134 = t153 * rSges(2,1) + rSges(2,2) * t155;
t133 = -rSges(6,1) * t152 - rSges(6,2) * t154;
t132 = Icges(2,1) * t155 - t183;
t131 = Icges(2,1) * t153 + t149;
t129 = -Icges(2,2) * t153 + t149;
t128 = Icges(2,2) * t155 + t183;
t123 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t122 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t121 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t119 = t147 * rSges(3,1) - rSges(3,2) * t173;
t118 = t147 * rSges(4,1) + rSges(4,3) * t173;
t116 = rSges(3,1) * t173 + t147 * rSges(3,2);
t115 = rSges(4,1) * t173 - t147 * rSges(4,3);
t98 = Icges(5,4) * t100;
t94 = V_base(5) * rSges(2,3) - t134 * t146 + t176;
t93 = t135 * t146 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t92 = t134 * V_base(4) - t135 * V_base(5) + V_base(3);
t91 = -pkin(4) * t99 + pkin(7) * t100;
t90 = -pkin(4) * t100 - pkin(7) * t99;
t89 = -rSges(5,1) * t99 - rSges(5,2) * t100;
t88 = -rSges(5,1) * t100 + rSges(5,2) * t99;
t87 = -Icges(5,1) * t99 - t98;
t86 = -Icges(5,1) * t100 + t186;
t85 = -Icges(5,2) * t100 - t186;
t84 = Icges(5,2) * t99 - t98;
t81 = V_base(5) * rSges(3,3) - t116 * t145 + t160;
t80 = t119 * t145 + (-rSges(3,3) + t191) * V_base(4) + t178;
t79 = V_base(4) * t116 + (-t119 - t188) * V_base(5) + t177;
t78 = rSges(6,3) * t100 + t168 * t99;
t77 = -rSges(6,3) * t99 + t100 * t168;
t70 = V_base(5) * rSges(4,2) + (-t114 - t115) * t145 + t158;
t69 = t118 * t145 + (-rSges(4,2) + t191) * V_base(4) + t161;
t68 = V_base(4) * t115 + (-t118 + t172) * V_base(5) + t171;
t67 = (-rSges(5,3) - qJ(4)) * V_base(5) + (-t88 + t162) * t145 + t158;
t66 = t145 * t89 + (rSges(5,3) + t191) * V_base(4) + t157;
t65 = V_base(4) * t88 + (t164 - t89) * V_base(5) + t163;
t64 = -V_base(5) * qJ(4) + t95 * t133 + (-t77 - t90 + t162) * t145 + t158;
t63 = -t133 * t96 + t191 * V_base(4) + (t78 + t91) * t145 + t157;
t62 = t96 * t77 - t95 * t78 + V_base(4) * t90 + (t164 - t91) * V_base(5) + t163;
t1 = m(1) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t96 * (t100 * t159 + t156 * t99) / 0.2e1 + t95 * (t100 * t156 - t159 * t99) / 0.2e1 + ((-t152 * t76 - t154 * t74) * t96 + (-t152 * t75 - t154 * t73) * t95 + (-t154 * t127 - t152 * t130 + Icges(4,2) + Icges(3,3) + Icges(5,3)) * t145) * t145 / 0.2e1 + ((-t100 * t84 - t153 * t128 + t131 * t155 + t195 * t147 + t173 * t197 - t86 * t99 + Icges(1,4)) * V_base(5) + (-t100 * t85 - t153 * t129 + t155 * t132 + t194 * t147 + t196 * t173 - t99 * t87 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t100 * t86 + t155 * t128 + t153 * t131 - t197 * t147 + t195 * t173 + t99 * t84 + Icges(1,2)) * V_base(5) + (-t100 * t87 + t129 * t155 + t153 * t132 - t147 * t196 + t173 * t194 + t85 * t99 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t145 * (Icges(5,5) * t100 - Icges(5,6) * t99 + t200 * t147 + t202 * t173) + V_base(4) * t145 * (Icges(5,5) * t99 + Icges(5,6) * t100 + t202 * t147 - t200 * t173) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t153 + Icges(2,6) * t155) * V_base(5) + (Icges(2,5) * t155 - Icges(2,6) * t153) * V_base(4) + Icges(2,3) * t146 / 0.2e1) * t146;
T = t1;
