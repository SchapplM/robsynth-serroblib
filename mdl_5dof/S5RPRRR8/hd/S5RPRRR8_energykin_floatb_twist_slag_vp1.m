% Calculate kinetic energy for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:41
% DurationCPUTime: 1.99s
% Computational Cost: add. (867->219), mult. (1312->298), div. (0->0), fcn. (1366->8), ass. (0->114)
t222 = Icges(2,4) - Icges(3,5);
t221 = Icges(2,1) + Icges(3,1);
t220 = Icges(3,4) + Icges(2,5);
t219 = Icges(2,2) + Icges(3,3);
t218 = Icges(2,6) - Icges(3,6);
t205 = sin(qJ(1));
t217 = t222 * t205;
t207 = cos(qJ(1));
t216 = t222 * t207;
t215 = -t207 * t219 - t217;
t214 = t205 * t219 - t216;
t211 = t205 * t221 + t216;
t210 = t207 * t221 - t217;
t206 = cos(qJ(3));
t204 = sin(qJ(3));
t164 = sin(qJ(4));
t203 = pkin(4) * t164;
t165 = cos(qJ(4));
t202 = pkin(4) * t165;
t121 = -t204 * t205 - t206 * t207;
t200 = Icges(4,4) * t121;
t199 = Icges(5,4) * t164;
t198 = Icges(5,4) * t165;
t163 = qJ(4) + qJ(5);
t157 = sin(t163);
t197 = Icges(6,4) * t157;
t158 = cos(t163);
t196 = Icges(6,4) * t158;
t142 = pkin(1) * t205 - qJ(2) * t207;
t195 = t142 * V_base(4) + V_base(3);
t194 = V_base(5) * pkin(5) + V_base(1);
t191 = t207 * pkin(2);
t190 = t205 * pkin(2);
t122 = t204 * t207 - t205 * t206;
t111 = qJD(4) * t122 + V_base(4);
t110 = -qJD(4) * t121 + V_base(5);
t155 = V_base(6) + qJD(1);
t187 = t190 * V_base(4) + t195;
t186 = qJD(2) * t205 + t194;
t145 = pkin(1) * t207 + qJ(2) * t205;
t185 = -t145 - t191;
t184 = -rSges(5,1) * t165 + rSges(5,2) * t164;
t183 = -rSges(6,1) * t158 + rSges(6,2) * t157;
t182 = -Icges(5,1) * t165 + t199;
t181 = -Icges(6,1) * t158 + t197;
t180 = Icges(5,2) * t164 - t198;
t179 = Icges(6,2) * t157 - t196;
t178 = -Icges(5,5) * t165 + Icges(5,6) * t164;
t177 = -Icges(6,5) * t158 + Icges(6,6) * t157;
t176 = -qJD(2) * t207 + t145 * t155 + V_base(2);
t150 = -qJD(3) + t155;
t94 = -qJD(5) * t121 + t110;
t95 = qJD(5) * t122 + t111;
t175 = (-Icges(6,5) * t157 - Icges(6,6) * t158) * t150 + (-Icges(6,3) * t121 + t122 * t177) * t94 + (Icges(6,3) * t122 + t121 * t177) * t95;
t174 = (-Icges(5,3) * t121 + t122 * t178) * t110 + (Icges(5,3) * t122 + t121 * t178) * t111 + (-Icges(5,5) * t164 - Icges(5,6) * t165) * t150;
t173 = V_base(4) * pkin(6) + t155 * t191 + t176;
t172 = (-t190 - t142) * t155 + t186;
t106 = -t121 * pkin(3) + t122 * pkin(7);
t171 = -V_base(4) * pkin(5) + t106 * t150 + t173;
t105 = -pkin(3) * t122 - t121 * pkin(7);
t170 = V_base(4) * t105 + (-t106 + t185) * V_base(5) + t187;
t169 = -V_base(5) * pkin(6) + t172;
t117 = -Icges(6,2) * t158 - t197;
t118 = -Icges(6,1) * t157 - t196;
t77 = -Icges(6,6) * t121 + t122 * t179;
t78 = Icges(6,6) * t122 + t121 * t179;
t79 = -Icges(6,5) * t121 + t122 * t181;
t80 = Icges(6,5) * t122 + t121 * t181;
t168 = (t157 * t78 - t158 * t80) * t95 + (t157 * t77 - t158 * t79) * t94 + (t117 * t157 - t118 * t158) * t150;
t131 = -Icges(5,2) * t165 - t199;
t136 = -Icges(5,1) * t164 - t198;
t85 = -Icges(5,6) * t121 + t122 * t180;
t86 = Icges(5,6) * t122 + t121 * t180;
t87 = -Icges(5,5) * t121 + t122 * t182;
t88 = Icges(5,5) * t122 + t121 * t182;
t167 = (t164 * t86 - t165 * t88) * t111 + (t164 * t85 - t165 * t87) * t110 + (t131 * t164 - t136 * t165) * t150;
t147 = rSges(2,1) * t207 - rSges(2,2) * t205;
t146 = rSges(3,1) * t207 + rSges(3,3) * t205;
t144 = rSges(2,1) * t205 + rSges(2,2) * t207;
t143 = rSges(3,1) * t205 - rSges(3,3) * t207;
t141 = -rSges(5,1) * t164 - rSges(5,2) * t165;
t125 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t124 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t123 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t119 = -rSges(6,1) * t157 - rSges(6,2) * t158;
t115 = Icges(4,4) * t122;
t109 = V_base(5) * rSges(2,3) - t144 * t155 + t194;
t108 = t147 * t155 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t107 = t144 * V_base(4) - t147 * V_base(5) + V_base(3);
t104 = -rSges(4,1) * t121 - rSges(4,2) * t122;
t103 = -rSges(4,1) * t122 + rSges(4,2) * t121;
t102 = -Icges(4,1) * t121 - t115;
t101 = -Icges(4,1) * t122 + t200;
t100 = -Icges(4,2) * t122 - t200;
t99 = Icges(4,2) * t121 - t115;
t92 = V_base(5) * rSges(3,2) + (-t142 - t143) * t155 + t186;
t91 = t155 * t146 + (-pkin(5) - rSges(3,2)) * V_base(4) + t176;
t90 = rSges(5,3) * t122 + t121 * t184;
t89 = -rSges(5,3) * t121 + t122 * t184;
t82 = rSges(6,3) * t122 + t121 * t183;
t81 = -rSges(6,3) * t121 + t122 * t183;
t74 = t143 * V_base(4) + (-t145 - t146) * V_base(5) + t195;
t73 = pkin(8) * t122 - t121 * t202;
t72 = -pkin(8) * t121 - t122 * t202;
t71 = -t150 * t103 + (-pkin(6) - rSges(4,3)) * V_base(5) + t172;
t70 = t150 * t104 + (rSges(4,3) - pkin(5)) * V_base(4) + t173;
t69 = V_base(4) * t103 + (-t104 + t185) * V_base(5) + t187;
t68 = t110 * t141 + (-t105 - t89) * t150 + t169;
t67 = -t111 * t141 + t150 * t90 + t171;
t66 = -t110 * t90 + t111 * t89 + t170;
t65 = -t110 * t203 + t94 * t119 + (-t105 - t72 - t81) * t150 + t169;
t64 = t111 * t203 - t95 * t119 + (t73 + t82) * t150 + t171;
t63 = -t110 * t73 + t111 * t72 + t81 * t95 - t82 * t94 + t170;
t1 = m(1) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t111 * (t167 * t121 + t174 * t122) / 0.2e1 + t110 * (-t174 * t121 + t167 * t122) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t95 * (t168 * t121 + t175 * t122) / 0.2e1 + t94 * (-t175 * t121 + t168 * t122) / 0.2e1 + ((-t164 * t88 - t165 * t86) * t111 + (-t164 * t87 - t165 * t85) * t110 + (-t157 * t80 - t158 * t78) * t95 + (-t157 * t79 - t158 * t77) * t94 + (-t158 * t117 - t157 * t118 - t165 * t131 - t164 * t136 + Icges(4,3)) * t150) * t150 / 0.2e1 + ((-t101 * t121 - t122 * t99 + t205 * t215 + t207 * t211 + Icges(1,4)) * V_base(5) + (-t122 * t100 - t121 * t102 + t214 * t205 + t210 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t122 * t101 + t121 * t99 + t211 * t205 - t215 * t207 + Icges(1,2)) * V_base(5) + (t100 * t121 - t102 * t122 + t205 * t210 - t207 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t150 * (Icges(4,5) * t121 + Icges(4,6) * t122) + V_base(5) * t150 * (Icges(4,5) * t122 - Icges(4,6) * t121) + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t155 + (t218 * V_base(5) + t220 * V_base(4)) * t207 + (-t218 * V_base(4) + t220 * V_base(5)) * t205) * t155 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
