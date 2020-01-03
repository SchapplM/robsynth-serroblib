% Calculate kinetic energy for
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:20
% EndTime: 2019-12-31 16:57:21
% DurationCPUTime: 1.38s
% Computational Cost: add. (662->173), mult. (821->227), div. (0->0), fcn. (653->6), ass. (0->91)
t201 = Icges(4,4) - Icges(5,5);
t200 = Icges(4,1) + Icges(5,1);
t199 = Icges(4,2) + Icges(5,3);
t129 = qJ(2) + pkin(6);
t124 = cos(t129);
t198 = t201 * t124;
t123 = sin(t129);
t197 = t201 * t123;
t196 = Icges(5,4) + Icges(4,5);
t195 = Icges(4,6) - Icges(5,6);
t194 = t199 * t123 - t198;
t193 = t200 * t124 - t197;
t192 = rSges(5,1) + pkin(3);
t191 = rSges(5,3) + qJ(4);
t132 = sin(qJ(1));
t134 = cos(qJ(1));
t190 = t194 * t132 + t195 * t134;
t189 = -t195 * t132 + t194 * t134;
t188 = t193 * t132 - t196 * t134;
t187 = t196 * t132 + t193 * t134;
t186 = -t199 * t124 - t197;
t185 = t200 * t123 + t198;
t184 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t131 = sin(qJ(2));
t133 = cos(qJ(2));
t183 = Icges(3,5) * t133 - Icges(3,6) * t131 - t195 * t123 + t196 * t124;
t182 = t191 * t123 + t192 * t124;
t168 = Icges(3,4) * t131;
t109 = Icges(3,2) * t133 + t168;
t167 = Icges(3,4) * t133;
t112 = Icges(3,1) * t131 + t167;
t120 = -qJD(2) * t134 + V_base(5);
t121 = qJD(2) * t132 + V_base(4);
t125 = V_base(6) + qJD(1);
t150 = -Icges(3,2) * t131 + t167;
t83 = -Icges(3,6) * t134 + t150 * t132;
t84 = Icges(3,6) * t132 + t150 * t134;
t153 = Icges(3,1) * t133 - t168;
t85 = -Icges(3,5) * t134 + t153 * t132;
t86 = Icges(3,5) * t132 + t153 * t134;
t181 = (-t109 * t131 + t112 * t133 + t186 * t123 + t185 * t124) * t125 + (t189 * t123 + t187 * t124 - t131 * t84 + t133 * t86) * t121 + (t190 * t123 + t188 * t124 - t131 * t83 + t133 * t85) * t120;
t180 = (Icges(3,5) * t131 + Icges(3,6) * t133 + t196 * t123 + t195 * t124) * t125 + (t184 * t132 + t183 * t134) * t121 + (t183 * t132 - t184 * t134) * t120;
t176 = pkin(2) * t131;
t175 = t133 * pkin(2);
t173 = -t134 * rSges(5,2) + t182 * t132;
t172 = t132 * rSges(5,2) + t182 * t134;
t171 = t192 * t123 - t191 * t124;
t118 = t132 * pkin(1) - t134 * pkin(5);
t63 = -qJ(3) * t134 + t175 * t132;
t170 = -t118 - t63;
t169 = Icges(2,4) * t132;
t162 = qJD(4) * t123;
t161 = V_base(5) * pkin(4) + V_base(1);
t158 = qJD(3) * t132 + t120 * t176 + t161;
t157 = rSges(3,1) * t133 - rSges(3,2) * t131;
t156 = rSges(4,1) * t124 - rSges(4,2) * t123;
t119 = t134 * pkin(1) + t132 * pkin(5);
t144 = -V_base(4) * pkin(4) + t125 * t119 + V_base(2);
t143 = V_base(4) * t118 - V_base(5) * t119 + V_base(3);
t142 = t121 * t63 + t143;
t64 = qJ(3) * t132 + t175 * t134;
t138 = -qJD(3) * t134 + t125 * t64 + t144;
t127 = Icges(2,4) * t134;
t117 = t134 * rSges(2,1) - t132 * rSges(2,2);
t116 = t132 * rSges(2,1) + t134 * rSges(2,2);
t115 = t131 * rSges(3,1) + t133 * rSges(3,2);
t114 = Icges(2,1) * t134 - t169;
t113 = Icges(2,1) * t132 + t127;
t111 = -Icges(2,2) * t132 + t127;
t110 = Icges(2,2) * t134 + t169;
t105 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t104 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t103 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t100 = t123 * rSges(4,1) + t124 * rSges(4,2);
t90 = t132 * rSges(3,3) + t157 * t134;
t89 = -t134 * rSges(3,3) + t157 * t132;
t80 = t132 * rSges(4,3) + t156 * t134;
t78 = -t134 * rSges(4,3) + t156 * t132;
t62 = V_base(5) * rSges(2,3) - t125 * t116 + t161;
t61 = t125 * t117 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t60 = V_base(4) * t116 - V_base(5) * t117 + V_base(3);
t57 = t120 * t115 + (-t118 - t89) * t125 + t161;
t56 = -t121 * t115 + t125 * t90 + t144;
t55 = -t120 * t90 + t121 * t89 + t143;
t54 = t120 * t100 + (-t78 + t170) * t125 + t158;
t53 = t125 * t80 + (-t100 - t176) * t121 + t138;
t52 = t121 * t78 + (-t64 - t80) * t120 + t142;
t51 = t134 * t162 + t171 * t120 + (t170 - t173) * t125 + t158;
t50 = t132 * t162 + t172 * t125 + (-t171 - t176) * t121 + t138;
t49 = -qJD(4) * t124 + t173 * t121 + (-t64 - t172) * t120 + t142;
t1 = m(1) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(2) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + ((-t132 * t110 + t134 * t113 + Icges(1,4)) * V_base(5) + (-t132 * t111 + t134 * t114 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t134 * t110 + t132 * t113 + Icges(1,2)) * V_base(5) + (t134 * t111 + t132 * t114 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t181 * t132 - t180 * t134) * t120 / 0.2e1 + (t180 * t132 + t181 * t134) * t121 / 0.2e1 + ((t187 * t123 - t189 * t124 + t131 * t86 + t133 * t84) * t121 + (t188 * t123 - t190 * t124 + t131 * t85 + t133 * t83) * t120 + (t133 * t109 + t131 * t112 + t185 * t123 - t186 * t124 + Icges(2,3)) * t125) * t125 / 0.2e1 + V_base(4) * t125 * (Icges(2,5) * t134 - Icges(2,6) * t132) + V_base(5) * t125 * (Icges(2,5) * t132 + Icges(2,6) * t134) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
