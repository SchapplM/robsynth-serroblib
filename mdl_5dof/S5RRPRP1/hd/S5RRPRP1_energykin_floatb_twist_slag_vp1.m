% Calculate kinetic energy for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:21:55
% EndTime: 2019-12-05 18:21:58
% DurationCPUTime: 2.20s
% Computational Cost: add. (1078->198), mult. (750->248), div. (0->0), fcn. (528->8), ass. (0->103)
t225 = Icges(5,4) + Icges(6,4);
t224 = Icges(5,1) + Icges(6,1);
t223 = Icges(5,2) + Icges(6,2);
t146 = cos(qJ(4));
t222 = t225 * t146;
t144 = sin(qJ(4));
t221 = t225 * t144;
t220 = Icges(5,5) + Icges(6,5);
t219 = Icges(5,6) + Icges(6,6);
t218 = -t223 * t144 + t222;
t217 = t224 * t146 - t221;
t216 = rSges(6,1) + pkin(4);
t142 = qJ(1) + qJ(2);
t135 = pkin(8) + t142;
t131 = sin(t135);
t132 = cos(t135);
t215 = -t218 * t131 + t219 * t132;
t214 = t219 * t131 + t218 * t132;
t213 = -t217 * t131 + t220 * t132;
t212 = t220 * t131 + t217 * t132;
t211 = Icges(5,3) + Icges(6,3);
t210 = t223 * t146 + t221;
t209 = t224 * t144 + t222;
t208 = -t219 * t144 + t220 * t146;
t207 = rSges(6,3) + qJ(5);
t206 = -rSges(6,2) * t144 + t216 * t146;
t106 = qJD(4) * t131 + V_base(6);
t107 = qJD(4) * t132 + V_base(5);
t136 = V_base(4) + qJD(1);
t133 = qJD(2) + t136;
t205 = t106 * (t214 * t144 - t212 * t146) + t107 * (t215 * t144 - t213 * t146) + (t210 * t144 - t209 * t146) * t133;
t202 = (t220 * t144 + t219 * t146) * t133 + (-t208 * t131 + t211 * t132) * t107 + (t211 * t131 + t208 * t132) * t106;
t198 = -pkin(5) - pkin(6);
t145 = sin(qJ(1));
t194 = pkin(1) * t145;
t147 = cos(qJ(1));
t193 = pkin(1) * t147;
t137 = sin(t142);
t192 = pkin(2) * t137;
t138 = cos(t142);
t191 = pkin(2) * t138;
t188 = t206 * t131 - t132 * t207;
t187 = t131 * t207 + t206 * t132;
t186 = Icges(2,4) * t145;
t185 = Icges(2,4) * t147;
t184 = Icges(3,4) * t137;
t183 = Icges(3,4) * t138;
t182 = Icges(4,4) * t131;
t181 = Icges(4,4) * t132;
t176 = -qJ(3) + t198;
t175 = V_base(6) * pkin(5) + V_base(2);
t97 = pkin(3) * t132 + pkin(7) * t131;
t172 = -t97 - t191;
t171 = rSges(6,2) * t146 + t216 * t144;
t170 = V_base(5) * t193 + V_base(6) * t194 + V_base(1);
t169 = rSges(5,1) * t146 - rSges(5,2) * t144;
t167 = -t136 * t194 + V_base(3);
t154 = V_base(5) * t191 + V_base(6) * t192 + qJD(3) + t170;
t153 = V_base(6) * pkin(6) - t136 * t193 + t175;
t150 = V_base(6) * qJ(3) + t153;
t96 = -pkin(3) * t131 + pkin(7) * t132;
t149 = -t96 * V_base(6) + V_base(5) * t97 + t154;
t148 = t133 * t96 + t176 * V_base(5) + t167;
t126 = rSges(2,1) * t147 - t145 * rSges(2,2);
t125 = -t145 * rSges(2,1) - rSges(2,2) * t147;
t124 = rSges(5,1) * t144 + rSges(5,2) * t146;
t122 = Icges(2,1) * t147 - t186;
t121 = -Icges(2,1) * t145 - t185;
t118 = -Icges(2,2) * t145 + t185;
t117 = -Icges(2,2) * t147 - t186;
t110 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t109 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t108 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t105 = rSges(3,1) * t138 - rSges(3,2) * t137;
t104 = -rSges(3,1) * t137 - rSges(3,2) * t138;
t103 = Icges(3,1) * t138 - t184;
t102 = -Icges(3,1) * t137 - t183;
t101 = -Icges(3,2) * t137 + t183;
t100 = -Icges(3,2) * t138 - t184;
t95 = rSges(4,1) * t132 - rSges(4,2) * t131;
t94 = -rSges(4,1) * t131 - rSges(4,2) * t132;
t93 = Icges(4,1) * t132 - t182;
t92 = -Icges(4,1) * t131 - t181;
t91 = -Icges(4,2) * t131 + t181;
t90 = -Icges(4,2) * t132 - t182;
t85 = V_base(6) * rSges(2,3) - t126 * t136 + t175;
t84 = t125 * t136 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t83 = -t125 * V_base(6) + t126 * V_base(5) + V_base(1);
t82 = rSges(5,3) * t131 + t132 * t169;
t80 = rSges(5,3) * t132 - t131 * t169;
t64 = V_base(6) * rSges(3,3) - t133 * t105 + t153;
t63 = t104 * t133 + (-rSges(3,3) + t198) * V_base(5) + t167;
t62 = -t104 * V_base(6) + t105 * V_base(5) + t170;
t61 = V_base(6) * rSges(4,3) + (-t95 - t191) * t133 + t150;
t60 = (t94 - t192) * t133 + (-rSges(4,3) + t176) * V_base(5) + t167;
t59 = -t94 * V_base(6) + t95 * V_base(5) + t154;
t58 = t106 * t124 + (t172 - t82) * t133 + t150;
t57 = -t107 * t124 + (t80 - t192) * t133 + t148;
t56 = -t106 * t80 + t107 * t82 + t149;
t55 = qJD(5) * t132 + t171 * t106 + (t172 - t187) * t133 + t150;
t54 = qJD(5) * t131 - t171 * t107 + (-t188 - t192) * t133 + t148;
t53 = t106 * t188 + t107 * t187 + t149;
t1 = m(1) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(2) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(6) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + (t202 * t131 - t205 * t132) * t106 / 0.2e1 + (t205 * t131 + t202 * t132) * t107 / 0.2e1 + ((t213 * t144 + t215 * t146) * t107 + (t212 * t144 + t214 * t146) * t106 + (t209 * t144 + t210 * t146 + Icges(3,3) + Icges(4,3)) * t133) * t133 / 0.2e1 + ((-t101 * t138 - t103 * t137 - t118 * t147 - t145 * t122 - t131 * t93 - t132 * t91 + Icges(1,6)) * V_base(6) + (-t138 * t100 - t137 * t102 - t147 * t117 - t145 * t121 - t131 * t92 - t132 * t90 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t137 * t101 + t138 * t103 - t145 * t118 + t147 * t122 - t131 * t91 + t132 * t93 + Icges(1,3)) * V_base(6) + (-t100 * t137 + t102 * t138 - t145 * t117 + t121 * t147 - t131 * t90 + t132 * t92 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + V_base(6) * t133 * (Icges(3,5) * t138 + Icges(4,5) * t132 - Icges(3,6) * t137 - Icges(4,6) * t131) + V_base(5) * t133 * (-Icges(3,5) * t137 - Icges(4,5) * t131 - Icges(3,6) * t138 - Icges(4,6) * t132) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t145 - Icges(2,6) * t147) * V_base(5) + (Icges(2,5) * t147 - Icges(2,6) * t145) * V_base(6) + Icges(2,3) * t136 / 0.2e1) * t136;
T = t1;
