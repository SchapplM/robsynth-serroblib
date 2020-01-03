% Calculate kinetic energy for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:19
% EndTime: 2019-12-31 17:46:20
% DurationCPUTime: 1.48s
% Computational Cost: add. (798->225), mult. (1240->282), div. (0->0), fcn. (1270->8), ass. (0->115)
t220 = Icges(2,4) - Icges(3,5);
t219 = Icges(2,1) + Icges(3,1);
t218 = Icges(3,4) + Icges(2,5);
t217 = Icges(2,2) + Icges(3,3);
t216 = Icges(2,6) - Icges(3,6);
t205 = sin(qJ(1));
t215 = t220 * t205;
t206 = cos(qJ(1));
t214 = t220 * t206;
t213 = -t217 * t206 - t215;
t212 = t217 * t205 - t214;
t211 = t219 * t205 + t214;
t210 = t219 * t206 - t215;
t201 = sin(pkin(7));
t202 = cos(pkin(7));
t118 = -t201 * t205 - t202 * t206;
t119 = t201 * t206 - t202 * t205;
t209 = Icges(4,5) * t119 - Icges(4,6) * t118 + t218 * t205 + t216 * t206;
t208 = Icges(4,5) * t118 + Icges(4,6) * t119 - t216 * t205 + t218 * t206;
t161 = sin(pkin(8));
t204 = pkin(4) * t161;
t162 = cos(pkin(8));
t203 = pkin(4) * t162;
t200 = Icges(4,4) * t118;
t199 = Icges(5,4) * t161;
t198 = Icges(5,4) * t162;
t160 = pkin(8) + qJ(5);
t152 = sin(t160);
t197 = Icges(6,4) * t152;
t153 = cos(t160);
t196 = Icges(6,4) * t153;
t140 = pkin(1) * t205 - qJ(2) * t206;
t194 = V_base(4) * t140 + V_base(3);
t193 = V_base(5) * pkin(5) + V_base(1);
t190 = t206 * pkin(2);
t189 = t205 * pkin(2);
t186 = qJD(2) * t205 + t193;
t185 = -t140 - t189;
t143 = pkin(1) * t206 + qJ(2) * t205;
t184 = -t143 - t190;
t183 = qJD(4) * t119 + t186;
t182 = V_base(4) * t189 - qJD(3) + t194;
t181 = -rSges(5,1) * t162 + rSges(5,2) * t161;
t180 = -rSges(6,1) * t153 + rSges(6,2) * t152;
t179 = -Icges(5,1) * t162 + t199;
t178 = -Icges(6,1) * t153 + t197;
t177 = Icges(5,2) * t161 - t198;
t176 = Icges(6,2) * t152 - t196;
t175 = -Icges(5,5) * t162 + Icges(5,6) * t161;
t174 = -Icges(6,5) * t153 + Icges(6,6) * t152;
t101 = -t119 * pkin(3) - t118 * qJ(4);
t173 = -t101 + t185;
t103 = -t118 * pkin(3) + t119 * qJ(4);
t172 = -t103 + t184;
t171 = V_base(4) * t101 + t182;
t154 = V_base(6) + qJD(1);
t170 = -qJD(2) * t206 + t154 * t143 + V_base(2);
t108 = -qJD(5) * t118 + V_base(5);
t109 = qJD(5) * t119 + V_base(4);
t169 = (-Icges(6,3) * t118 + t119 * t174) * t108 + (Icges(6,3) * t119 + t118 * t174) * t109 + (-Icges(6,5) * t152 - Icges(6,6) * t153) * t154;
t168 = V_base(4) * qJ(3) + t154 * t190 + t170;
t167 = (-Icges(5,5) * t161 - Icges(5,6) * t162) * t154 + (-Icges(5,3) * t118 + t119 * t175) * V_base(5) + (Icges(5,3) * t119 + t118 * t175) * V_base(4);
t166 = -qJD(4) * t118 + t154 * t103 + t168;
t115 = -Icges(6,2) * t153 - t197;
t116 = -Icges(6,1) * t152 - t196;
t76 = -Icges(6,6) * t118 + t119 * t176;
t77 = Icges(6,6) * t119 + t118 * t176;
t78 = -Icges(6,5) * t118 + t119 * t178;
t79 = Icges(6,5) * t119 + t118 * t178;
t165 = (t152 * t77 - t153 * t79) * t109 + (t152 * t76 - t153 * t78) * t108 + (t115 * t152 - t116 * t153) * t154;
t125 = -Icges(5,2) * t162 - t199;
t126 = -Icges(5,1) * t161 - t198;
t85 = -Icges(5,6) * t118 + t119 * t177;
t86 = Icges(5,6) * t119 + t118 * t177;
t87 = -Icges(5,5) * t118 + t119 * t179;
t88 = Icges(5,5) * t119 + t118 * t179;
t164 = (t161 * t86 - t162 * t88) * V_base(4) + (t161 * t85 - t162 * t87) * V_base(5) + (t125 * t161 - t126 * t162) * t154;
t145 = rSges(2,1) * t206 - rSges(2,2) * t205;
t144 = rSges(3,1) * t206 + rSges(3,3) * t205;
t142 = rSges(2,1) * t205 + rSges(2,2) * t206;
t141 = rSges(3,1) * t205 - rSges(3,3) * t206;
t127 = -rSges(5,1) * t161 - rSges(5,2) * t162;
t123 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t122 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t121 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t117 = -rSges(6,1) * t152 - rSges(6,2) * t153;
t112 = Icges(4,4) * t119;
t107 = V_base(5) * rSges(2,3) - t142 * t154 + t193;
t106 = t145 * t154 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t105 = t142 * V_base(4) - t145 * V_base(5) + V_base(3);
t104 = -rSges(4,1) * t118 - rSges(4,2) * t119;
t102 = -rSges(4,1) * t119 + rSges(4,2) * t118;
t100 = -Icges(4,1) * t118 - t112;
t99 = -Icges(4,1) * t119 + t200;
t98 = -Icges(4,2) * t119 - t200;
t97 = Icges(4,2) * t118 - t112;
t92 = V_base(5) * rSges(3,2) + (-t140 - t141) * t154 + t186;
t91 = t154 * t144 + (-pkin(5) - rSges(3,2)) * V_base(4) + t170;
t90 = rSges(5,3) * t119 + t118 * t181;
t89 = -rSges(5,3) * t118 + t119 * t181;
t82 = t141 * V_base(4) + (-t143 - t144) * V_base(5) + t194;
t81 = rSges(6,3) * t119 + t118 * t180;
t80 = -rSges(6,3) * t118 + t119 * t180;
t73 = pkin(6) * t119 - t118 * t203;
t72 = -pkin(6) * t118 - t119 * t203;
t71 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t102 + t185) * t154 + t186;
t70 = t154 * t104 + (rSges(4,3) - pkin(5)) * V_base(4) + t168;
t69 = V_base(4) * t102 + (-t104 + t184) * V_base(5) + t182;
t68 = (-qJ(3) + t127) * V_base(5) + (-t89 + t173) * t154 + t183;
t67 = t154 * t90 + (-pkin(5) - t127) * V_base(4) + t166;
t66 = V_base(4) * t89 + (-t90 + t172) * V_base(5) + t171;
t65 = t108 * t117 + (-qJ(3) - t204) * V_base(5) + (-t72 - t80 + t173) * t154 + t183;
t64 = -t109 * t117 + (-pkin(5) + t204) * V_base(4) + (t73 + t81) * t154 + t166;
t63 = -t108 * t81 + t109 * t80 + V_base(4) * t72 + (-t73 + t172) * V_base(5) + t171;
t1 = m(1) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t82 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t109 * (t165 * t118 + t169 * t119) / 0.2e1 + t108 * (-t169 * t118 + t165 * t119) / 0.2e1 + ((-t152 * t79 - t153 * t77) * t109 + (-t152 * t78 - t153 * t76) * t108 + (-t161 * t87 - t162 * t85 + t209) * V_base(5) + (-t161 * t88 - t162 * t86 + t208) * V_base(4) + (-t153 * t115 - t152 * t116 - t162 * t125 - t161 * t126 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t154) * t154 / 0.2e1 + (t118 * t164 + t119 * t167 + t208 * t154 + (-t118 * t99 - t119 * t97 + t213 * t205 + t211 * t206 + Icges(1,4)) * V_base(5) + (-t118 * t100 - t119 * t98 + t212 * t205 + t210 * t206 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (-t118 * t167 + t119 * t164 + t209 * t154 + (t118 * t97 - t119 * t99 + t211 * t205 - t213 * t206 + Icges(1,2)) * V_base(5) + (-t100 * t119 + t118 * t98 + t210 * t205 - t212 * t206 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
