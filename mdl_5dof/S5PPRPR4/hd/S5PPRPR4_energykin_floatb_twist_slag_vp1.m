% Calculate kinetic energy for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:14
% EndTime: 2019-12-31 17:32:16
% DurationCPUTime: 1.73s
% Computational Cost: add. (774->224), mult. (1240->287), div. (0->0), fcn. (1270->8), ass. (0->115)
t221 = Icges(2,4) - Icges(3,5);
t220 = Icges(2,1) + Icges(3,1);
t219 = Icges(3,4) + Icges(2,5);
t218 = Icges(2,2) + Icges(3,3);
t217 = Icges(2,6) - Icges(3,6);
t202 = sin(pkin(7));
t216 = t221 * t202;
t203 = cos(pkin(7));
t215 = t221 * t203;
t214 = -t218 * t203 - t216;
t213 = t218 * t202 - t215;
t212 = t220 * t202 + t215;
t211 = t220 * t203 - t216;
t207 = cos(qJ(3));
t206 = sin(qJ(3));
t161 = sin(pkin(8));
t205 = pkin(4) * t161;
t162 = cos(pkin(8));
t204 = pkin(4) * t162;
t119 = -t202 * t206 - t203 * t207;
t201 = Icges(4,4) * t119;
t200 = Icges(5,4) * t161;
t199 = Icges(5,4) * t162;
t160 = pkin(8) + qJ(5);
t153 = sin(t160);
t198 = Icges(6,4) * t153;
t154 = cos(t160);
t197 = Icges(6,4) * t154;
t195 = V_base(5) * qJ(1) + V_base(1);
t191 = qJD(1) + V_base(3);
t190 = t203 * pkin(2);
t189 = t202 * pkin(2);
t186 = qJD(2) * t202 + t195;
t140 = pkin(1) * t202 - qJ(2) * t203;
t185 = V_base(4) * t140 + t191;
t143 = pkin(1) * t203 + qJ(2) * t202;
t184 = -t143 - t190;
t183 = V_base(4) * t189 + t185;
t182 = -rSges(5,1) * t162 + rSges(5,2) * t161;
t181 = -rSges(6,1) * t154 + rSges(6,2) * t153;
t180 = -Icges(5,1) * t162 + t200;
t179 = -Icges(6,1) * t154 + t198;
t178 = Icges(5,2) * t161 - t199;
t177 = Icges(6,2) * t153 - t197;
t176 = -Icges(5,5) * t162 + Icges(5,6) * t161;
t175 = -Icges(6,5) * t154 + Icges(6,6) * t153;
t120 = -t202 * t207 + t203 * t206;
t102 = -t120 * pkin(3) - t119 * qJ(4);
t174 = V_base(4) * t102 + t183;
t104 = -t119 * pkin(3) + t120 * qJ(4);
t173 = -t104 + t184;
t172 = -qJD(2) * t203 + V_base(6) * t143 + V_base(2);
t108 = -qJD(5) * t119 + V_base(5);
t109 = qJD(5) * t120 + V_base(4);
t157 = V_base(6) - qJD(3);
t171 = t108 * (-Icges(6,3) * t119 + t120 * t175) + t109 * (Icges(6,3) * t120 + t119 * t175) + (-Icges(6,5) * t153 - Icges(6,6) * t154) * t157;
t170 = V_base(4) * pkin(5) + V_base(6) * t190 + t172;
t169 = (-Icges(5,5) * t161 - Icges(5,6) * t162) * t157 + (-Icges(5,3) * t119 + t120 * t176) * V_base(5) + (Icges(5,3) * t120 + t119 * t176) * V_base(4);
t168 = -qJD(4) * t119 + t157 * t104 + t170;
t167 = (-t189 - t140) * V_base(6) + t186;
t166 = qJD(4) * t120 + t167;
t114 = -Icges(6,2) * t154 - t198;
t115 = -Icges(6,1) * t153 - t197;
t77 = -Icges(6,6) * t119 + t120 * t177;
t78 = Icges(6,6) * t120 + t119 * t177;
t79 = -Icges(6,5) * t119 + t120 * t179;
t80 = Icges(6,5) * t120 + t119 * t179;
t165 = (t153 * t78 - t154 * t80) * t109 + (t153 * t77 - t154 * t79) * t108 + (t114 * t153 - t115 * t154) * t157;
t129 = -Icges(5,2) * t162 - t200;
t134 = -Icges(5,1) * t161 - t199;
t85 = -Icges(5,6) * t119 + t120 * t178;
t86 = Icges(5,6) * t120 + t119 * t178;
t87 = -Icges(5,5) * t119 + t120 * t180;
t88 = Icges(5,5) * t120 + t119 * t180;
t164 = (t161 * t86 - t162 * t88) * V_base(4) + (t161 * t85 - t162 * t87) * V_base(5) + (t129 * t161 - t134 * t162) * t157;
t145 = rSges(2,1) * t203 - rSges(2,2) * t202;
t144 = rSges(3,1) * t203 + rSges(3,3) * t202;
t142 = rSges(2,1) * t202 + rSges(2,2) * t203;
t141 = rSges(3,1) * t202 - rSges(3,3) * t203;
t139 = -rSges(5,1) * t161 - rSges(5,2) * t162;
t123 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t122 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t121 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t116 = -rSges(6,1) * t153 - rSges(6,2) * t154;
t112 = Icges(4,4) * t120;
t107 = V_base(5) * rSges(2,3) - t142 * V_base(6) + t195;
t106 = t145 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t105 = -rSges(4,1) * t119 - rSges(4,2) * t120;
t103 = -rSges(4,1) * t120 + rSges(4,2) * t119;
t101 = t142 * V_base(4) - t145 * V_base(5) + t191;
t100 = -Icges(4,1) * t119 - t112;
t99 = -Icges(4,1) * t120 + t201;
t98 = -Icges(4,2) * t120 - t201;
t97 = Icges(4,2) * t119 - t112;
t96 = Icges(4,5) * t119 + Icges(4,6) * t120;
t95 = Icges(4,5) * t120 - Icges(4,6) * t119;
t92 = V_base(5) * rSges(3,2) + (-t140 - t141) * V_base(6) + t186;
t91 = V_base(6) * t144 + (-qJ(1) - rSges(3,2)) * V_base(4) + t172;
t90 = rSges(5,3) * t120 + t119 * t182;
t89 = -rSges(5,3) * t119 + t120 * t182;
t82 = rSges(6,3) * t120 + t119 * t181;
t81 = -rSges(6,3) * t119 + t120 * t181;
t74 = t141 * V_base(4) + (-t143 - t144) * V_base(5) + t185;
t73 = pkin(6) * t120 - t119 * t204;
t72 = -pkin(6) * t119 - t120 * t204;
t71 = -t157 * t103 + (-pkin(5) - rSges(4,3)) * V_base(5) + t167;
t70 = t157 * t105 + (rSges(4,3) - qJ(1)) * V_base(4) + t170;
t69 = V_base(4) * t103 + (-t105 + t184) * V_base(5) + t183;
t68 = (-pkin(5) + t139) * V_base(5) + (-t102 - t89) * t157 + t166;
t67 = t157 * t90 + (-qJ(1) - t139) * V_base(4) + t168;
t66 = V_base(4) * t89 + (-t90 + t173) * V_base(5) + t174;
t65 = t108 * t116 + (-pkin(5) - t205) * V_base(5) + (-t102 - t72 - t81) * t157 + t166;
t64 = -t109 * t116 + (-qJ(1) + t205) * V_base(4) + (t73 + t82) * t157 + t168;
t63 = -t108 * t82 + t109 * t81 + V_base(4) * t72 + (-t73 + t173) * V_base(5) + t174;
t1 = m(1) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t109 * (t119 * t165 + t120 * t171) / 0.2e1 + t108 * (-t119 * t171 + t120 * t165) / 0.2e1 + ((-t153 * t80 - t154 * t78) * t109 + (-t153 * t79 - t154 * t77) * t108 + (-t161 * t87 - t162 * t85 + t95) * V_base(5) + (-t161 * t88 - t162 * t86 + t96) * V_base(4) + (-t154 * t114 - t153 * t115 - t162 * t129 - t161 * t134 + Icges(4,3)) * t157) * t157 / 0.2e1 + (t119 * t164 + t120 * t169 + t96 * t157 + (-t119 * t99 - t120 * t97 + t202 * t214 + t212 * t203 + Icges(1,4)) * V_base(5) + (-t119 * t100 - t120 * t98 + t213 * t202 + t211 * t203 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (-t119 * t169 + t120 * t164 + t95 * t157 + (t119 * t97 - t120 * t99 + t212 * t202 - t214 * t203 + Icges(1,2)) * V_base(5) + (-t100 * t120 + t119 * t98 + t202 * t211 - t203 * t213 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + (t217 * V_base(5) + t219 * V_base(4)) * t203 + (-t217 * V_base(4) + t219 * V_base(5)) * t202) * V_base(6);
T = t1;
