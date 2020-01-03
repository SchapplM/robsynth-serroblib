% Calculate joint inertia matrix for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:40
% EndTime: 2019-12-31 21:50:45
% DurationCPUTime: 1.88s
% Computational Cost: add. (3089->201), mult. (2416->289), div. (0->0), fcn. (2153->8), ass. (0->116)
t209 = Icges(5,4) - Icges(6,5);
t208 = Icges(5,1) + Icges(6,1);
t207 = -Icges(5,2) - Icges(6,3);
t120 = qJ(3) + qJ(4);
t115 = sin(t120);
t206 = t209 * t115;
t117 = cos(t120);
t205 = t209 * t117;
t204 = Icges(6,4) + Icges(5,5);
t203 = Icges(5,6) - Icges(6,6);
t202 = t207 * t115 + t205;
t201 = t208 * t117 - t206;
t121 = qJ(1) + qJ(2);
t116 = sin(t121);
t118 = cos(t121);
t200 = t202 * t116 - t203 * t118;
t199 = t203 * t116 + t202 * t118;
t198 = t201 * t116 - t204 * t118;
t197 = t204 * t116 + t201 * t118;
t196 = Icges(6,2) + Icges(5,3);
t195 = -t203 * t115 + t204 * t117;
t194 = rSges(6,1) + pkin(4);
t193 = t207 * t117 - t206;
t192 = t208 * t115 + t205;
t185 = rSges(6,3) + qJ(5);
t191 = -t195 * t116 + t196 * t118;
t190 = t196 * t116 + t195 * t118;
t189 = t199 * t115 - t197 * t117;
t188 = t200 * t115 - t198 * t117;
t186 = t204 * t115 + t203 * t117;
t184 = t193 * t115 + t192 * t117;
t158 = t117 * t118;
t159 = t115 * t118;
t183 = t116 * rSges(6,2) + t194 * t158 + t185 * t159;
t114 = t118 ^ 2;
t182 = t191 * t114 + (t189 * t116 + (-t188 + t190) * t118) * t116;
t113 = t116 ^ 2;
t181 = (t190 * t113 + ((-t189 + t191) * t116 + t188 * t118) * t118) * t116;
t122 = sin(qJ(3));
t124 = cos(qJ(3));
t96 = t122 * rSges(4,1) + t124 * rSges(4,2);
t180 = m(4) * t96;
t80 = t115 * rSges(5,1) + t117 * rSges(5,2);
t179 = m(5) * t80;
t178 = t116 / 0.2e1;
t177 = -t118 / 0.2e1;
t176 = m(6) * t115;
t175 = pkin(3) * t122;
t123 = sin(qJ(1));
t174 = t123 * pkin(1);
t109 = t118 * pkin(7);
t111 = t124 * pkin(3) + pkin(2);
t126 = -pkin(8) - pkin(7);
t152 = t118 * t111 - t116 * t126;
t156 = -t118 * pkin(2) - t116 * pkin(7);
t157 = t118 * t126;
t173 = t116 * (t157 + t109 + (-pkin(2) + t111) * t116) + t118 * (t152 + t156);
t131 = rSges(5,1) * t158 - rSges(5,2) * t159 + t116 * rSges(5,3);
t167 = t116 * t115 * rSges(5,2) + t118 * rSges(5,3);
t169 = rSges(5,1) * t117;
t20 = t116 * (t116 * t169 - t167) + t118 * t131;
t172 = -t115 * t194 + t185 * t117;
t170 = rSges(4,1) * t124;
t168 = rSges(4,2) * t122;
t166 = t118 * rSges(4,3) + t116 * t168;
t165 = Icges(4,4) * t122;
t164 = Icges(4,4) * t124;
t155 = t113 + t114;
t153 = -t80 - t175;
t106 = t118 * rSges(6,2);
t7 = t113 * (pkin(4) * t117 + qJ(5) * t115) + t116 * (-t106 + (rSges(6,1) * t117 + rSges(6,3) * t115) * t116) + t183 * t118;
t82 = t118 * rSges(3,1) - t116 * rSges(3,2);
t151 = t172 - t175;
t81 = -t116 * rSges(3,1) - t118 * rSges(3,2);
t94 = Icges(4,2) * t124 + t165;
t95 = Icges(4,1) * t122 + t164;
t142 = -t122 * t94 + t124 * t95;
t141 = t182 * t118 + t181;
t140 = Icges(4,1) * t124 - t165;
t137 = -Icges(4,2) * t122 + t164;
t134 = Icges(4,5) * t124 - Icges(4,6) * t122;
t130 = t116 * rSges(4,3) + (-t168 + t170) * t118;
t129 = (t197 * t115 + t186 * t116 + t199 * t117 + t184 * t118) * t178 + (t198 * t115 + t184 * t116 + t200 * t117 - t186 * t118) * t177;
t37 = t130 - t156;
t128 = t192 * t115 - t193 * t117 + t122 * t95 + t124 * t94 + Icges(3,3);
t19 = t152 + t183;
t36 = t109 + (-pkin(2) - t170) * t116 + t166;
t29 = t131 + t152;
t28 = -t157 + (-t111 - t169) * t116 + t167;
t93 = Icges(4,5) * t122 + Icges(4,6) * t124;
t127 = t129 + (t116 * t93 + t142 * t118 + t122 * (Icges(4,5) * t116 + t140 * t118) + t124 * (Icges(4,6) * t116 + t137 * t118)) * t178 + (t142 * t116 - t118 * t93 + t122 * (-Icges(4,5) * t118 + t140 * t116) + t124 * (-Icges(4,6) * t118 + t137 * t116)) * t177;
t18 = -t157 + t106 + (-t185 * t115 - t117 * t194 - t111) * t116;
t125 = cos(qJ(1));
t119 = t125 * pkin(1);
t98 = t125 * rSges(2,1) - t123 * rSges(2,2);
t97 = -t123 * rSges(2,1) - t125 * rSges(2,2);
t70 = t119 + t82;
t69 = t81 - t174;
t61 = Icges(4,3) * t116 + t134 * t118;
t60 = -Icges(4,3) * t118 + t134 * t116;
t57 = t153 * t118;
t56 = t153 * t116;
t35 = t172 * t118;
t34 = t172 * t116;
t33 = t119 + t37;
t32 = t36 - t174;
t31 = t151 * t118;
t30 = t151 * t116;
t27 = t119 + t29;
t26 = t28 - t174;
t23 = t116 * (t116 * t170 - t166) + t118 * t130;
t13 = t119 + t19;
t12 = t18 - t174;
t6 = t20 + t173;
t5 = t7 + t173;
t1 = [Icges(2,3) + m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t69 ^ 2 + t70 ^ 2) + m(2) * (t97 ^ 2 + t98 ^ 2) + t128; m(6) * (t18 * t12 + t19 * t13) + m(5) * (t28 * t26 + t29 * t27) + m(4) * (t36 * t32 + t37 * t33) + m(3) * (t81 * t69 + t82 * t70) + t128; m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t28 ^ 2 + t29 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t81 ^ 2 + t82 ^ 2) + t128; t127 + m(6) * (t31 * t12 + t30 * t13) + m(5) * (t57 * t26 + t56 * t27) + (-t116 * t33 - t118 * t32) * t180; t127 + m(6) * (t31 * t18 + t30 * t19) + m(5) * (t57 * t28 + t56 * t29) + (-t116 * t37 - t118 * t36) * t180; m(6) * (t30 ^ 2 + t31 ^ 2 + t5 ^ 2) + m(5) * (t56 ^ 2 + t57 ^ 2 + t6 ^ 2) + t116 * t113 * t61 + m(4) * (t155 * t96 ^ 2 + t23 ^ 2) + t181 + (-t114 * t60 + (-t116 * t60 + t118 * t61) * t116 + t182) * t118; m(6) * (t35 * t12 + t34 * t13) + (-t116 * t27 - t118 * t26) * t179 + t129; m(6) * (t35 * t18 + t34 * t19) + (-t116 * t29 - t118 * t28) * t179 + t129; m(6) * (t34 * t30 + t35 * t31 + t7 * t5) + m(5) * (t20 * t6 + (-t116 * t56 - t118 * t57) * t80) + t141; m(5) * (t155 * t80 ^ 2 + t20 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2 + t7 ^ 2) + t141; (t116 * t13 + t118 * t12) * t176; (t116 * t19 + t118 * t18) * t176; m(6) * (-t117 * t5 + (t116 * t30 + t118 * t31) * t115); m(6) * (-t117 * t7 + (t116 * t34 + t118 * t35) * t115); m(6) * (t155 * t115 ^ 2 + t117 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
