% Calculate joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:26
% EndTime: 2022-01-20 11:42:27
% DurationCPUTime: 1.04s
% Computational Cost: add. (2683->207), mult. (1866->293), div. (0->0), fcn. (1632->10), ass. (0->112)
t185 = Icges(4,3) + Icges(5,3);
t117 = qJ(3) + pkin(9);
t109 = sin(t117);
t110 = cos(t117);
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t184 = Icges(4,5) * t122 + Icges(5,5) * t110 - Icges(4,6) * t120 - Icges(5,6) * t109;
t183 = Icges(4,5) * t120 + Icges(5,5) * t109 + Icges(4,6) * t122 + Icges(5,6) * t110;
t118 = qJ(1) + qJ(2);
t112 = sin(t118);
t113 = cos(t118);
t182 = t112 * t113;
t159 = Icges(5,4) * t109;
t70 = Icges(5,2) * t110 + t159;
t158 = Icges(5,4) * t110;
t71 = Icges(5,1) * t109 + t158;
t161 = Icges(4,4) * t120;
t86 = Icges(4,2) * t122 + t161;
t160 = Icges(4,4) * t122;
t87 = Icges(4,1) * t120 + t160;
t181 = -t109 * t70 + t110 * t71 - t120 * t86 + t122 * t87;
t180 = -t112 * t184 + t113 * t185;
t179 = t112 * t185 + t113 * t184;
t107 = t112 ^ 2;
t108 = t113 ^ 2;
t88 = rSges(4,1) * t120 + rSges(4,2) * t122;
t178 = m(4) * t88;
t111 = qJ(5) + t117;
t101 = sin(t111);
t102 = cos(t111);
t68 = rSges(6,1) * t101 + rSges(6,2) * t102;
t177 = m(6) * t68;
t176 = t112 / 0.2e1;
t175 = -t113 / 0.2e1;
t174 = pkin(3) * t120;
t121 = sin(qJ(1));
t173 = t121 * pkin(1);
t119 = -qJ(4) - pkin(7);
t106 = t122 * pkin(3) + pkin(2);
t104 = t113 * pkin(7);
t83 = t113 * t106;
t148 = -t112 * t119 + t83;
t153 = -pkin(2) * t113 - pkin(7) * t112;
t154 = t113 * t119;
t172 = t112 * (t154 + t104 + (-pkin(2) + t106) * t112) + t113 * (t148 + t153);
t162 = rSges(6,2) * t101;
t165 = rSges(6,1) * t102;
t127 = t112 * rSges(6,3) + (-t162 + t165) * t113;
t170 = rSges(6,3) * t113 + t112 * t162;
t13 = t112 * (t112 * t165 - t170) + t113 * t127;
t116 = pkin(8) - t119;
t79 = pkin(4) * t110 + t106;
t171 = t112 * t116 + t113 * t79;
t163 = rSges(5,2) * t109;
t169 = rSges(5,3) * t113 + t112 * t163;
t164 = rSges(4,2) * t120;
t168 = rSges(4,3) * t113 + t112 * t164;
t167 = rSges(4,1) * t122;
t166 = rSges(5,1) * t110;
t157 = Icges(6,4) * t101;
t156 = Icges(6,4) * t102;
t155 = t113 * t116;
t152 = t107 + t108;
t133 = -Icges(6,2) * t101 + t156;
t136 = Icges(6,1) * t102 - t157;
t66 = Icges(6,2) * t102 + t157;
t67 = Icges(6,1) * t101 + t156;
t145 = -t101 * t66 + t102 * t67;
t65 = Icges(6,5) * t101 + Icges(6,6) * t102;
t151 = (t101 * (Icges(6,5) * t112 + t113 * t136) + t102 * (Icges(6,6) * t112 + t113 * t133) + t112 * t65 + t145 * t113) * t176 + (t101 * (-Icges(6,5) * t113 + t112 * t136) + t102 * (-Icges(6,6) * t113 + t112 * t133) + t145 * t112 - t113 * t65) * t175;
t130 = Icges(6,5) * t102 - Icges(6,6) * t101;
t37 = -Icges(6,3) * t113 + t112 * t130;
t38 = Icges(6,3) * t112 + t113 * t130;
t150 = -t113 * (t108 * t37 - t182 * t38) + t112 * (t107 * t38 - t182 * t37);
t149 = -rSges(5,1) * t109 - rSges(5,2) * t110 - t174;
t74 = rSges(3,1) * t113 - rSges(3,2) * t112;
t73 = -rSges(3,1) * t112 - rSges(3,2) * t113;
t138 = Icges(4,1) * t122 - t161;
t137 = Icges(5,1) * t110 - t159;
t135 = -Icges(4,2) * t120 + t160;
t134 = -Icges(5,2) * t109 + t158;
t129 = t112 * rSges(4,3) + (-t164 + t167) * t113;
t128 = t112 * rSges(5,3) + (-t163 + t166) * t113;
t126 = -pkin(4) * t109 - t174 - t68;
t125 = t101 * t67 + t102 * t66 + t109 * t71 + t110 * t70 + t120 * t87 + t122 * t86 + Icges(3,3);
t20 = t127 + t171;
t34 = t129 - t153;
t33 = t104 + (-pkin(2) - t167) * t112 + t168;
t26 = t128 + t148;
t124 = t151 + (t109 * (Icges(5,5) * t112 + t113 * t137) + t110 * (Icges(5,6) * t112 + t113 * t134) + t120 * (Icges(4,5) * t112 + t113 * t138) + t122 * (Icges(4,6) * t112 + t113 * t135) + t181 * t113 + t183 * t112) * t176 + (t109 * (-Icges(5,5) * t113 + t112 * t137) + t110 * (-Icges(5,6) * t113 + t112 * t134) + t120 * (-Icges(4,5) * t113 + t112 * t138) + t122 * (-Icges(4,6) * t113 + t112 * t135) - t183 * t113 + t181 * t112) * t175;
t19 = t155 + (-t79 - t165) * t112 + t170;
t25 = -t154 + (-t106 - t166) * t112 + t169;
t123 = cos(qJ(1));
t115 = t123 * pkin(1);
t90 = rSges(2,1) * t123 - rSges(2,2) * t121;
t89 = -rSges(2,1) * t121 - rSges(2,2) * t123;
t62 = t115 + t74;
t61 = t73 - t173;
t50 = t149 * t113;
t49 = t149 * t112;
t30 = t115 + t34;
t29 = t33 - t173;
t28 = t126 * t113;
t27 = t126 * t112;
t24 = t115 + t26;
t23 = t25 - t173;
t18 = t112 * (t112 * t167 - t168) + t113 * t129;
t17 = t115 + t20;
t16 = t19 - t173;
t4 = t112 * (t112 * t166 - t169) + t113 * t128 + t172;
t3 = t113 * (-t83 + t171) + (-t155 + (-t106 + t79) * t112) * t112 + t13 + t172;
t1 = [Icges(2,3) + m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t61 ^ 2 + t62 ^ 2) + m(2) * (t89 ^ 2 + t90 ^ 2) + t125; m(6) * (t16 * t19 + t17 * t20) + m(5) * (t23 * t25 + t24 * t26) + m(4) * (t29 * t33 + t30 * t34) + m(3) * (t61 * t73 + t62 * t74) + t125; m(6) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t73 ^ 2 + t74 ^ 2) + t125; t124 + m(6) * (t16 * t28 + t17 * t27) + m(5) * (t23 * t50 + t24 * t49) + (-t112 * t30 - t113 * t29) * t178; t124 + m(6) * (t19 * t28 + t20 * t27) + m(5) * (t25 * t50 + t26 * t49) + (-t112 * t34 - t113 * t33) * t178; m(6) * (t27 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(4) * (t152 * t88 ^ 2 + t18 ^ 2) + t150 + t179 * t112 * t107 + (t180 * t108 + (t112 * t180 + t113 * t179) * t112) * t113; m(6) * (t112 * t16 - t113 * t17) + m(5) * (t112 * t23 - t113 * t24); m(6) * (t112 * t19 - t113 * t20) + m(5) * (t112 * t25 - t113 * t26); m(6) * (t112 * t28 - t113 * t27) + m(5) * (t112 * t50 - t113 * t49); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t152; (-t112 * t17 - t113 * t16) * t177 + t151; (-t112 * t20 - t113 * t19) * t177 + t151; m(6) * (t13 * t3 + (-t112 * t27 - t113 * t28) * t68) + t150; 0; m(6) * (t152 * t68 ^ 2 + t13 ^ 2) + t150;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
