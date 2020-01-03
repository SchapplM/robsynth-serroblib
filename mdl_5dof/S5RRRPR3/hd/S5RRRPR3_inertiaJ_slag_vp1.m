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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:08:50
% EndTime: 2020-01-03 12:08:54
% DurationCPUTime: 1.22s
% Computational Cost: add. (2683->207), mult. (1866->292), div. (0->0), fcn. (1632->10), ass. (0->114)
t187 = Icges(4,3) + Icges(5,3);
t121 = qJ(3) + pkin(9);
t112 = sin(t121);
t113 = cos(t121);
t124 = sin(qJ(3));
t126 = cos(qJ(3));
t186 = Icges(4,5) * t126 + Icges(5,5) * t113 - Icges(4,6) * t124 - Icges(5,6) * t112;
t185 = -Icges(4,5) * t124 - Icges(5,5) * t112 - Icges(4,6) * t126 - Icges(5,6) * t113;
t122 = qJ(1) + qJ(2);
t115 = sin(t122);
t116 = cos(t122);
t184 = t115 * t116;
t163 = Icges(5,4) * t112;
t71 = Icges(5,2) * t113 + t163;
t162 = Icges(5,4) * t113;
t72 = Icges(5,1) * t112 + t162;
t165 = Icges(4,4) * t124;
t88 = Icges(4,2) * t126 + t165;
t164 = Icges(4,4) * t126;
t89 = Icges(4,1) * t124 + t164;
t183 = t112 * t71 - t113 * t72 + t124 * t88 - t126 * t89;
t114 = qJ(5) + t121;
t104 = sin(t114);
t105 = cos(t114);
t182 = -rSges(6,1) * t105 + rSges(6,2) * t104;
t181 = rSges(5,1) * t113 - rSges(5,2) * t112;
t180 = -t186 * t115 + t187 * t116;
t179 = t187 * t115 + t186 * t116;
t42 = -t115 * rSges(6,3) + t182 * t116;
t109 = t126 * pkin(3) + pkin(2);
t80 = pkin(4) * t113 + t109;
t178 = t116 * t80 - t42;
t110 = t115 ^ 2;
t111 = t116 ^ 2;
t90 = t124 * rSges(4,1) + t126 * rSges(4,2);
t177 = m(4) * t90;
t69 = t104 * rSges(6,1) + t105 * rSges(6,2);
t176 = m(6) * t69;
t175 = -t115 / 0.2e1;
t174 = -t116 / 0.2e1;
t173 = pkin(3) * t124;
t123 = -qJ(4) - pkin(7);
t120 = -pkin(8) + t123;
t172 = t115 * t80 + t116 * t120;
t171 = t115 * t109 + t116 * t123;
t170 = rSges(4,1) * t126;
t99 = t115 * rSges(5,3);
t161 = Icges(6,4) * t104;
t160 = Icges(6,4) * t105;
t159 = t116 * t124;
t74 = t115 * rSges(3,1) + t116 * rSges(3,2);
t158 = t116 * pkin(2) + t115 * pkin(7);
t157 = t111 + t110;
t136 = -Icges(6,2) * t104 + t160;
t139 = Icges(6,1) * t105 - t161;
t67 = Icges(6,2) * t105 + t161;
t68 = Icges(6,1) * t104 + t160;
t148 = t104 * t67 - t105 * t68;
t66 = Icges(6,5) * t104 + Icges(6,6) * t105;
t156 = (t104 * (-Icges(6,5) * t115 - t139 * t116) + t105 * (-Icges(6,6) * t115 - t136 * t116) - t115 * t66 + t148 * t116) * t175 + (t104 * (-Icges(6,5) * t116 + t139 * t115) + t105 * (-Icges(6,6) * t116 + t136 * t115) - t148 * t115 - t116 * t66) * t174;
t155 = pkin(4) * t112 + t69;
t75 = t116 * rSges(3,1) - t115 * rSges(3,2);
t85 = t116 * t109;
t154 = t115 * t123 - t85;
t153 = (-rSges(4,2) * t124 + t170) * t115;
t152 = t181 * t115;
t133 = Icges(6,5) * t105 - Icges(6,6) * t104;
t35 = -Icges(6,3) * t116 + t133 * t115;
t36 = -Icges(6,3) * t115 - t133 * t116;
t151 = -t116 * (t111 * t35 + t36 * t184) - t115 * (t110 * t36 + t35 * t184);
t141 = Icges(4,1) * t126 - t165;
t140 = Icges(5,1) * t113 - t163;
t138 = -Icges(4,2) * t124 + t164;
t137 = -Icges(5,2) * t112 + t162;
t132 = t181 * t116 + t99;
t131 = rSges(4,2) * t159 - t115 * rSges(4,3) - t116 * t170;
t130 = t104 * t68 + t105 * t67 + t112 * t72 + t113 * t71 + t124 * t89 + t126 * t88 + Icges(3,3);
t129 = -t116 * rSges(6,3) - t182 * t115;
t33 = -t131 + t158;
t25 = -t116 * rSges(5,3) + t152 + t171;
t19 = t129 + t172;
t26 = t132 - t154;
t20 = -t115 * t120 + t178;
t107 = t115 * pkin(2);
t32 = t107 + (-rSges(4,3) - pkin(7)) * t116 + t153;
t128 = t156 + (t112 * (-Icges(5,5) * t115 - t140 * t116) + t113 * (-Icges(5,6) * t115 - t137 * t116) + t124 * (-Icges(4,5) * t115 - t141 * t116) + t126 * (-Icges(4,6) * t115 - t138 * t116) + t183 * t116 + t185 * t115) * t175 + (t112 * (-Icges(5,5) * t116 + t140 * t115) + t113 * (-Icges(5,6) * t116 + t137 * t115) + t124 * (-Icges(4,5) * t116 + t141 * t115) + t126 * (-Icges(4,6) * t116 + t138 * t115) + t185 * t116 - t183 * t115) * t174;
t127 = cos(qJ(1));
t125 = sin(qJ(1));
t119 = t127 * pkin(1);
t117 = t125 * pkin(1);
t96 = pkin(3) * t159;
t92 = t127 * rSges(2,1) - t125 * rSges(2,2);
t91 = t125 * rSges(2,1) + t127 * rSges(2,2);
t73 = t112 * rSges(5,1) + t113 * rSges(5,2);
t62 = t119 + t75;
t61 = t117 + t74;
t50 = t116 * t73 + t96;
t49 = (-t73 - t173) * t115;
t41 = t154 + t158;
t34 = t115 * t129;
t31 = t115 * (t116 * pkin(7) - t107 + t171);
t30 = t119 + t33;
t29 = t117 + t32;
t28 = t155 * t116 + t96;
t27 = (-t155 - t173) * t115;
t24 = t119 + t26;
t23 = t117 + t25;
t18 = -t116 * t131 + t115 * (-t116 * rSges(4,3) + t153);
t17 = t119 + t20;
t16 = t117 + t19;
t13 = -t116 * t42 + t34;
t4 = t31 + t115 * t152 + (t132 - t41 - t99) * t116;
t3 = t31 + t115 * (-t171 + t172) + t34 + (-t41 - t85 + (-t120 + t123) * t115 + t178) * t116;
t1 = [Icges(2,3) + m(2) * (t91 ^ 2 + t92 ^ 2) + m(3) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2) + t130; m(3) * (t74 * t61 + t75 * t62) + m(4) * (t32 * t29 + t33 * t30) + m(5) * (t25 * t23 + t26 * t24) + m(6) * (t19 * t16 + t20 * t17) + t130; m(6) * (t19 ^ 2 + t20 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t74 ^ 2 + t75 ^ 2) + t130; t128 + m(5) * (t50 * t23 + t49 * t24) + m(6) * (t28 * t16 + t27 * t17) + (-t115 * t30 + t116 * t29) * t177; t128 + m(6) * (t28 * t19 + t27 * t20) + m(5) * (t50 * t25 + t49 * t26) + (-t115 * t33 + t116 * t32) * t177; m(6) * (t27 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(5) * (t4 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(4) * (t157 * t90 ^ 2 + t18 ^ 2) + t151 + t179 * t115 * t110 + (t180 * t111 + (t180 * t115 + t179 * t116) * t115) * t116; m(5) * (-t115 * t23 - t116 * t24) + m(6) * (-t115 * t16 - t116 * t17); m(6) * (-t115 * t19 - t116 * t20) + m(5) * (-t115 * t25 - t116 * t26); m(6) * (-t115 * t28 - t116 * t27) + m(5) * (-t115 * t50 - t116 * t49); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t157; (-t115 * t17 + t116 * t16) * t176 + t156; (-t115 * t20 + t116 * t19) * t176 + t156; m(6) * (t13 * t3 + (-t115 * t27 + t116 * t28) * t69) + t151; 0; m(6) * (t157 * t69 ^ 2 + t13 ^ 2) + t151;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
