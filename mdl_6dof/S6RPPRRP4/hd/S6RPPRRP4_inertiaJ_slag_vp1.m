% Calculate joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:10
% EndTime: 2019-03-09 02:05:15
% DurationCPUTime: 1.92s
% Computational Cost: add. (3874->298), mult. (8532->441), div. (0->0), fcn. (10746->8), ass. (0->135)
t186 = rSges(7,1) + pkin(5);
t185 = rSges(7,3) + qJ(6);
t131 = cos(qJ(4));
t127 = sin(qJ(5));
t128 = sin(qJ(4));
t153 = t127 * t128;
t130 = cos(qJ(5));
t88 = Icges(6,3) * t131 + (-Icges(6,5) * t130 + Icges(6,6) * t127) * t128;
t89 = Icges(7,2) * t131 + (-Icges(7,4) * t130 - Icges(7,6) * t127) * t128;
t90 = Icges(6,6) * t131 + (-Icges(6,4) * t130 + Icges(6,2) * t127) * t128;
t184 = t90 * t153 + (t88 + t89) * t131;
t87 = Icges(7,6) * t131 + (-Icges(7,5) * t130 - Icges(7,3) * t127) * t128;
t91 = Icges(7,4) * t131 + (-Icges(7,1) * t130 - Icges(7,5) * t127) * t128;
t92 = Icges(6,5) * t131 + (-Icges(6,1) * t130 + Icges(6,4) * t127) * t128;
t183 = -t127 * t87 + (-t91 - t92) * t130;
t129 = sin(qJ(1));
t132 = cos(qJ(1));
t159 = sin(pkin(9));
t160 = cos(pkin(9));
t108 = -t129 * t159 - t132 * t160;
t109 = -t129 * t160 + t132 * t159;
t152 = t127 * t131;
t78 = t108 * t130 - t109 * t152;
t151 = t130 * t131;
t79 = -t108 * t127 - t109 * t151;
t182 = -t185 * t78 - t186 * t79;
t154 = t109 * t128;
t38 = Icges(7,5) * t79 - Icges(7,6) * t154 + Icges(7,3) * t78;
t42 = Icges(7,4) * t79 - Icges(7,2) * t154 + Icges(7,6) * t78;
t46 = Icges(7,1) * t79 - Icges(7,4) * t154 + Icges(7,5) * t78;
t11 = -t154 * t42 + t78 * t38 + t79 * t46;
t156 = t108 * t128;
t80 = -t108 * t152 - t109 * t130;
t81 = -t108 * t151 + t109 * t127;
t39 = Icges(7,5) * t81 - Icges(7,6) * t156 + Icges(7,3) * t80;
t43 = Icges(7,4) * t81 - Icges(7,2) * t156 + Icges(7,6) * t80;
t47 = Icges(7,1) * t81 - Icges(7,4) * t156 + Icges(7,5) * t80;
t12 = -t154 * t43 + t78 * t39 + t79 * t47;
t40 = Icges(6,5) * t79 - Icges(6,6) * t78 - Icges(6,3) * t154;
t44 = Icges(6,4) * t79 - Icges(6,2) * t78 - Icges(6,6) * t154;
t48 = Icges(6,1) * t79 - Icges(6,4) * t78 - Icges(6,5) * t154;
t13 = -t154 * t40 - t78 * t44 + t79 * t48;
t41 = Icges(6,5) * t81 - Icges(6,6) * t80 - Icges(6,3) * t156;
t45 = Icges(6,4) * t81 - Icges(6,2) * t80 - Icges(6,6) * t156;
t49 = Icges(6,1) * t81 - Icges(6,4) * t80 - Icges(6,5) * t156;
t14 = -t154 * t41 - t78 * t45 + t79 * t49;
t31 = -t154 * t89 + t78 * t87 + t79 * t91;
t32 = -t154 * t88 - t78 * t90 + t79 * t92;
t181 = (t32 + t31) * t131 + ((-t11 - t13) * t109 + (-t12 - t14) * t108) * t128;
t15 = -t156 * t42 + t80 * t38 + t81 * t46;
t16 = -t156 * t43 + t80 * t39 + t81 * t47;
t17 = -t156 * t40 - t80 * t44 + t81 * t48;
t18 = -t156 * t41 - t80 * t45 + t81 * t49;
t33 = -t156 * t89 + t80 * t87 + t81 * t91;
t34 = -t156 * t88 - t80 * t90 + t81 * t92;
t180 = (t34 + t33) * t131 + ((-t15 - t17) * t109 + (-t16 - t18) * t108) * t128;
t179 = -t128 / 0.2e1;
t178 = -t131 / 0.2e1;
t20 = t131 * t42 + (-t127 * t38 - t130 * t46) * t128;
t22 = t131 * t40 + (t127 * t44 - t130 * t48) * t128;
t177 = -t20 - t22;
t21 = t131 * t43 + (-t127 * t39 - t130 * t47) * t128;
t23 = t131 * t41 + (t127 * t45 - t130 * t49) * t128;
t176 = t21 + t23;
t175 = t108 ^ 2;
t174 = t109 ^ 2;
t171 = Icges(5,5) * t179 + Icges(5,6) * t178;
t170 = t131 / 0.2e1;
t115 = -t128 * rSges(5,1) - t131 * rSges(5,2);
t169 = m(5) * t115;
t168 = pkin(4) * t131;
t167 = rSges(7,2) * t154 + t182;
t166 = -rSges(7,2) * t156 + t185 * t80 + t186 * t81;
t165 = (t128 * t183 + t184) * t131;
t163 = t108 * rSges(5,3);
t161 = t131 * rSges(7,2) + (-t185 * t127 - t130 * t186) * t128;
t158 = Icges(5,4) * t128;
t157 = Icges(5,4) * t131;
t155 = t108 * t131;
t150 = -pkin(4) * t155 - pkin(8) * t156;
t149 = t132 * pkin(1) + t129 * qJ(2);
t53 = t81 * rSges(6,1) - t80 * rSges(6,2) - rSges(6,3) * t156;
t148 = t132 * pkin(2) + t149;
t147 = pkin(3) + t168;
t124 = t132 * qJ(2);
t146 = t124 + (-pkin(1) - pkin(2)) * t129;
t145 = -t79 * rSges(6,1) + t78 * rSges(6,2);
t144 = -rSges(5,1) * t131 + rSges(5,2) * t128;
t141 = -Icges(5,1) * t131 + t158;
t140 = Icges(5,2) * t128 - t157;
t139 = -Icges(5,5) * t131 + Icges(5,6) * t128;
t138 = -rSges(5,1) * t155 + rSges(5,2) * t156 + t109 * rSges(5,3);
t137 = t108 * pkin(7) + t146;
t136 = -t22 / 0.2e1 - t20 / 0.2e1 - t32 / 0.2e1 - t31 / 0.2e1;
t135 = -t23 / 0.2e1 - t21 / 0.2e1 - t34 / 0.2e1 - t33 / 0.2e1;
t134 = -t108 * pkin(3) + t109 * pkin(7) + t148;
t133 = t134 + t150;
t118 = -t128 * pkin(4) + t131 * pkin(8);
t117 = t132 * rSges(2,1) - t129 * rSges(2,2);
t116 = -t129 * rSges(2,1) - t132 * rSges(2,2);
t96 = t132 * rSges(3,1) + t129 * rSges(3,3) + t149;
t95 = t132 * rSges(3,3) + t124 + (-rSges(3,1) - pkin(1)) * t129;
t94 = t131 * rSges(6,3) + (-rSges(6,1) * t130 + rSges(6,2) * t127) * t128;
t84 = t108 * t118;
t82 = (-pkin(8) * t128 - t168) * t109;
t72 = -t108 * rSges(4,1) - t109 * rSges(4,2) + t148;
t71 = t109 * rSges(4,1) - t108 * rSges(4,2) + t146;
t70 = t108 * t150;
t65 = Icges(5,3) * t109 + t108 * t139;
t64 = -Icges(5,3) * t108 + t109 * t139;
t63 = -t108 * t94 - t84;
t62 = (-t118 - t94) * t109;
t57 = -t108 * t161 - t84;
t56 = (-t118 - t161) * t109;
t55 = t134 + t138;
t54 = t163 + (pkin(3) - t144) * t109 + t137;
t51 = -rSges(6,3) * t154 - t145;
t37 = t108 * t138 + (t109 * t144 - t163) * t109;
t36 = t131 * t53 + t156 * t94;
t35 = -t131 * t51 - t154 * t94;
t30 = t133 + t53;
t29 = ((rSges(6,3) + pkin(8)) * t128 + t147) * t109 + t137 + t145;
t28 = (-t108 * t51 + t109 * t53) * t128;
t27 = t131 * t166 + t156 * t161;
t26 = t131 * t167 - t154 * t161;
t25 = t133 + t166;
t24 = ((rSges(7,2) + pkin(8)) * t128 + t147) * t109 + t137 + t182;
t19 = t108 * t53 + t70 + (t51 + t82) * t109;
t10 = (t108 * t167 + t109 * t166) * t128;
t9 = t70 + t166 * t108 + (t82 - t167) * t109;
t8 = -t17 * t108 + t18 * t109;
t7 = -t15 * t108 + t16 * t109;
t6 = -t13 * t108 + t14 * t109;
t5 = -t11 * t108 + t12 * t109;
t1 = [-t131 * (-Icges(5,2) * t131 - t158) + Icges(3,2) + Icges(2,3) + Icges(4,3) + (Icges(5,1) * t128 + t157 + t183) * t128 + m(7) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t29 ^ 2 + t30 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2) + m(3) * (t95 ^ 2 + t96 ^ 2) + m(4) * (t71 ^ 2 + t72 ^ 2) + m(2) * (t116 ^ 2 + t117 ^ 2) + t184; m(7) * (t129 * t24 - t132 * t25) + m(6) * (t129 * t29 - t132 * t30) + m(5) * (t129 * t54 - t132 * t55) + m(3) * (t129 * t95 - t132 * t96) + m(4) * (t129 * t71 - t132 * t72); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t129 ^ 2 + t132 ^ 2); 0; 0; m(4) + m(5) + m(6) + m(7); m(7) * (t57 * t24 + t56 * t25) + m(6) * (t63 * t29 + t62 * t30) + (t109 * t171 - t55 * t169 + (Icges(5,6) * t109 + t108 * t140) * t178 + (Icges(5,5) * t109 + t108 * t141) * t179 - t135) * t109 + (-t54 * t169 + t128 * (-Icges(5,5) * t108 + t109 * t141) / 0.2e1 + (-Icges(5,6) * t108 + t109 * t140) * t170 + t108 * t171 + t136) * t108; m(6) * (t63 * t129 - t62 * t132) + m(7) * (t57 * t129 - t56 * t132) + (-t108 * t129 + t109 * t132) * t169; -m(5) * t37 - m(6) * t19 - m(7) * t9; m(7) * (t56 ^ 2 + t57 ^ 2 + t9 ^ 2) + m(6) * (t19 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t37 ^ 2 + (t174 + t175) * t115 ^ 2) + (t174 * t65 + t7 + t8) * t109 + (-t175 * t64 - t5 - t6 + (t108 * t65 - t109 * t64) * t109) * t108; m(7) * (t26 * t24 + t27 * t25) + m(6) * (t35 * t29 + t36 * t30) + (t108 * t135 + t109 * t136) * t128 + t165; m(6) * (t35 * t129 - t36 * t132) + m(7) * (t26 * t129 - t27 * t132); -m(6) * t28 - m(7) * t10; m(7) * (t10 * t9 + t26 * t57 + t27 * t56) + m(6) * (t28 * t19 + t35 * t63 + t36 * t62) + ((-t6 / 0.2e1 - t5 / 0.2e1) * t109 + (-t8 / 0.2e1 - t7 / 0.2e1) * t108) * t128 - t181 * t108 / 0.2e1 + t180 * t109 / 0.2e1 + (t108 * t177 + t109 * t176) * t170; t165 * t131 + m(7) * (t10 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + ((t131 * t177 - t181) * t109 + (-t131 * t176 - t180) * t108) * t128; m(7) * (t80 * t24 + t78 * t25); m(7) * (t80 * t129 - t78 * t132); m(7) * t153; m(7) * (-t153 * t9 + t78 * t56 + t80 * t57); m(7) * (-t10 * t153 + t80 * t26 + t78 * t27); m(7) * (t128 ^ 2 * t127 ^ 2 + t78 ^ 2 + t80 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
