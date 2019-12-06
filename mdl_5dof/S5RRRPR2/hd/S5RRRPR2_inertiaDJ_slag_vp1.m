% Calculate time derivative of joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:48
% DurationCPUTime: 1.57s
% Computational Cost: add. (5296->221), mult. (3006->297), div. (0->0), fcn. (2068->10), ass. (0->136)
t120 = sin(qJ(5));
t156 = qJD(5) * t120;
t122 = cos(qJ(5));
t155 = qJD(5) * t122;
t118 = qJD(1) + qJD(2);
t114 = qJD(3) + t118;
t133 = Icges(6,5) * t122 - Icges(6,6) * t120;
t166 = Icges(6,4) * t122;
t134 = -Icges(6,2) * t120 + t166;
t167 = Icges(6,4) * t120;
t135 = Icges(6,1) * t122 - t167;
t96 = Icges(6,2) * t122 + t167;
t97 = Icges(6,1) * t120 + t166;
t136 = t120 * t96 - t122 * t97;
t194 = t133 * qJD(5) + (-t120 * t135 - t122 * t134 + t136) * t114;
t191 = rSges(6,1) * t156 + rSges(6,2) * t155;
t119 = qJ(1) + qJ(2);
t117 = qJ(3) + t119;
t110 = pkin(9) + t117;
t108 = cos(t110);
t107 = sin(t110);
t46 = Icges(6,6) * t108 - t107 * t134;
t48 = Icges(6,5) * t108 - t107 * t135;
t140 = t120 * t46 - t122 * t48;
t190 = t140 * t108;
t95 = Icges(6,5) * t120 + Icges(6,6) * t122;
t189 = -Icges(6,3) * t114 + qJD(5) * t95;
t185 = 2 * m(3);
t184 = 2 * m(4);
t183 = 2 * m(5);
t182 = 2 * m(6);
t179 = -rSges(6,3) - pkin(8);
t121 = sin(qJ(1));
t178 = pkin(1) * t121;
t123 = cos(qJ(1));
t177 = pkin(1) * t123;
t115 = sin(t119);
t176 = pkin(2) * t115;
t116 = cos(t119);
t175 = pkin(2) * t116;
t111 = sin(t117);
t174 = pkin(3) * t111;
t112 = cos(t117);
t173 = pkin(3) * t112;
t159 = t112 * t114;
t160 = t111 * t114;
t60 = rSges(4,1) * t160 + rSges(4,2) * t159;
t172 = rSges(6,1) * t122;
t171 = rSges(6,2) * t120;
t170 = pkin(1) * qJD(1);
t169 = t107 * rSges(6,3);
t104 = t108 * rSges(6,3);
t93 = t107 * t171;
t168 = t104 + t93;
t157 = t116 * t118;
t158 = t115 * t118;
t64 = rSges(3,1) * t158 + rSges(3,2) * t157;
t162 = t107 * t114;
t161 = t108 * t114;
t151 = t107 * t172;
t154 = -t108 * t191 - t114 * t151;
t94 = t108 * t171;
t153 = t107 * t191 + t114 * t94;
t92 = pkin(3) * t160;
t36 = rSges(5,1) * t162 + rSges(5,2) * t161 + t92;
t152 = pkin(2) * t157;
t101 = pkin(2) * t158;
t42 = t101 + t60;
t150 = t123 * t170;
t145 = -pkin(4) - t172;
t80 = -rSges(3,1) * t116 + rSges(3,2) * t115;
t73 = -rSges(4,1) * t112 + rSges(4,2) * t111;
t34 = t101 + t36;
t65 = -rSges(3,1) * t157 + rSges(3,2) * t158;
t61 = -rSges(4,1) * t159 + rSges(4,2) * t160;
t144 = -rSges(5,1) * t108 - t173;
t79 = -rSges(3,1) * t115 - rSges(3,2) * t116;
t72 = -rSges(4,1) * t111 - rSges(4,2) * t112;
t141 = t120 * t48 + t122 * t46;
t47 = Icges(6,6) * t107 + t108 * t134;
t49 = Icges(6,5) * t107 + t108 * t135;
t139 = t120 * t49 + t122 * t47;
t138 = t120 * t47 - t122 * t49;
t63 = t73 - t175;
t55 = rSges(5,2) * t107 + t144;
t132 = (t135 - t96) * t156 + (t134 + t97) * t155;
t131 = (-t138 * qJD(5) + t194 * t107) * t107 / 0.2e1 + (-t140 * qJD(5) + t194 * t108) * t108 / 0.2e1 - (t107 * t136 + t108 * t95 + t141) * t162 / 0.2e1 + (t107 * t95 - t108 * t136 + t139) * t161 / 0.2e1;
t130 = t138 * t107;
t127 = t133 * t114;
t62 = t72 - t176;
t54 = -rSges(5,1) * t107 - rSges(5,2) * t108 - t174;
t43 = t61 - t152;
t37 = rSges(5,2) * t162 + t114 * t144;
t53 = t55 - t175;
t52 = t54 - t176;
t30 = pkin(8) * t108 + t107 * t145 + t168 - t174;
t35 = t37 - t152;
t124 = t107 * t179 + t108 * t145 - t173;
t31 = t94 + t124;
t28 = t30 - t176;
t14 = pkin(4) * t162 + t92 + (t108 * t179 - t93) * t114 - t154;
t29 = t31 - t175;
t12 = t101 + t14;
t15 = t114 * t124 + t153;
t13 = t15 - t152;
t113 = t121 * t170;
t102 = rSges(6,1) * t120 + rSges(6,2) * t122;
t91 = (-t171 + t172) * qJD(5);
t67 = t80 - t177;
t66 = t79 - t178;
t59 = t65 - t150;
t58 = t113 + t64;
t57 = t63 - t177;
t56 = t62 - t178;
t51 = t108 * t172 + t169 - t94;
t50 = -t151 + t168;
t45 = Icges(6,3) * t107 + t108 * t133;
t44 = Icges(6,3) * t108 - t107 * t133;
t41 = t53 - t177;
t40 = t52 - t178;
t39 = t43 - t150;
t38 = t113 + t42;
t33 = t35 - t150;
t32 = t113 + t34;
t27 = t29 - t177;
t26 = t28 - t178;
t19 = t107 * t189 - t108 * t127;
t18 = -t107 * t127 - t108 * t189;
t11 = t13 - t150;
t10 = t113 + t12;
t9 = t107 * t45 - t108 * t138;
t8 = t107 * t44 - t190;
t7 = t108 * t45 + t130;
t6 = t107 * t140 + t108 * t44;
t1 = t108 * t154 - t107 * t153 + ((-t50 + t104) * t108 + (t169 - t51 + (t171 + t172) * t108) * t107) * t114;
t2 = [(t58 * t67 + t59 * t66) * t185 + (t38 * t57 + t39 * t56) * t184 + (t32 * t41 + t33 * t40) * t183 + (t10 * t27 + t11 * t26) * t182 + t132; m(3) * (t58 * t80 + t59 * t79 + t64 * t67 + t65 * t66) + m(4) * (t38 * t63 + t39 * t62 + t42 * t57 + t43 * t56) + m(5) * (t32 * t53 + t33 * t52 + t34 * t41 + t35 * t40) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t132; (t12 * t29 + t13 * t28) * t182 + (t34 * t53 + t35 * t52) * t183 + (t42 * t63 + t43 * t62) * t184 + (t64 * t80 + t65 * t79) * t185 + t132; m(4) * (t38 * t73 + t39 * t72 + t56 * t61 + t57 * t60) + m(5) * (t32 * t55 + t33 * t54 + t36 * t41 + t37 * t40) + m(6) * (t10 * t31 + t11 * t30 + t14 * t27 + t15 * t26) + t132; m(6) * (t12 * t31 + t13 * t30 + t14 * t29 + t15 * t28) + m(5) * (t34 * t55 + t35 * t54 + t36 * t53 + t37 * t52) + m(4) * (t42 * t73 + t43 * t72 + t60 * t63 + t61 * t62) + t132; (t60 * t73 + t61 * t72) * t184 + (t36 * t55 + t37 * t54) * t183 + (t14 * t31 + t15 * t30) * t182 + t132; 0; 0; 0; 0; m(6) * ((t107 * t27 - t108 * t26) * t91 + ((t114 * t27 - t11) * t108 + (t114 * t26 + t10) * t107) * t102) + t131; m(6) * ((t107 * t29 - t108 * t28) * t91 + ((t114 * t29 - t13) * t108 + (t114 * t28 + t12) * t107) * t102) + t131; m(6) * ((t107 * t31 - t108 * t30) * t91 + ((t114 * t31 - t15) * t108 + (t114 * t30 + t14) * t107) * t102) + t131; m(6) * t1; ((-t107 * t50 + t108 * t51) * t1 + (t107 ^ 2 + t108 ^ 2) * t91 * t102) * t182 - (t107 * t7 + t108 * t6) * t162 + t108 * ((t108 * t19 + (t7 + t190) * t114) * t108 + (-t6 * t114 + (t155 * t47 + t156 * t49) * t107 + (t141 * qJD(5) + t114 * t138 + t18) * t108) * t107) + (t107 * t9 + t108 * t8) * t161 + t107 * ((t107 * t18 + (-t8 + t130) * t114) * t107 + (t9 * t114 + (-t155 * t46 - t156 * t48) * t108 + (-t139 * qJD(5) + t114 * t140 + t19) * t107) * t108);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
