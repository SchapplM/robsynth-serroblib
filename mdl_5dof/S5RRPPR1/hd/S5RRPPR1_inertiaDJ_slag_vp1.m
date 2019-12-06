% Calculate time derivative of joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:16
% DurationCPUTime: 1.63s
% Computational Cost: add. (3820->200), mult. (2408->275), div. (0->0), fcn. (1712->10), ass. (0->125)
t110 = pkin(9) + qJ(5);
t105 = sin(t110);
t151 = qJD(5) * t105;
t106 = cos(t110);
t150 = qJD(5) * t106;
t111 = qJD(1) + qJD(2);
t128 = Icges(6,5) * t106 - Icges(6,6) * t105;
t160 = Icges(6,4) * t106;
t129 = -Icges(6,2) * t105 + t160;
t161 = Icges(6,4) * t105;
t130 = Icges(6,1) * t106 - t161;
t69 = Icges(6,2) * t106 + t161;
t70 = Icges(6,1) * t105 + t160;
t131 = t105 * t69 - t106 * t70;
t194 = t128 * qJD(5) + (-t105 * t130 - t106 * t129 + t131) * t111;
t112 = qJ(1) + qJ(2);
t107 = pkin(8) + t112;
t101 = sin(t107);
t102 = cos(t107);
t168 = rSges(5,2) * sin(pkin(9));
t188 = rSges(5,3) + qJ(4);
t193 = -t101 * t168 - t188 * t102;
t114 = cos(pkin(9));
t190 = rSges(5,1) * t114 + pkin(3);
t109 = cos(t112);
t172 = pkin(2) * t109;
t35 = -t188 * t101 - t172 + (-t190 + t168) * t102;
t189 = t111 * t35;
t187 = rSges(6,1) * t151 + rSges(6,2) * t150;
t40 = Icges(6,6) * t102 - t129 * t101;
t42 = Icges(6,5) * t102 - t130 * t101;
t134 = t105 * t40 - t106 * t42;
t186 = t134 * t102;
t68 = Icges(6,5) * t105 + Icges(6,6) * t106;
t185 = -Icges(6,3) * t111 + t68 * qJD(5);
t181 = 2 * m(3);
t180 = 2 * m(4);
t179 = 2 * m(5);
t178 = 2 * m(6);
t116 = sin(qJ(1));
t175 = pkin(1) * t116;
t117 = cos(qJ(1));
t174 = pkin(1) * t117;
t108 = sin(t112);
t173 = pkin(2) * t108;
t167 = rSges(6,2) * t105;
t98 = t102 * rSges(6,3);
t171 = t101 * t167 + t98;
t153 = t109 * t111;
t154 = t108 * t111;
t57 = rSges(3,1) * t154 + rSges(3,2) * t153;
t169 = rSges(6,1) * t106;
t166 = pkin(1) * qJD(1);
t165 = t101 * rSges(6,3);
t156 = t101 * t111;
t155 = t102 * t111;
t115 = -pkin(7) - qJ(4);
t152 = t111 * t115;
t147 = t101 * t169;
t149 = -t187 * t102 - t111 * t147;
t83 = t102 * t167;
t148 = t187 * t101 + t111 * t83;
t93 = pkin(2) * t154;
t46 = rSges(4,1) * t156 + rSges(4,2) * t155 + t93;
t146 = t117 * t166;
t104 = t116 * t166;
t103 = pkin(4) * t114 + pkin(3);
t12 = -rSges(6,3) * t155 + t103 * t156 + t102 * t152 + t93 + (-t111 * t167 - qJD(4)) * t101 - t149;
t10 = t104 + t12;
t136 = -t103 - t169;
t30 = t136 * t101 - t102 * t115 + t171 - t173;
t28 = t30 - t175;
t140 = t111 * t28 + t10;
t119 = t136 * t102 - t165 - t172;
t95 = qJD(4) * t102;
t13 = t101 * t152 + t119 * t111 + t148 + t95;
t11 = t13 - t146;
t31 = t101 * t115 + t119 + t83;
t29 = t31 - t174;
t139 = t111 * t29 - t11;
t138 = t111 * t30 + t12;
t137 = t111 * t31 - t13;
t79 = -t109 * rSges(3,1) + t108 * rSges(3,2);
t58 = -rSges(3,1) * t153 + rSges(3,2) * t154;
t135 = -t102 * rSges(4,1) - t172;
t78 = -rSges(3,1) * t108 - rSges(3,2) * t109;
t41 = Icges(6,6) * t101 + t129 * t102;
t43 = Icges(6,5) * t101 + t130 * t102;
t133 = t105 * t41 - t106 * t43;
t55 = t101 * rSges(4,2) + t135;
t127 = (t130 - t69) * t151 + (t129 + t70) * t150;
t126 = (-t133 * qJD(5) + t194 * t101) * t101 / 0.2e1 + (-t134 * qJD(5) + t194 * t102) * t102 / 0.2e1 - (t131 * t101 + t102 * t68 + t105 * t42 + t106 * t40) * t156 / 0.2e1 + (t101 * t68 - t131 * t102 + t105 * t43 + t106 * t41) * t155 / 0.2e1;
t125 = t133 * t101;
t122 = t128 * t111;
t54 = -rSges(4,1) * t101 - rSges(4,2) * t102 - t173;
t47 = rSges(4,2) * t156 + t135 * t111;
t34 = -t101 * t190 - t173 - t193;
t20 = -qJD(4) * t101 + t193 * t111 + t190 * t156 + t93;
t21 = t95 + t189;
t77 = rSges(6,1) * t105 + rSges(6,2) * t106;
t64 = (-t167 + t169) * qJD(5);
t60 = t79 - t174;
t59 = t78 - t175;
t51 = t58 - t146;
t50 = t104 + t57;
t49 = t55 - t174;
t48 = t54 - t175;
t45 = t102 * t169 + t165 - t83;
t44 = -t147 + t171;
t39 = Icges(6,3) * t101 + t128 * t102;
t38 = Icges(6,3) * t102 - t128 * t101;
t37 = t47 - t146;
t36 = t104 + t46;
t33 = t35 - t174;
t32 = t34 - t175;
t23 = t185 * t101 - t102 * t122;
t22 = -t101 * t122 - t185 * t102;
t17 = t21 - t146;
t16 = t104 + t20;
t9 = t101 * t39 - t133 * t102;
t8 = t101 * t38 - t186;
t7 = t102 * t39 + t125;
t6 = t134 * t101 + t102 * t38;
t3 = t102 * t149 - t101 * t148 + ((-t44 + t98) * t102 + (t165 - t45 + (t167 + t169) * t102) * t101) * t111;
t1 = [(t50 * t60 + t51 * t59) * t181 + (t36 * t49 + t37 * t48) * t180 + (t16 * t33 + t17 * t32) * t179 + (t10 * t29 + t11 * t28) * t178 + t127; m(3) * (t50 * t79 + t51 * t78 + t57 * t60 + t58 * t59) + m(4) * (t36 * t55 + t37 * t54 + t46 * t49 + t47 * t48) + m(5) * (t16 * t35 + t17 * t34 + t20 * t33 + t21 * t32) + m(6) * (t10 * t31 + t11 * t30 + t12 * t29 + t13 * t28) + t127; (t12 * t31 + t13 * t30) * t178 + (t46 * t55 + t47 * t54) * t180 + (t20 * t35 + t21 * t34) * t179 + (t57 * t79 + t58 * t78) * t181 + t127; 0; 0; 0; m(5) * ((t111 * t32 + t16) * t102 + (-t111 * t33 + t17) * t101) + m(6) * (-t139 * t101 + t140 * t102); m(6) * (-t137 * t101 + t138 * t102) + m(5) * ((t111 * t34 + t20) * t102 + (t21 - t189) * t101); 0; 0; m(6) * ((t101 * t29 - t102 * t28) * t64 + (t140 * t101 + t139 * t102) * t77) + t126; m(6) * ((t101 * t31 - t102 * t30) * t64 + (t138 * t101 + t137 * t102) * t77) + t126; m(6) * t3; 0; ((-t101 * t44 + t102 * t45) * t3 + (t101 ^ 2 + t102 ^ 2) * t77 * t64) * t178 - (t101 * t7 + t6 * t102) * t156 + t102 * ((t102 * t23 + (t7 + t186) * t111) * t102 + (-t6 * t111 + (t41 * t150 + t43 * t151) * t101 + (t22 + (qJD(5) * t40 - t111 * t43) * t106 + (qJD(5) * t42 + t111 * t41) * t105) * t102) * t101) + (t9 * t101 + t102 * t8) * t155 + t101 * ((t101 * t22 + (-t8 + t125) * t111) * t101 + (t9 * t111 + (-t40 * t150 - t42 * t151) * t102 + (t23 + (-qJD(5) * t41 - t111 * t42) * t106 + (-qJD(5) * t43 + t111 * t40) * t105) * t101) * t102);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
