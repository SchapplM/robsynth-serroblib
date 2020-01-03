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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:48
% DurationCPUTime: 1.86s
% Computational Cost: add. (3820->201), mult. (2408->276), div. (0->0), fcn. (1712->10), ass. (0->121)
t118 = pkin(9) + qJ(5);
t111 = sin(t118);
t153 = qJD(5) * t111;
t112 = cos(t118);
t152 = qJD(5) * t112;
t119 = qJD(1) + qJD(2);
t136 = Icges(6,5) * t112 - Icges(6,6) * t111;
t160 = Icges(6,4) * t112;
t137 = -Icges(6,2) * t111 + t160;
t161 = Icges(6,4) * t111;
t138 = Icges(6,1) * t112 - t161;
t69 = Icges(6,2) * t112 + t161;
t70 = Icges(6,1) * t111 + t160;
t139 = t111 * t69 - t112 * t70;
t188 = -t136 * qJD(5) + (t111 * t138 + t112 * t137 - t139) * t119;
t166 = rSges(5,2) * sin(pkin(9));
t185 = pkin(3) - t166;
t184 = rSges(5,3) + qJ(4);
t165 = rSges(6,2) * t111;
t168 = rSges(6,1) * t112;
t183 = -t165 + t168;
t120 = qJ(1) + qJ(2);
t113 = pkin(8) + t120;
t106 = cos(t113);
t114 = sin(t120);
t108 = pkin(2) * t114;
t123 = -pkin(7) - qJ(4);
t182 = -t106 * t123 - t108;
t105 = sin(t113);
t41 = -Icges(6,6) * t105 - t106 * t137;
t43 = -Icges(6,5) * t105 - t106 * t138;
t141 = t111 * t41 - t112 * t43;
t181 = t105 * t141;
t115 = cos(t120);
t74 = t114 * rSges(3,1) + t115 * rSges(3,2);
t54 = t105 * rSges(4,1) + t106 * rSges(4,2) + t108;
t68 = Icges(6,5) * t111 + Icges(6,6) * t112;
t180 = -Icges(6,3) * t119 + qJD(5) * t68;
t176 = 2 * m(3);
t175 = 2 * m(4);
t174 = 2 * m(5);
t173 = 2 * m(6);
t156 = t105 * t119;
t79 = t106 * t168;
t170 = rSges(6,3) * t156 + t119 * t79;
t122 = cos(pkin(9));
t169 = rSges(5,1) * t122;
t167 = rSges(3,2) * t114;
t164 = pkin(1) * qJD(1);
t155 = t106 * t119;
t154 = t115 * t119;
t151 = t119 * t166;
t150 = t119 * t165;
t124 = sin(qJ(1));
t149 = t124 * t164;
t89 = t106 * t169;
t107 = pkin(4) * t122 + pkin(3);
t73 = rSges(6,1) * t111 + rSges(6,2) * t112;
t132 = t73 * qJD(5);
t126 = rSges(6,3) * t155 + t105 * t150 - t106 * t132;
t93 = qJD(4) * t105;
t12 = t93 + ((-t107 - t168) * t105 + t182) * t119 + t126;
t10 = t12 - t149;
t116 = t124 * pkin(1);
t44 = -rSges(6,3) * t106 + t105 * t183;
t30 = t105 * t107 - t182 + t44;
t28 = t116 + t30;
t146 = -t119 * t28 - t10;
t125 = cos(qJ(1));
t110 = t125 * t164;
t91 = pkin(2) * t154;
t13 = t107 * t155 + t91 + (-qJD(4) - t150) * t106 + (-t119 * t123 - t132) * t105 + t170;
t11 = t110 + t13;
t117 = t125 * pkin(1);
t109 = pkin(2) * t115;
t45 = -t105 * rSges(6,3) + t106 * t165 - t79;
t31 = -t105 * t123 + t106 * t107 + t109 - t45;
t29 = t117 + t31;
t145 = t119 * t29 - t11;
t144 = -t119 * t30 - t12;
t143 = t119 * t31 - t13;
t75 = t115 * rSges(3,1) - t167;
t58 = rSges(3,1) * t154 - t119 * t167;
t55 = t106 * rSges(4,1) - rSges(4,2) * t105 + t109;
t40 = -Icges(6,6) * t106 + t105 * t137;
t42 = -Icges(6,5) * t106 + t105 * t138;
t142 = t111 * t40 - t112 * t42;
t47 = rSges(4,1) * t155 - rSges(4,2) * t156 + t91;
t135 = (t138 - t69) * t153 + (t137 + t70) * t152;
t134 = -(-qJD(5) * t141 + t188 * t105) * t105 / 0.2e1 - (-qJD(5) * t142 + t188 * t106) * t106 / 0.2e1 + (-t105 * t139 - t106 * t68 + t111 * t42 + t112 * t40) * t156 / 0.2e1 - (-t105 * t68 + t106 * t139 + t111 * t43 + t112 * t41) * t155 / 0.2e1;
t57 = t74 * t119;
t133 = t142 * t106;
t129 = t136 * t119;
t35 = t184 * t105 + t106 * t185 + t109 + t89;
t46 = t54 * t119;
t34 = t108 - t184 * t106 + (t169 + t185) * t105;
t21 = t119 * t89 + pkin(3) * t155 + t91 + (-qJD(4) - t151) * t106 + t184 * t156;
t20 = t105 * t151 + t93 + (-t108 + (-pkin(3) - t169) * t105) * t119 + t184 * t155;
t64 = t183 * qJD(5);
t60 = t117 + t75;
t59 = t116 + t74;
t51 = t110 + t58;
t50 = -t57 - t149;
t49 = t117 + t55;
t48 = t116 + t54;
t39 = -Icges(6,3) * t105 - t106 * t136;
t38 = -Icges(6,3) * t106 + t105 * t136;
t37 = t110 + t47;
t36 = -t46 - t149;
t33 = t117 + t35;
t32 = t116 + t34;
t23 = -t180 * t105 + t106 * t129;
t22 = t105 * t129 + t180 * t106;
t17 = t110 + t21;
t16 = t20 - t149;
t9 = -t105 * t39 + t106 * t141;
t8 = -t105 * t38 + t133;
t7 = -t106 * t39 - t181;
t6 = -t105 * t142 - t106 * t38;
t3 = (t119 * t44 + t126) * t106 + (-t105 * t132 + (t45 + (-t165 - t168) * t106) * t119 + t170) * t105;
t1 = [(t50 * t60 + t51 * t59) * t176 + (t36 * t49 + t37 * t48) * t175 + (t16 * t33 + t17 * t32) * t174 + (t10 * t29 + t11 * t28) * t173 + t135; m(3) * (t50 * t75 + t51 * t74 - t57 * t60 + t58 * t59) + m(4) * (t36 * t55 + t37 * t54 - t46 * t49 + t47 * t48) + m(5) * (t16 * t35 + t17 * t34 + t20 * t33 + t21 * t32) + m(6) * (t10 * t31 + t11 * t30 + t12 * t29 + t13 * t28) + t135; (t12 * t31 + t13 * t30) * t173 + (-t46 * t55 + t47 * t54) * t175 + (t20 * t35 + t21 * t34) * t174 + (-t57 * t75 + t58 * t74) * t176 + t135; 0; 0; 0; m(5) * ((-t119 * t32 - t16) * t106 + (t119 * t33 - t17) * t105) + m(6) * (t105 * t145 + t106 * t146); m(6) * (t105 * t143 + t106 * t144) + m(5) * ((-t119 * t34 - t20) * t106 + (t119 * t35 - t21) * t105); 0; 0; m(6) * ((-t105 * t29 + t106 * t28) * t64 + (t105 * t146 - t106 * t145) * t73) + t134; m(6) * ((-t105 * t31 + t106 * t30) * t64 + (t105 * t144 - t106 * t143) * t73) + t134; m(6) * t3; 0; ((t105 * t44 - t106 * t45) * t3 + (t105 ^ 2 + t106 ^ 2) * t73 * t64) * t173 + (-t105 * t7 - t106 * t6) * t156 - t106 * ((t106 * t23 + (-t7 + t133) * t119) * t106 + (t6 * t119 + (t41 * t152 + t43 * t153) * t105 + (t22 + (qJD(5) * t40 - t119 * t43) * t112 + (qJD(5) * t42 + t119 * t41) * t111) * t106) * t105) - (-t105 * t9 - t106 * t8) * t155 - t105 * ((t105 * t22 + (t8 + t181) * t119) * t105 + (-t9 * t119 + (-t40 * t152 - t42 * t153) * t106 + (t23 + (-qJD(5) * t41 - t119 * t42) * t112 + (-qJD(5) * t43 + t119 * t40) * t111) * t105) * t106);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
