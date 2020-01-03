% Calculate time derivative of joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:07
% DurationCPUTime: 1.61s
% Computational Cost: add. (4786->196), mult. (2744->282), div. (0->0), fcn. (1920->10), ass. (0->109)
t111 = qJ(1) + pkin(9);
t149 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t111);
t112 = sin(qJ(5));
t142 = qJD(5) * t112;
t114 = cos(qJ(5));
t141 = qJD(5) * t114;
t107 = qJ(3) + t111;
t102 = qJ(4) + t107;
t96 = sin(t102);
t97 = cos(t102);
t56 = t96 * rSges(5,1) + t97 * rSges(5,2);
t148 = Icges(6,4) * t112;
t125 = Icges(6,1) * t114 - t148;
t170 = Icges(6,5) * t96 + t125 * t97;
t147 = Icges(6,4) * t114;
t124 = -Icges(6,2) * t112 + t147;
t171 = Icges(6,6) * t96 + t124 * t97;
t128 = -t112 * t171 + t114 * t170;
t174 = t128 * t96;
t150 = pkin(2) * sin(t111) + sin(qJ(1)) * pkin(1);
t123 = Icges(6,5) * t114 - Icges(6,6) * t112;
t172 = Icges(6,3) * t96 + t123 * t97;
t168 = 2 * m(4);
t167 = 2 * m(5);
t166 = 2 * m(6);
t110 = qJD(1) + qJD(3);
t106 = qJD(4) + t110;
t160 = rSges(6,2) * t112;
t140 = t96 * t160;
t154 = t106 * t97;
t163 = rSges(6,3) * t154 + t106 * t140;
t155 = t106 * t96;
t161 = rSges(6,1) * t114;
t76 = t97 * t161;
t162 = rSges(6,3) * t155 + t106 * t76;
t100 = sin(t107);
t101 = cos(t107);
t60 = t100 * rSges(4,1) + t101 * rSges(4,2);
t40 = -Icges(6,3) * t97 + t123 * t96;
t156 = t106 * t40;
t153 = t112 * t97;
t152 = t114 * t96;
t151 = t149 * qJD(1);
t146 = qJD(5) * t96;
t145 = qJD(5) * t97;
t144 = t100 * t110;
t143 = t101 * t110;
t94 = pkin(3) * t100;
t52 = t94 + t56;
t139 = pkin(3) * t144;
t57 = t97 * rSges(5,1) - rSges(5,2) * t96;
t61 = t101 * rSges(4,1) - rSges(4,2) * t100;
t95 = pkin(3) * t101;
t53 = t57 + t95;
t49 = rSges(5,1) * t154 - rSges(5,2) * t155;
t136 = rSges(6,1) * t152 - t140;
t55 = rSges(4,1) * t143 - rSges(4,2) * t144;
t82 = rSges(6,1) * t112 + rSges(6,2) * t114;
t42 = -Icges(6,6) * t97 + t124 * t96;
t44 = -Icges(6,5) * t97 + t125 * t96;
t131 = t112 * t44 + t114 * t42;
t130 = t112 * t42 - t114 * t44;
t129 = -t112 * t170 - t114 * t171;
t80 = Icges(6,2) * t114 + t148;
t81 = Icges(6,1) * t112 + t147;
t126 = t112 * t80 - t114 * t81;
t78 = pkin(3) * t143;
t35 = t49 + t78;
t47 = rSges(6,2) * t153 - t96 * rSges(6,3) - t76;
t79 = Icges(6,5) * t112 + Icges(6,6) * t114;
t48 = t56 * t106;
t122 = (t125 - t80) * t142 + (t124 + t81) * t141;
t121 = t150 * qJD(1);
t120 = t130 * t97;
t116 = -t123 * qJD(5) - t106 * t126;
t119 = -(-qJD(5) * t128 + t112 * (t106 * t44 + t81 * t145) + t114 * (t106 * t42 + t80 * t145) + t116 * t96) * t96 / 0.2e1 - (-qJD(5) * t130 + t112 * (t170 * t106 - t81 * t146) + t114 * (t171 * t106 - t80 * t146) + t116 * t97) * t97 / 0.2e1 + (-t126 * t96 - t79 * t97 + t131) * t155 / 0.2e1 - (t126 * t97 - t79 * t96 + t129) * t154 / 0.2e1;
t54 = t60 * t110;
t118 = qJD(5) * t82;
t33 = t97 * pkin(4) + t96 * pkin(8) - t47;
t29 = t33 + t95;
t32 = t96 * pkin(4) + (-rSges(6,3) - pkin(8)) * t97 + t136;
t28 = t94 + t32;
t34 = -t48 - t139;
t15 = -t96 * rSges(6,1) * t142 + pkin(8) * t155 + pkin(4) * t154 + (-t106 * t153 - t96 * t141) * rSges(6,2) + t162;
t13 = t78 + t15;
t14 = -rSges(6,2) * t97 * t141 - pkin(4) * t155 + pkin(8) * t154 + (-t106 * t152 - t97 * t142) * rSges(6,1) + t163;
t12 = t14 - t139;
t74 = (-t160 + t161) * qJD(5);
t51 = t61 + t149;
t50 = t150 + t60;
t46 = -rSges(6,3) * t97 + t136;
t39 = t55 + t151;
t38 = -t54 - t121;
t37 = t53 + t149;
t36 = t52 + t150;
t31 = t35 + t151;
t30 = -t121 + t34;
t27 = t29 + t149;
t26 = t28 + t150;
t19 = t106 * t172 - t79 * t146;
t18 = t79 * t145 + t156;
t11 = t13 + t151;
t10 = -t121 + t12;
t9 = t128 * t97 + t172 * t96;
t8 = -t40 * t96 + t120;
t7 = t172 * t97 - t174;
t6 = -t130 * t96 - t40 * t97;
t1 = (t106 * t46 - t118 * t97 + t163) * t97 + (-t96 * t118 + (t47 + (-t160 - t161) * t97) * t106 + t162) * t96;
t2 = [(t38 * t51 + t39 * t50) * t168 + (t30 * t37 + t31 * t36) * t167 + (t10 * t27 + t11 * t26) * t166 + t122; 0; 0; m(4) * (t38 * t61 + t39 * t60 + t50 * t55 - t51 * t54) + m(5) * (t30 * t53 + t31 * t52 + t34 * t37 + t35 * t36) + m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + t122; 0; (-t54 * t61 + t55 * t60) * t168 + (t34 * t53 + t35 * t52) * t167 + (t12 * t29 + t13 * t28) * t166 + t122; m(5) * (t30 * t57 + t31 * t56 + t36 * t49 - t37 * t48) + m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + t122; 0; m(5) * (t34 * t57 + t35 * t56 - t48 * t53 + t49 * t52) + m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + t122; (-t48 * t57 + t49 * t56) * t167 + (t14 * t33 + t15 * t32) * t166 + t122; m(6) * ((t26 * t97 - t27 * t96) * t74 + (-t10 * t96 + t11 * t97 + (-t26 * t96 - t27 * t97) * t106) * t82) + t119; m(6) * t1; m(6) * ((t28 * t97 - t29 * t96) * t74 + (-t12 * t96 + t13 * t97 + (-t28 * t96 - t29 * t97) * t106) * t82) + t119; m(6) * ((t32 * t97 - t33 * t96) * t74 + (-t14 * t96 + t15 * t97 + (-t32 * t96 - t33 * t97) * t106) * t82) + t119; ((t46 * t96 - t47 * t97) * t1 + (t96 ^ 2 + t97 ^ 2) * t82 * t74) * t166 + (-t6 * t97 - t7 * t96) * t155 - t97 * ((t97 * t19 + (-t7 + t120) * t106) * t97 + (t6 * t106 + (-t141 * t171 - t142 * t170) * t96 + (t18 + t131 * qJD(5) + (t128 - t40) * t106) * t97) * t96) - (-t8 * t97 - t9 * t96) * t154 - t96 * ((t96 * t18 + (t8 + t174) * t106) * t96 + (-t9 * t106 + (-t42 * t141 - t44 * t142 + t156) * t97 + (-t129 * qJD(5) + t130 * t106 + t19) * t96) * t97);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
