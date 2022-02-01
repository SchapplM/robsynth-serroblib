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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:31
% DurationCPUTime: 1.45s
% Computational Cost: add. (4786->192), mult. (2744->275), div. (0->0), fcn. (1920->10), ass. (0->109)
t102 = sin(qJ(5));
t135 = qJD(5) * t102;
t104 = cos(qJ(5));
t134 = qJD(5) * t104;
t101 = qJ(1) + pkin(9);
t160 = -pkin(2) * cos(t101) - cos(qJ(1)) * pkin(1);
t136 = Icges(6,4) * t104;
t119 = -Icges(6,2) * t102 + t136;
t98 = qJ(3) + t101;
t94 = qJ(4) + t98;
t90 = cos(t94);
t112 = t119 * t90;
t89 = sin(t94);
t43 = Icges(6,6) * t89 + t112;
t137 = Icges(6,4) * t102;
t120 = Icges(6,1) * t104 - t137;
t113 = t120 * t90;
t45 = Icges(6,5) * t89 + t113;
t122 = t102 * t43 - t104 * t45;
t159 = t122 * t89;
t42 = -Icges(6,6) * t90 + t119 * t89;
t44 = -Icges(6,5) * t90 + t120 * t89;
t123 = t102 * t42 - t104 * t44;
t158 = t123 * t90;
t118 = Icges(6,5) * t104 - Icges(6,6) * t102;
t78 = Icges(6,2) * t104 + t137;
t79 = Icges(6,1) * t102 + t136;
t121 = t102 * t78 - t104 * t79;
t100 = qJD(1) + qJD(3);
t97 = qJD(4) + t100;
t157 = t118 * qJD(5) + t121 * t97;
t156 = 2 * m(4);
t155 = 2 * m(5);
t154 = 2 * m(6);
t92 = sin(t98);
t151 = pkin(3) * t92;
t81 = t89 * rSges(6,3);
t147 = t89 * t97;
t146 = t90 * t97;
t141 = rSges(6,2) * t102;
t74 = t89 * t141;
t145 = rSges(6,3) * t146 + t97 * t74;
t144 = t90 * rSges(6,3) + t74;
t142 = rSges(6,1) * t104;
t140 = t100 * t92;
t93 = cos(t98);
t139 = t100 * t93;
t138 = t104 * t89;
t133 = pkin(3) * t140;
t132 = pkin(3) * t139;
t127 = rSges(6,2) * t134;
t130 = t90 * t141;
t131 = -t97 * t130 + (-t135 * rSges(6,1) - t127) * t89;
t126 = -pkin(4) - t142;
t61 = t93 * rSges(4,1) - rSges(4,2) * t92;
t57 = t90 * rSges(5,1) - rSges(5,2) * t89;
t49 = -rSges(5,1) * t146 + rSges(5,2) * t147;
t88 = pkin(3) * t93;
t53 = t57 + t88;
t55 = -rSges(4,1) * t139 + rSges(4,2) * t140;
t125 = -pkin(2) * sin(t101) - sin(qJ(1)) * pkin(1);
t60 = -rSges(4,1) * t92 - rSges(4,2) * t93;
t56 = -rSges(5,1) * t89 - rSges(5,2) * t90;
t80 = rSges(6,1) * t102 + rSges(6,2) * t104;
t47 = t90 * t142 - t130 + t81;
t77 = Icges(6,5) * t102 + Icges(6,6) * t104;
t48 = t56 * t97;
t54 = t60 * t100;
t117 = t125 * qJD(1);
t116 = t160 * qJD(1);
t115 = (t120 - t78) * t135 + (t119 + t79) * t134;
t114 = (-qJD(5) * t122 + t157 * t89 + (-t102 * t120 - t104 * t119) * t147) * t89 / 0.2e1 - (-qJD(5) * t123 - t157 * t90 + (t102 * t113 + t104 * t112) * t97) * t90 / 0.2e1 + (t102 * t44 + t104 * t42 - t121 * t89 - t77 * t90) * t147 / 0.2e1 + (t102 * t45 + t104 * t43 - t121 * t90 + t77 * t89) * t146 / 0.2e1;
t111 = t118 * t90;
t52 = t56 - t151;
t33 = t90 * pkin(4) + t89 * pkin(8) + t47;
t35 = t49 - t132;
t29 = t33 + t88;
t32 = t90 * pkin(8) + t126 * t89 + t144;
t34 = t48 - t133;
t107 = Icges(6,3) * t97 - qJD(5) * t77;
t28 = t32 - t151;
t15 = (t126 * t90 + (-rSges(6,3) - pkin(8)) * t89) * t97 - t131;
t13 = t15 - t132;
t14 = -t90 * t127 - pkin(4) * t147 + pkin(8) * t146 + (-t135 * t90 - t138 * t97) * rSges(6,1) + t145;
t12 = t14 - t133;
t73 = (-t141 + t142) * qJD(5);
t51 = t61 - t160;
t50 = t125 + t60;
t46 = rSges(6,1) * t138 - t144;
t41 = Icges(6,3) * t89 + t111;
t40 = -Icges(6,3) * t90 + t118 * t89;
t39 = t116 + t55;
t38 = t54 + t117;
t37 = t53 - t160;
t36 = t125 + t52;
t31 = t116 + t35;
t30 = t117 + t34;
t27 = t29 - t160;
t26 = t125 + t28;
t19 = t107 * t89 + t111 * t97;
t18 = t107 * t90 - t118 * t147;
t11 = t116 + t13;
t10 = t117 + t12;
t9 = -t122 * t90 + t41 * t89;
t8 = t40 * t89 - t158;
t7 = -t41 * t90 - t159;
t6 = -t123 * t89 - t40 * t90;
t1 = ((-t47 + t81) * t97 + t131) * t89 + (-qJD(5) * t80 * t90 + t97 * t46 + t145) * t90;
t2 = [(t10 * t27 + t11 * t26) * t154 + (t30 * t37 + t31 * t36) * t155 + (t38 * t51 + t39 * t50) * t156 + t115; 0; 0; m(6) * (t10 * t29 + t11 * t28 + t12 * t27 + t13 * t26) + m(5) * (t30 * t53 + t31 * t52 + t34 * t37 + t35 * t36) + m(4) * (t38 * t61 + t39 * t60 + t50 * t55 + t51 * t54) + t115; 0; (t54 * t61 + t55 * t60) * t156 + (t34 * t53 + t35 * t52) * t155 + (t12 * t29 + t13 * t28) * t154 + t115; m(6) * (t10 * t33 + t11 * t32 + t14 * t27 + t15 * t26) + m(5) * (t30 * t57 + t31 * t56 + t36 * t49 + t37 * t48) + t115; 0; m(5) * (t34 * t57 + t35 * t56 + t48 * t53 + t49 * t52) + m(6) * (t12 * t33 + t13 * t32 + t14 * t29 + t15 * t28) + t115; (t48 * t57 + t49 * t56) * t155 + (t14 * t33 + t15 * t32) * t154 + t115; m(6) * ((-t26 * t90 - t27 * t89) * t73 + ((-t27 * t97 - t11) * t90 + (t26 * t97 - t10) * t89) * t80) + t114; m(6) * t1; m(6) * ((-t28 * t90 - t29 * t89) * t73 + ((-t29 * t97 - t13) * t90 + (t28 * t97 - t12) * t89) * t80) + t114; m(6) * ((-t32 * t90 - t33 * t89) * t73 + ((-t33 * t97 - t15) * t90 + (t32 * t97 - t14) * t89) * t80) + t114; ((t46 * t89 + t47 * t90) * t1 + (t89 ^ 2 + t90 ^ 2) * t80 * t73) * t154 + (-t8 * t90 + t89 * t9) * t146 + t89 * ((t89 * t18 + (t8 + t159) * t97) * t89 + (t9 * t97 + (t134 * t42 + t135 * t44) * t90 + (-t19 + (-qJD(5) * t43 + t44 * t97) * t104 + (-qJD(5) * t45 - t42 * t97) * t102) * t89) * t90) + (-t6 * t90 + t7 * t89) * t147 - t90 * ((t90 * t19 + (t7 + t158) * t97) * t90 + (t6 * t97 + (-t134 * t43 - t135 * t45) * t89 + (-t18 + (qJD(5) * t42 + t45 * t97) * t104 + (qJD(5) * t44 - t43 * t97) * t102) * t90) * t89);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
