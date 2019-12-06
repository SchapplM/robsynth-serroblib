% Calculate time derivative of joint inertia matrix for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:59
% EndTime: 2019-12-05 17:31:05
% DurationCPUTime: 1.56s
% Computational Cost: add. (4249->307), mult. (12011->451), div. (0->0), fcn. (13377->10), ass. (0->141)
t120 = cos(pkin(7));
t122 = sin(qJ(1));
t159 = cos(pkin(8));
t143 = t122 * t159;
t117 = sin(pkin(8));
t124 = cos(qJ(1));
t153 = t124 * t117;
t106 = t120 * t153 - t143;
t167 = 2 * m(6);
t118 = sin(pkin(7));
t166 = 0.2e1 * t118;
t165 = m(5) / 0.2e1;
t164 = m(6) / 0.2e1;
t163 = -rSges(6,3) - pkin(6);
t162 = pkin(2) * t120;
t161 = -rSges(3,3) - qJ(2);
t160 = -rSges(5,3) - qJ(4);
t158 = t117 * t118;
t157 = t118 * t122;
t156 = t118 * t124;
t155 = t122 * qJ(2);
t154 = t122 * t117;
t115 = t124 * qJ(2);
t152 = qJD(1) * t122;
t151 = qJD(1) * t124;
t150 = qJD(3) * t118;
t149 = t165 + t164;
t121 = sin(qJ(5));
t123 = cos(qJ(5));
t142 = t124 * t159;
t104 = t120 * t154 + t142;
t105 = -t120 * t143 + t153;
t116 = sin(pkin(9));
t119 = cos(pkin(9));
t88 = t105 * t119 - t116 * t157;
t135 = t104 * t121 - t88 * t123;
t145 = t118 * t151;
t107 = t120 * t142 + t154;
t98 = t107 * qJD(1);
t81 = -t116 * t145 - t98 * t119;
t97 = t106 * qJD(1);
t49 = t135 * qJD(5) - t81 * t121 - t97 * t123;
t71 = -t104 * t123 - t88 * t121;
t50 = t71 * qJD(5) - t97 * t121 + t81 * t123;
t80 = -t98 * t116 + t119 * t145;
t28 = t50 * rSges(6,1) + t49 * rSges(6,2) + t80 * rSges(6,3);
t87 = t105 * t116 + t119 * t157;
t39 = -rSges(6,1) * t135 + t71 * rSges(6,2) + t87 * rSges(6,3);
t147 = -pkin(1) - t162;
t146 = t118 * t152;
t144 = t118 * t159;
t141 = pkin(1) * t152 - qJD(2) * t122;
t102 = t116 * t144 + t120 * t119;
t103 = -t120 * t116 + t119 * t144;
t85 = -t103 * t121 + t123 * t158;
t86 = t103 * t123 + t121 * t158;
t59 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t102;
t60 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t102;
t82 = t85 * qJD(5);
t83 = t86 * qJD(5);
t63 = Icges(6,5) * t82 - Icges(6,6) * t83;
t64 = Icges(6,4) * t82 - Icges(6,2) * t83;
t65 = Icges(6,1) * t82 - Icges(6,4) * t83;
t140 = t102 * t63 - t83 * t59 + t82 * t60 + t85 * t64 + t86 * t65;
t90 = t107 * t119 + t116 * t156;
t74 = t106 * t121 + t90 * t123;
t96 = t105 * qJD(1);
t79 = -t116 * t146 + t96 * t119;
t95 = t104 * qJD(1);
t47 = -t74 * qJD(5) - t79 * t121 - t95 * t123;
t73 = t106 * t123 - t90 * t121;
t48 = t73 * qJD(5) - t95 * t121 + t79 * t123;
t139 = -t48 * rSges(6,1) - t47 * rSges(6,2);
t138 = -t74 * rSges(6,1) - t73 * rSges(6,2);
t114 = qJD(2) * t124;
t137 = -t122 * t150 + t114;
t136 = rSges(3,1) * t120 - rSges(3,2) * t118;
t134 = -pkin(1) - t136;
t133 = -qJ(3) * t118 + t147;
t132 = qJ(3) * t146 + t152 * t162 + t141;
t131 = (-rSges(4,3) - qJ(3)) * t118 + t147;
t130 = t105 * pkin(3) + t133 * t122 + t115;
t129 = t133 * t124 - t155;
t128 = t131 * t124 - t155;
t92 = t161 * t122 + t134 * t124;
t127 = -t107 * pkin(3) - t106 * qJ(4) + t129;
t126 = t95 * qJ(4) + (-qJ(2) * qJD(1) - t150) * t124 - t106 * qJD(4) - t96 * pkin(3) + t132;
t125 = -t98 * pkin(3) + t129 * qJD(1) - t104 * qJD(4) + t137;
t91 = t124 * rSges(3,3) + t134 * t122 + t115;
t89 = t107 * t116 - t119 * t156;
t78 = t96 * t116 + t119 * t146;
t76 = t92 * qJD(1) + t114;
t75 = (t136 * t122 + t161 * t124) * qJD(1) + t141;
t68 = -t107 * rSges(4,1) + t106 * rSges(4,2) + t128;
t67 = t105 * rSges(4,1) + t104 * rSges(4,2) + t131 * t122 + t115;
t66 = t82 * rSges(6,1) - t83 * rSges(6,2);
t61 = t86 * rSges(6,1) + t85 * rSges(6,2) + t102 * rSges(6,3);
t58 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t102;
t55 = -t98 * rSges(4,1) + t97 * rSges(4,2) + t128 * qJD(1) + t137;
t54 = -t124 * t150 - t96 * rSges(4,1) - t95 * rSges(4,2) + (rSges(4,3) * t157 - t115) * qJD(1) + t132;
t42 = -t90 * rSges(5,1) + t89 * rSges(5,2) - t106 * rSges(5,3) + t127;
t41 = t88 * rSges(5,1) - t87 * rSges(5,2) + t160 * t104 + t130;
t40 = t89 * rSges(6,3) - t138;
t38 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t89;
t37 = -Icges(6,1) * t135 + Icges(6,4) * t71 + Icges(6,5) * t87;
t36 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t89;
t35 = -Icges(6,4) * t135 + Icges(6,2) * t71 + Icges(6,6) * t87;
t34 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t89;
t33 = -Icges(6,5) * t135 + Icges(6,6) * t71 + Icges(6,3) * t87;
t32 = t81 * rSges(5,1) - t80 * rSges(5,2) + t160 * t97 + t125;
t31 = -t79 * rSges(5,1) + t78 * rSges(5,2) + t95 * rSges(5,3) + t126;
t30 = -t90 * pkin(4) + t163 * t89 + t127 + t138;
t29 = t88 * pkin(4) + t87 * pkin(6) - t104 * qJ(4) + t130 + t39;
t27 = t78 * rSges(6,3) - t139;
t26 = Icges(6,1) * t50 + Icges(6,4) * t49 + Icges(6,5) * t80;
t25 = Icges(6,1) * t48 + Icges(6,4) * t47 + Icges(6,5) * t78;
t24 = Icges(6,4) * t50 + Icges(6,2) * t49 + Icges(6,6) * t80;
t23 = Icges(6,4) * t48 + Icges(6,2) * t47 + Icges(6,6) * t78;
t22 = Icges(6,5) * t50 + Icges(6,6) * t49 + Icges(6,3) * t80;
t21 = Icges(6,5) * t48 + Icges(6,6) * t47 + Icges(6,3) * t78;
t20 = -t102 * t40 + t89 * t61;
t19 = t102 * t39 - t87 * t61;
t18 = t89 * t58 + t73 * t59 + t74 * t60;
t17 = -t135 * t60 + t87 * t58 + t71 * t59;
t16 = t81 * pkin(4) + t80 * pkin(6) - t97 * qJ(4) + t125 + t28;
t15 = -t79 * pkin(4) + t163 * t78 + t126 + t139;
t14 = t102 * t28 - t80 * t61 - t87 * t66;
t13 = -t102 * t27 + t78 * t61 + t89 * t66;
t12 = t102 * t34 + t85 * t36 + t86 * t38;
t11 = t102 * t33 + t85 * t35 + t86 * t37;
t10 = t140 * t102;
t9 = t89 * t34 + t73 * t36 + t74 * t38;
t8 = t89 * t33 + t73 * t35 + t74 * t37;
t7 = -t135 * t38 + t87 * t34 + t71 * t36;
t6 = -t135 * t37 + t87 * t33 + t71 * t35;
t5 = -t135 * t65 + t49 * t59 + t50 * t60 + t80 * t58 + t87 * t63 + t71 * t64;
t4 = t47 * t59 + t48 * t60 + t78 * t58 + t89 * t63 + t73 * t64 + t74 * t65;
t3 = t87 * t27 - t89 * t28 - t78 * t39 + t80 * t40;
t2 = t102 * t21 + t85 * t23 + t86 * t25 - t83 * t36 + t82 * t38;
t1 = t102 * t22 + t85 * t24 + t86 * t26 - t83 * t35 + t82 * t37;
t43 = [0.2e1 * m(3) * (t92 * t75 + t91 * t76) + 0.2e1 * m(4) * (t68 * t54 + t67 * t55) + 0.2e1 * m(5) * (t42 * t31 + t41 * t32) + (t30 * t15 + t29 * t16) * t167 + t140; m(3) * (t122 * t76 + t124 * t75 + (-t122 * t92 + t124 * t91) * qJD(1)) + m(4) * (t122 * t55 + t124 * t54 + (-t122 * t68 + t124 * t67) * qJD(1)) + m(5) * (t122 * t32 + t124 * t31 + (-t122 * t42 + t124 * t41) * qJD(1)) + m(6) * (t122 * t16 + t124 * t15 + (-t122 * t30 + t124 * t29) * qJD(1)); 0; (m(4) * (-t122 * t54 + t124 * t55 - t68 * t151 - t67 * t152) / 0.2e1 + (-t122 * t31 + t124 * t32 - t42 * t151 - t41 * t152) * t165 + (-t122 * t15 + t124 * t16 - t30 * t151 - t29 * t152) * t164) * t166; 0; 0; m(5) * (-t104 * t31 + t106 * t32 - t95 * t41 - t97 * t42) + m(6) * (-t104 * t15 + t106 * t16 - t95 * t29 - t97 * t30); 0.2e1 * t149 * (-t95 * t122 - t97 * t124 + (t104 * t122 + t106 * t124) * qJD(1)); t149 * (t122 * t97 - t124 * t95 + (t104 * t124 - t106 * t122) * qJD(1)) * t166; 0.4e1 * t149 * (t104 * t97 - t106 * t95); t10 + m(6) * (t13 * t30 + t14 * t29 + t20 * t15 + t19 * t16) + (t2 / 0.2e1 + t4 / 0.2e1) * t89 + (t1 / 0.2e1 + t5 / 0.2e1) * t87 + (t11 / 0.2e1 + t17 / 0.2e1) * t80 + (t12 / 0.2e1 + t18 / 0.2e1) * t78; m(6) * (t14 * t122 + t13 * t124 + (-t122 * t20 + t124 * t19) * qJD(1)); m(6) * (-t3 * t120 + (-t122 * t13 + t124 * t14 + (-t122 * t19 - t124 * t20) * qJD(1)) * t118); m(6) * (-t13 * t104 + t14 * t106 + t3 * t158 - t19 * t95 - t20 * t97); ((-t89 * t39 + t87 * t40) * t3 + t20 * t13 + t19 * t14) * t167 + t102 * (t1 * t87 + t11 * t80 + t12 * t78 + t2 * t89 + t10) + t80 * (t17 * t102 + t6 * t87 + t7 * t89) + t87 * (t5 * t102 + (-t135 * t26 + t87 * t22 + t71 * t24 + t80 * t33 + t49 * t35 + t50 * t37) * t87 + t6 * t80 + (-t135 * t25 + t87 * t21 + t71 * t23 + t80 * t34 + t49 * t36 + t50 * t38) * t89 + t7 * t78) + t78 * (t18 * t102 + t8 * t87 + t9 * t89) + t89 * (t4 * t102 + (t89 * t22 + t73 * t24 + t74 * t26 + t78 * t33 + t47 * t35 + t48 * t37) * t87 + t8 * t80 + (t89 * t21 + t73 * t23 + t74 * t25 + t78 * t34 + t47 * t36 + t48 * t38) * t89 + t9 * t78);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t43(1), t43(2), t43(4), t43(7), t43(11); t43(2), t43(3), t43(5), t43(8), t43(12); t43(4), t43(5), t43(6), t43(9), t43(13); t43(7), t43(8), t43(9), t43(10), t43(14); t43(11), t43(12), t43(13), t43(14), t43(15);];
Mq = res;
