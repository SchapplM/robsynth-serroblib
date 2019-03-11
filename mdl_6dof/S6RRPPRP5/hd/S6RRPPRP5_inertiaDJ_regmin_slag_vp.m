% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:07
% EndTime: 2019-03-09 08:44:10
% DurationCPUTime: 1.26s
% Computational Cost: add. (1818->191), mult. (3935->338), div. (0->0), fcn. (3320->6), ass. (0->97)
t127 = pkin(3) + pkin(7);
t125 = cos(qJ(5));
t82 = sin(qJ(2));
t116 = t82 * qJ(3);
t80 = -pkin(2) - qJ(4);
t83 = cos(qJ(2));
t94 = -t80 * t83 + t116;
t47 = -pkin(1) - t94;
t59 = t127 * t82;
t79 = cos(pkin(9));
t51 = t79 * t59;
t78 = sin(pkin(9));
t22 = t82 * pkin(4) + t51 + (pkin(8) * t83 - t47) * t78;
t121 = t79 * t83;
t31 = t79 * t47 + t78 * t59;
t26 = -pkin(8) * t121 + t31;
t81 = sin(qJ(5));
t131 = t125 * t26 + t81 * t22;
t100 = qJD(5) * t125;
t115 = qJD(5) * t81;
t44 = t79 * t100 - t78 * t115;
t103 = t125 * t79;
t119 = t81 * t78;
t48 = -t103 + t119;
t55 = (t78 ^ 2 + t79 ^ 2) * qJD(4);
t112 = t83 * qJD(3);
t130 = t94 * qJD(2) + qJD(4) * t82 - t112;
t122 = t78 * t82;
t114 = t82 * qJD(2);
t99 = pkin(2) * t114 - t82 * qJD(3);
t35 = -t83 * qJD(4) + (-qJ(3) * t83 + qJ(4) * t82) * qJD(2) + t99;
t72 = t83 * qJD(2);
t70 = pkin(7) * t72;
t53 = pkin(3) * t72 + t70;
t19 = -t78 * t35 + t79 * t53;
t11 = (pkin(4) * t83 - pkin(8) * t122) * qJD(2) + t19;
t105 = t79 * t114;
t20 = t79 * t35 + t78 * t53;
t14 = pkin(8) * t105 + t20;
t4 = -qJD(5) * t131 + t125 * t11 - t81 * t14;
t129 = 0.2e1 * qJD(3);
t128 = 2 * qJD(6);
t126 = -pkin(8) + t80;
t43 = -t78 * t100 - t79 * t115;
t124 = t48 * t43;
t52 = t127 * t114;
t123 = t52 * t78;
t60 = t127 * t83;
t68 = t78 * pkin(4) + qJ(3);
t113 = t82 * qJD(6);
t111 = qJ(3) * qJD(3);
t110 = -0.2e1 * pkin(1) * qJD(2);
t42 = pkin(4) * t121 + t60;
t109 = pkin(5) * t72;
t108 = pkin(7) * t114;
t106 = t78 * t114;
t104 = t81 * t126;
t102 = qJ(6) * t72;
t97 = -t83 * pkin(2) - t116;
t6 = t19 * t79 + t20 * t78;
t49 = t125 * t78 + t81 * t79;
t24 = t49 * t114 - t44 * t83;
t38 = t49 * t83;
t96 = -t48 * t24 - t38 * t43;
t95 = -t49 * t44 + t124;
t93 = t126 * t103;
t54 = t126 * t78;
t16 = t49 * qJD(4) - qJD(5) * t93 + t54 * t115;
t34 = t79 * t104 + t125 * t54;
t91 = t16 * t82 - t34 * t72;
t17 = t54 * t100 - qJD(4) * t119 + (t125 * qJD(4) + qJD(5) * t104) * t79;
t33 = t81 * t54 - t93;
t90 = -t17 * t82 - t33 * t72;
t27 = t43 * t82 - t48 * t72;
t28 = t44 * t82 + t49 * t72;
t89 = t125 * t22 - t81 * t26;
t3 = -t22 * t100 - t81 * t11 + t26 * t115 - t125 * t14;
t36 = (-pkin(4) * t79 - t127) * t114;
t87 = pkin(5) * t43 + t44 * qJ(6) + t49 * qJD(6);
t1 = t102 - t3 + t113;
t2 = -t109 - t4;
t7 = t82 * qJ(6) + t131;
t8 = -t82 * pkin(5) - t89;
t86 = t1 * t49 + t2 * t48 - t8 * t43 + t7 * t44;
t85 = t16 * t49 - t17 * t48 + t33 * t43 - t34 * t44;
t84 = t97 * qJD(2) + t112;
t62 = 0.2e1 * t82 * t72;
t56 = -pkin(1) + t97;
t40 = -qJ(3) * t72 + t99;
t37 = t48 * t83;
t30 = -t78 * t47 + t51;
t29 = t49 * pkin(5) + t48 * qJ(6) + t68;
t25 = qJD(5) * t38 + t103 * t114 - t81 * t106;
t15 = t44 * pkin(5) - t43 * qJ(6) + t48 * qJD(6) + qJD(3);
t13 = -t37 * pkin(5) + t38 * qJ(6) + t42;
t5 = -t25 * pkin(5) - t24 * qJ(6) + t38 * qJD(6) + t36;
t9 = [0, 0, 0, t62, 0.2e1 * (-t82 ^ 2 + t83 ^ 2) * qJD(2), 0, 0, 0, t82 * t110, t83 * t110, 0, -0.2e1 * t56 * t114 + 0.2e1 * t40 * t83, -0.2e1 * t40 * t82 - 0.2e1 * t56 * t72, 0.2e1 * t56 * t40, -0.2e1 * t52 * t121 + 0.2e1 * t19 * t82 + 0.2e1 * (-t60 * t79 * t82 + t30 * t83) * qJD(2), 0.2e1 * t83 * t123 - 0.2e1 * t20 * t82 + 0.2e1 * (t60 * t122 - t31 * t83) * qJD(2), 0.2e1 * (t19 * t78 - t20 * t79) * t83 + 0.2e1 * (-t30 * t78 + t31 * t79) * t114, 0.2e1 * t30 * t19 + 0.2e1 * t31 * t20 - 0.2e1 * t60 * t52, -0.2e1 * t38 * t24, 0.2e1 * t24 * t37 - 0.2e1 * t38 * t25, 0.2e1 * t24 * t82 - 0.2e1 * t38 * t72, 0.2e1 * t25 * t82 + 0.2e1 * t37 * t72, t62, -0.2e1 * t42 * t25 - 0.2e1 * t36 * t37 + 0.2e1 * t4 * t82 + 0.2e1 * t72 * t89, -0.2e1 * t131 * t72 + 0.2e1 * t42 * t24 + 0.2e1 * t3 * t82 - 0.2e1 * t36 * t38, -0.2e1 * t13 * t25 - 0.2e1 * t2 * t82 - 0.2e1 * t5 * t37 - 0.2e1 * t72 * t8, 0.2e1 * t1 * t37 - 0.2e1 * t2 * t38 + 0.2e1 * t8 * t24 + 0.2e1 * t7 * t25, 0.2e1 * t1 * t82 - 0.2e1 * t13 * t24 + 0.2e1 * t5 * t38 + 0.2e1 * t7 * t72, 0.2e1 * t7 * t1 + 0.2e1 * t13 * t5 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, t72, -t114, 0, -t70, t108, t84, t70, -t108, t84 * pkin(7), -t130 * t79 - t123, t130 * t78 - t52 * t79, -t6, -t52 * qJ(3) + t60 * qJD(3) + t6 * t80 + (-t30 * t79 - t31 * t78) * qJD(4), t96, -t24 * t49 - t48 * t25 + t43 * t37 + t38 * t44, t27, -t28, 0, -qJD(3) * t37 - t68 * t25 + t36 * t49 + t42 * t44 + t90, -qJD(3) * t38 + t68 * t24 - t36 * t48 + t42 * t43 + t91, t13 * t44 - t15 * t37 - t29 * t25 + t5 * t49 + t90, -t16 * t37 - t17 * t38 + t33 * t24 + t34 * t25 - t86, -t13 * t43 + t15 * t38 - t29 * t24 + t5 * t48 - t91, t1 * t34 + t13 * t15 - t7 * t16 + t8 * t17 + t2 * t33 + t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0.2e1 * t111, t78 * t129, t79 * t129, 0.2e1 * t55, -0.2e1 * t80 * t55 + 0.2e1 * t111, -0.2e1 * t124, -0.2e1 * t43 * t49 + 0.2e1 * t48 * t44, 0, 0, 0, 0.2e1 * qJD(3) * t49 + 0.2e1 * t68 * t44, -0.2e1 * qJD(3) * t48 + 0.2e1 * t68 * t43, 0.2e1 * t15 * t49 + 0.2e1 * t29 * t44, 0.2e1 * t85, 0.2e1 * t15 * t48 - 0.2e1 * t29 * t43, 0.2e1 * t29 * t15 - 0.2e1 * t34 * t16 + 0.2e1 * t33 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, t70, t79 * t72, -t78 * t72, 0, t6, 0, 0, 0, 0, 0, t27, -t28, t27, t49 * t25 + t44 * t37 - t96, t28, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t95, 0, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t106, 0, -t52, 0, 0, 0, 0, 0, -t25, t24, -t25, 0, -t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t44, t43, t44, 0, -t43, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, t72, t4, t3, t4 + 0.2e1 * t109, -pkin(5) * t24 + t25 * qJ(6) + t37 * qJD(6), 0.2e1 * t102 - t3 + 0.2e1 * t113, -t2 * pkin(5) + t1 * qJ(6) + t7 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, 0, -t17, t16, -t17, -t87, -t16, -t17 * pkin(5) - t16 * qJ(6) + t34 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, t43, 0, t44, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, qJ(6) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
