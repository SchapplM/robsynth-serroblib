% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:37
% EndTime: 2019-03-09 03:09:40
% DurationCPUTime: 1.20s
% Computational Cost: add. (1704->188), mult. (3984->330), div. (0->0), fcn. (3512->8), ass. (0->94)
t124 = cos(qJ(5));
t68 = -cos(pkin(9)) * pkin(1) - pkin(2);
t79 = sin(qJ(3));
t80 = cos(qJ(3));
t90 = -t80 * pkin(3) - t79 * qJ(4);
t50 = t68 + t90;
t76 = cos(pkin(10));
t45 = t76 * t50;
t67 = sin(pkin(9)) * pkin(1) + pkin(7);
t75 = sin(pkin(10));
t22 = -t76 * t79 * pkin(8) + t45 + (-t67 * t75 - pkin(4)) * t80;
t122 = t75 * t79;
t120 = t76 * t80;
t55 = t67 * t120;
t34 = t75 * t50 + t55;
t26 = -pkin(8) * t122 + t34;
t78 = sin(qJ(5));
t129 = t124 * t26 + t78 * t22;
t112 = qJD(5) * t78;
t96 = qJD(5) * t124;
t48 = t75 * t112 - t76 * t96;
t101 = t124 * t76;
t128 = -t78 * t75 + t101;
t70 = t79 * qJD(3);
t102 = t67 * t70;
t125 = pkin(3) * t79;
t47 = -t79 * qJD(4) + (-qJ(4) * t80 + t125) * qJD(3);
t31 = t75 * t102 + t76 * t47;
t15 = (pkin(4) * t79 - pkin(8) * t120) * qJD(3) + t31;
t121 = t75 * t80;
t39 = t75 * t47;
t60 = t79 * t67;
t21 = t39 + (-pkin(8) * t121 - t76 * t60) * qJD(3);
t4 = -qJD(5) * t129 + t124 * t15 - t78 * t21;
t127 = -0.2e1 * t48;
t126 = 2 * qJD(6);
t110 = t80 * qJD(3);
t104 = t75 * t110;
t54 = t124 * t75 + t78 * t76;
t49 = t54 * qJD(5);
t24 = -t101 * t110 + t78 * t104 + t79 * t49;
t43 = t128 * t79;
t123 = t43 * t24;
t118 = pkin(8) + qJ(4);
t117 = -t128 * t24 - t43 * t49;
t56 = t67 * t110;
t41 = pkin(4) * t104 + t56;
t46 = pkin(4) * t122 + t60;
t115 = t75 ^ 2 + t76 ^ 2;
t114 = t79 ^ 2 - t80 ^ 2;
t113 = qJD(4) * t75;
t111 = t76 * qJD(4);
t109 = t80 * qJD(6);
t108 = t67 * t121;
t107 = 0.2e1 * qJD(3) * t68;
t106 = pkin(5) * t70;
t103 = t79 * t110;
t69 = -t76 * pkin(4) - pkin(3);
t100 = qJ(6) * t70;
t99 = t118 * t75;
t98 = t115 * t80;
t95 = t124 * qJD(4);
t94 = t115 * qJD(4);
t93 = 0.2e1 * t103;
t92 = 0.2e1 * t94;
t25 = t54 * t110 - t48 * t79;
t42 = t54 * t79;
t89 = -t25 * t54 + t42 * t48;
t32 = -t76 * t102 + t39;
t88 = -t31 * t75 + t32 * t76;
t33 = t45 - t108;
t87 = -t33 * t75 + t34 * t76;
t86 = t124 * t99;
t59 = t118 * t76;
t19 = qJD(5) * t86 - t76 * t95 + (qJD(5) * t59 + t113) * t78;
t36 = t124 * t59 - t78 * t99;
t85 = -t19 * t80 - t36 * t70;
t20 = t59 * t96 + t78 * t111 + (-t118 * t112 + t95) * t75;
t35 = t78 * t59 + t86;
t84 = t20 * t80 - t35 * t70;
t28 = t80 * t48 + t54 * t70;
t27 = -t128 * t70 - t80 * t49;
t83 = t124 * t22 - t78 * t26;
t3 = t26 * t112 - t124 * t21 - t78 * t15 - t22 * t96;
t81 = -t25 * pkin(5) - t24 * qJ(6) + t43 * qJD(6);
t29 = -pkin(5) * t128 - t54 * qJ(6) + t69;
t14 = t49 * pkin(5) + t48 * qJ(6) - t54 * qJD(6);
t10 = t42 * pkin(5) - t43 * qJ(6) + t46;
t7 = t80 * pkin(5) - t83;
t6 = -t80 * qJ(6) + t129;
t5 = -t81 + t41;
t2 = -t106 - t4;
t1 = t100 - t3 - t109;
t8 = [0, 0, 0, 0, t93, -0.2e1 * t114 * qJD(3), 0, 0, 0, t79 * t107, t80 * t107, -0.2e1 * t31 * t80 + 0.2e1 * (t33 + 0.2e1 * t108) * t70, 0.2e1 * t32 * t80 + 0.2e1 * (-t34 + 0.2e1 * t55) * t70, 0.2e1 * (-t31 * t76 - t32 * t75) * t79 + 0.2e1 * (-t33 * t76 - t34 * t75) * t110, 0.2e1 * t67 ^ 2 * t103 + 0.2e1 * t33 * t31 + 0.2e1 * t34 * t32, -0.2e1 * t123, 0.2e1 * t42 * t24 - 0.2e1 * t25 * t43, 0.2e1 * t80 * t24 + 0.2e1 * t43 * t70, 0.2e1 * t80 * t25 - 0.2e1 * t42 * t70, -0.2e1 * t103, 0.2e1 * t46 * t25 - 0.2e1 * t4 * t80 + 0.2e1 * t41 * t42 + 0.2e1 * t70 * t83, -0.2e1 * t129 * t70 - 0.2e1 * t46 * t24 - 0.2e1 * t3 * t80 + 0.2e1 * t41 * t43, 0.2e1 * t10 * t25 + 0.2e1 * t2 * t80 + 0.2e1 * t5 * t42 - 0.2e1 * t7 * t70, -0.2e1 * t1 * t42 + 0.2e1 * t2 * t43 - 0.2e1 * t7 * t24 - 0.2e1 * t6 * t25, -0.2e1 * t1 * t80 + 0.2e1 * t10 * t24 - 0.2e1 * t5 * t43 + 0.2e1 * t6 * t70, 0.2e1 * t6 * t1 + 0.2e1 * t10 * t5 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t79 + (t114 * t67 + t87 * t80) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t43 + t10 * t70 + t2 * t42 - t6 * t24 + t7 * t25 - t5 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t115) * t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 * t25 - 0.2e1 * t103 - 0.2e1 * t123; 0, 0, 0, 0, 0, 0, t110, -t70, 0, -t56, t102, t80 * t113 + (t90 * t75 - t55) * qJD(3), t80 * t111 + (t90 * t76 + t108) * qJD(3), t88, -pkin(3) * t56 + qJ(4) * t88 + qJD(4) * t87, -t24 * t54 - t43 * t48, t89 + t117, t28, -t27, 0, -t128 * t41 + t69 * t25 + t46 * t49 + t84, -t69 * t24 + t41 * t54 - t46 * t48 + t85, t10 * t49 - t128 * t5 + t14 * t42 + t29 * t25 + t84, t1 * t128 + t19 * t42 + t2 * t54 + t20 * t43 - t35 * t24 - t36 * t25 - t7 * t48 - t6 * t49, t10 * t48 - t14 * t43 + t29 * t24 - t5 * t54 - t85, t1 * t36 + t10 * t14 - t6 * t19 + t2 * t35 + t7 * t20 + t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t110, -t76 * t70, t75 * t70, qJD(3) * t98, t79 * t94 + (qJ(4) * t98 - t125) * qJD(3), 0, 0, 0, 0, 0, t27, t28, t27, -t89 + t117, -t28, -t80 * t14 - t43 * t19 + t42 * t20 - t24 * t36 + t25 * t35 + t29 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, qJ(4) * t92, t54 * t127, -0.2e1 * t128 * t48 - 0.2e1 * t54 * t49, 0, 0, 0, 0.2e1 * t69 * t49, t69 * t127, -0.2e1 * t128 * t14 + 0.2e1 * t29 * t49, -0.2e1 * t128 * t19 + 0.2e1 * t20 * t54 - 0.2e1 * t35 * t48 - 0.2e1 * t36 * t49, -0.2e1 * t14 * t54 + 0.2e1 * t29 * t48, 0.2e1 * t29 * t14 - 0.2e1 * t36 * t19 + 0.2e1 * t35 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t76 * t110, 0, t56, 0, 0, 0, 0, 0, t25, -t24, t25, 0, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, t49, 0, t48, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, t70, t4, t3, t4 + 0.2e1 * t106, pkin(5) * t24 - t25 * qJ(6) - t42 * qJD(6), 0.2e1 * t100 - t3 - 0.2e1 * t109, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, -t25, 0, -t24, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49, 0, -t20, t19, -t20, pkin(5) * t48 - t49 * qJ(6) + qJD(6) * t128, -t19, -t20 * pkin(5) - t19 * qJ(6) + t36 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, qJ(6) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
