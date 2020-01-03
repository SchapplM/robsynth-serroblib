% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:05:58
% DurationCPUTime: 1.30s
% Computational Cost: add. (824->181), mult. (2042->315), div. (0->0), fcn. (1414->4), ass. (0->110)
t53 = cos(qJ(3));
t102 = t53 * qJD(4);
t51 = sin(qJ(3));
t111 = t51 * qJ(4);
t124 = pkin(3) + pkin(4);
t130 = t124 * t53 + t111;
t133 = t130 * qJD(3) - t102;
t52 = sin(qJ(2));
t103 = t52 * qJD(2);
t54 = cos(qJ(2));
t106 = qJD(3) * t54;
t90 = t53 * t106;
t132 = t51 * t103 - t90;
t131 = t124 * qJD(2);
t129 = -0.4e1 * t52;
t49 = t53 ^ 2;
t115 = t51 ^ 2 - t49;
t80 = t115 * qJD(3);
t104 = t51 * qJD(5);
t100 = t54 * qJD(2);
t82 = qJ(5) * t100;
t44 = qJD(3) * t53;
t91 = t52 * t44;
t128 = qJ(5) * t91 + t52 * t104 + t51 * t82;
t122 = pkin(7) * t54;
t74 = pkin(2) * t52 - t122;
t27 = t74 * qJD(2);
t119 = t52 * pkin(7);
t75 = -t54 * pkin(2) - t119;
t29 = -pkin(1) + t75;
t92 = t51 * t106;
t62 = t53 * t103 + t92;
t7 = t62 * pkin(6) - t51 * t27 - t29 * t44;
t113 = qJ(4) * t53;
t71 = pkin(3) * t51 - t113;
t68 = pkin(6) + t71;
t15 = t68 * t52;
t105 = t51 * qJD(4);
t17 = t71 * qJD(3) - t105;
t72 = t53 * pkin(3) + t111;
t28 = -pkin(2) - t72;
t127 = qJD(2) * (-t28 * t54 + t119) - qJD(3) * t15 - t17 * t52;
t125 = t72 * qJD(3) - t102;
t56 = 0.2e1 * qJD(4);
t123 = pkin(6) * t54;
t66 = -t124 * t51 + t113;
t63 = -pkin(6) + t66;
t3 = t63 * t100 - t133 * t52;
t121 = t3 * t51;
t120 = t3 * t53;
t118 = pkin(7) - qJ(5);
t40 = t53 * t123;
t116 = t51 * t29 + t40;
t48 = t52 ^ 2;
t114 = -t54 ^ 2 + t48;
t112 = qJ(5) * t52;
t110 = qJD(2) * t51;
t109 = qJD(2) * t53;
t108 = qJD(3) * t51;
t107 = qJD(3) * t52;
t101 = t53 * qJD(5);
t99 = t54 * qJD(4);
t98 = qJ(4) * qJD(2);
t39 = t51 * t123;
t97 = -0.2e1 * pkin(1) * qJD(2);
t96 = -0.2e1 * pkin(2) * qJD(3);
t95 = pkin(3) * t103;
t94 = pkin(7) * t108;
t93 = pkin(7) * t44;
t25 = pkin(2) + t130;
t89 = t25 * t44;
t87 = t51 * t44;
t86 = t52 * t100;
t85 = t53 * t100;
t84 = t54 * t98;
t83 = qJ(4) * t107;
t31 = t118 * t53;
t81 = t53 * t29 - t39;
t79 = t114 * qJD(2);
t78 = 0.2e1 * t86;
t77 = -t132 * pkin(6) + t29 * t108 - t53 * t27;
t76 = t51 * t85;
t12 = -t54 * qJ(4) + t116;
t45 = t54 * pkin(3);
t13 = t45 - t81;
t70 = -t12 * t51 + t13 * t53;
t14 = t66 * qJD(3) + t105;
t67 = -t25 * t108 + t14 * t53;
t11 = t63 * t52;
t65 = qJD(3) * t11 + t25 * t100;
t64 = -qJ(5) * t108 + t101;
t61 = -t53 * t82 + t77;
t42 = t52 * t98;
t4 = t42 - t7 - t99;
t5 = t77 - t95;
t58 = t70 * qJD(3) + t4 * t53 + t5 * t51;
t57 = 0.2e1 * t42 - 0.2e1 * t99 - t7;
t46 = qJ(4) * t56;
t36 = pkin(7) * t90;
t30 = t118 * t51;
t20 = -t51 * t100 - t91;
t19 = t51 * t107 - t85;
t18 = qJD(3) * t31 - t104;
t16 = -t118 * t108 - t101;
t10 = t51 * t112 + t12;
t9 = t54 * pkin(4) + t39 + t45 + (-t29 - t112) * t53;
t6 = t68 * t100 + t125 * t52;
t2 = t4 + t128;
t1 = (-t64 - t131) * t52 + t61;
t8 = [0, 0, 0, t78, -0.2e1 * t79, 0, 0, 0, t52 * t97, t54 * t97, -0.2e1 * t48 * t87 + 0.2e1 * t49 * t86, t76 * t129 + 0.2e1 * t48 * t80, 0.2e1 * t114 * t109 + 0.2e1 * t52 * t92, -0.2e1 * t51 * t79 + 0.2e1 * t52 * t90, -0.2e1 * t86, 0.2e1 * t77 * t54 + 0.2e1 * t81 * t103 + 0.2e1 * (t44 * t48 + t51 * t78) * pkin(6), -0.2e1 * t7 * t54 - 0.2e1 * t116 * t103 + 0.2e1 * (-t108 * t48 + t53 * t78) * pkin(6), 0.2e1 * (t15 * t110 + t5) * t54 + 0.2e1 * (-qJD(2) * t13 + t15 * t44 + t6 * t51) * t52, 0.2e1 * t70 * t100 + 0.2e1 * (-t4 * t51 + t5 * t53 + (-t12 * t53 - t13 * t51) * qJD(3)) * t52, 0.2e1 * (-t15 * t109 - t4) * t54 + 0.2e1 * (qJD(2) * t12 + t108 * t15 - t6 * t53) * t52, 0.2e1 * t12 * t4 + 0.2e1 * t13 * t5 + 0.2e1 * t15 * t6, 0.2e1 * (-t11 * t110 + t1) * t54 + 0.2e1 * (-qJD(2) * t9 - t11 * t44 - t121) * t52, 0.2e1 * (t11 * t109 - t2) * t54 + 0.2e1 * (qJD(2) * t10 - t108 * t11 + t120) * t52, 0.2e1 * (t10 * t51 - t53 * t9) * t100 + 0.2e1 * (-t1 * t53 + t2 * t51 + (t10 * t53 + t51 * t9) * qJD(3)) * t52, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, t100, -t103, 0, -pkin(6) * t100, pkin(6) * t103, -t52 * t80 + t76, -t115 * t100 + t87 * t129, t132, t62, 0, t36 + (-pkin(2) * t53 + pkin(6) * t51) * t107 + (t51 * t75 - t40) * qJD(2), (pkin(6) * t52 * t53 + t51 * t74) * qJD(3) + (t53 * t75 + t39) * qJD(2), t36 + (t107 * t28 - t6) * t53 - t127 * t51, t58, (-t6 + (t28 * t52 + t122) * qJD(3)) * t51 + t127 * t53, pkin(7) * t58 + t15 * t17 + t6 * t28, t18 * t54 + t120 + (-qJD(2) * t30 - t89) * t52 + (-t14 * t52 - t65) * t51, -t16 * t54 + t121 + t65 * t53 + (qJD(2) * t31 + t67) * t52, (-t30 * t100 - t18 * t52 - t2 + (t31 * t52 - t9) * qJD(3)) * t53 + (t31 * t100 + t16 * t52 - t1 + (t30 * t52 + t10) * qJD(3)) * t51, t1 * t30 + t10 * t16 + t11 * t14 + t9 * t18 + t2 * t31 + t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t80, 0, 0, 0, t51 * t96, t53 * t96, 0.2e1 * t108 * t28 - 0.2e1 * t17 * t53, 0, -0.2e1 * t17 * t51 - 0.2e1 * t28 * t44, 0.2e1 * t28 * t17, 0.2e1 * t67, 0.2e1 * t14 * t51 + 0.2e1 * t89, -0.2e1 * t16 * t53 - 0.2e1 * t18 * t51 + 0.2e1 * (-t30 * t53 + t31 * t51) * qJD(3), 0.2e1 * t25 * t14 + 0.2e1 * t31 * t16 + 0.2e1 * t30 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t20, t103, -t77, t7, -t77 + 0.2e1 * t95, (-pkin(3) * t100 - t83) * t53 + (-t84 + (pkin(3) * qJD(3) - qJD(4)) * t52) * t51, t57, -t5 * pkin(3) + t4 * qJ(4) + t12 * qJD(4), (t64 + 0.2e1 * t131) * t52 - t61, t57 + t128, (t100 * t124 + t83) * t53 + (t84 + (-qJD(3) * t124 + qJD(4)) * t52) * t51, t2 * qJ(4) + t10 * qJD(4) - t1 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t108, 0, -t93, t94, -t93, -t125, -t94, -t125 * pkin(7), -t18, t16, t133, t16 * qJ(4) + t31 * qJD(4) - t124 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t46, 0, t56, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t19, 0, t5, -t103, 0, t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t93, 0, 0, -t44, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t44, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
