% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:04
% DurationCPUTime: 1.38s
% Computational Cost: add. (1210->149), mult. (2677->295), div. (0->0), fcn. (2208->6), ass. (0->97)
t118 = pkin(3) + pkin(6);
t109 = pkin(2) + qJ(4);
t97 = pkin(7) + t109;
t61 = cos(qJ(2));
t124 = t109 * t61;
t106 = qJ(3) * t61;
t55 = sin(pkin(8));
t59 = sin(qJ(2));
t56 = cos(pkin(8));
t73 = t55 * qJ(3) + t56 * t118;
t69 = pkin(4) + t73;
t122 = t69 * t59 + (t97 * t61 + pkin(1)) * t55;
t101 = t59 * qJD(3);
t75 = -t61 * qJD(4) - t101;
t71 = t56 * t75;
t94 = t55 * t118;
t123 = -t71 - (t61 * t94 + (t97 * t59 - t106) * t56) * qJD(2) - qJD(5) * t122;
t60 = cos(qJ(5));
t103 = qJD(5) * t60;
t58 = sin(qJ(5));
t104 = qJD(5) * t58;
t29 = t56 * t103 - t55 * t104;
t111 = t58 * t55;
t31 = t60 * t56 - t111;
t121 = qJD(2) * (t59 ^ 2 - t61 ^ 2);
t51 = t55 ^ 2;
t52 = t56 ^ 2;
t36 = (t51 + t52) * qJD(4);
t100 = t61 * qJD(3);
t105 = t59 * qJ(3);
t120 = (t105 + t124) * qJD(2) + qJD(4) * t59 - t100;
t119 = 0.2e1 * qJD(3);
t30 = t60 * t55 + t58 * t56;
t116 = t30 * t29;
t28 = t30 * qJD(5);
t115 = t31 * t28;
t102 = t59 * qJD(2);
t34 = t118 * t102;
t114 = t34 * t55;
t113 = t55 * t59;
t112 = t56 * t61;
t82 = -pkin(1) - t124;
t18 = t56 * (t82 - t105) + t59 * t94;
t40 = t118 * t61;
t48 = t61 * qJD(2);
t99 = qJ(3) * qJD(3);
t98 = -0.2e1 * pkin(1) * qJD(2);
t96 = pkin(6) * t102;
t95 = pkin(6) * t48;
t91 = t55 * t102;
t90 = t55 * t48;
t89 = t56 * t102;
t88 = t59 * t48;
t87 = t97 * t56;
t86 = t109 * t59;
t43 = -0.2e1 * t88;
t85 = t55 * t89;
t84 = t60 * t87;
t83 = 0.2e1 * t121;
t81 = -t61 * pkin(2) - t105;
t72 = t55 * t75;
t11 = -t72 + (-t55 * t86 + t73 * t61) * qJD(2);
t12 = t71 + (t56 * t86 + (-t56 * qJ(3) + t94) * t61) * qJD(2);
t3 = t11 * t56 + t12 * t55;
t13 = t30 * t102 - t29 * t61;
t23 = t30 * t61;
t80 = t13 * t31 + t28 * t23;
t14 = t61 * t28 - t58 * t91 + t60 * t89;
t22 = t31 * t61;
t79 = -t30 * t14 + t29 * t22;
t78 = t115 - t116;
t74 = -t29 * t59 - t30 * t48;
t15 = -pkin(7) * t112 + t18;
t62 = -t72 + (-t97 * t113 + t69 * t61) * qJD(2);
t1 = t15 * t104 + t123 * t60 - t58 * t62;
t2 = -t15 * t103 + t123 * t58 + t60 * t62;
t4 = t122 * t60 - t58 * t15;
t5 = t122 * t58 + t60 * t15;
t67 = -t1 * t30 + t2 * t31 - t4 * t28 + t5 * t29;
t35 = t97 * t55;
t19 = t58 * t35 - t84;
t20 = -t60 * t35 - t58 * t87;
t8 = t30 * qJD(4) + qJD(5) * t84 - t35 * t104;
t9 = qJD(4) * t111 + t35 * t103 + (-t60 * qJD(4) + t104 * t97) * t56;
t66 = t19 * t28 - t20 * t29 + t8 * t30 - t9 * t31;
t65 = t81 * qJD(2) + t100;
t46 = t55 * pkin(4) + qJ(3);
t44 = t56 * t48;
t42 = 0.2e1 * t88;
t37 = -pkin(1) + t81;
t33 = -0.2e1 * t121;
t27 = pkin(4) * t112 + t40;
t25 = -t101 + (pkin(2) * t59 - t106) * qJD(2);
t21 = (-pkin(4) * t56 - t118) * t102;
t17 = -t55 * t82 + t73 * t59;
t16 = -t28 * t59 + t31 * t48;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t33, 0, t43, 0, 0, t59 * t98, t61 * t98, 0, 0, 0, 0, 0, t42, t33, t43, 0, -0.2e1 * t37 * t102 + 0.2e1 * t25 * t61, -0.2e1 * t25 * t59 - 0.2e1 * t37 * t48, 0.2e1 * t37 * t25, t51 * t43, -0.4e1 * t61 * t85, t55 * t83, t52 * t43, t56 * t83, t42, -0.2e1 * t34 * t112 + 0.2e1 * t11 * t59 + 0.2e1 * (-t40 * t56 * t59 + t17 * t61) * qJD(2), 0.2e1 * t61 * t114 - 0.2e1 * t12 * t59 + 0.2e1 * (t40 * t113 - t18 * t61) * qJD(2), 0.2e1 * (t11 * t55 - t12 * t56) * t61 + 0.2e1 * (-t17 * t55 + t18 * t56) * t102, 0.2e1 * t17 * t11 + 0.2e1 * t18 * t12 - 0.2e1 * t40 * t34, -0.2e1 * t23 * t13, -0.2e1 * t13 * t22 - 0.2e1 * t23 * t14, 0.2e1 * t13 * t59 - 0.2e1 * t23 * t48, -0.2e1 * t22 * t14, 0.2e1 * t14 * t59 - 0.2e1 * t22 * t48, t42, -0.2e1 * t27 * t14 + 0.2e1 * t2 * t59 + 0.2e1 * t21 * t22 + 0.2e1 * t4 * t48, 0.2e1 * t1 * t59 + 0.2e1 * t27 * t13 - 0.2e1 * t21 * t23 - 0.2e1 * t5 * t48, 0.2e1 * t1 * t22 - 0.2e1 * t4 * t13 + 0.2e1 * t5 * t14 + 0.2e1 * t2 * t23, -0.2e1 * t5 * t1 + 0.2e1 * t4 * t2 + 0.2e1 * t27 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t102, 0, -t95, t96, 0, 0, 0, -t48, t102, 0, 0, 0, t65, t95, -t96, t65 * pkin(6), t85, (-t51 + t52) * t102, t44, -t85, -t90, 0, -t120 * t56 - t114, t120 * t55 - t34 * t56, -t3, -t34 * qJ(3) + t40 * qJD(3) - t3 * t109 + (-t17 * t56 - t18 * t55) * qJD(4), t80, -t13 * t30 + t31 * t14 + t28 * t22 + t23 * t29, t16, t79, t74, 0, qJD(3) * t22 - t46 * t14 + t19 * t48 + t21 * t30 + t27 * t29 + t9 * t59, -qJD(3) * t23 + t46 * t13 - t20 * t48 + t21 * t31 - t27 * t28 + t8 * t59, -t19 * t13 + t20 * t14 + t8 * t22 + t9 * t23 - t67, t27 * qJD(3) - t1 * t20 + t2 * t19 + t21 * t46 + t4 * t9 - t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, 0.2e1 * t99, 0, 0, 0, 0, 0, 0, t55 * t119, t56 * t119, 0.2e1 * t36, 0.2e1 * t109 * t36 + 0.2e1 * t99, -0.2e1 * t115, 0.2e1 * t28 * t30 - 0.2e1 * t31 * t29, 0, 0.2e1 * t116, 0, 0, 0.2e1 * qJD(3) * t30 + 0.2e1 * t46 * t29, 0.2e1 * qJD(3) * t31 - 0.2e1 * t46 * t28, 0.2e1 * t66, 0.2e1 * t46 * qJD(3) + 0.2e1 * t19 * t9 - 0.2e1 * t20 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t95, 0, 0, 0, 0, 0, 0, t44, -t90, 0, t3, 0, 0, 0, 0, 0, 0, t16, t74, -t79 - t80, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t91, 0, -t34, 0, 0, 0, 0, 0, 0, -t14, t13, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, 0, t29, -t28, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t14, t48, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t29, 0, t9, t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
