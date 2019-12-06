% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:10
% EndTime: 2019-12-05 17:21:14
% DurationCPUTime: 0.97s
% Computational Cost: add. (735->154), mult. (2155->299), div. (0->0), fcn. (1951->10), ass. (0->114)
t59 = sin(qJ(3));
t127 = -0.4e1 * t59;
t57 = sin(qJ(5));
t58 = sin(qJ(4));
t61 = cos(qJ(5));
t62 = cos(qJ(4));
t35 = t57 * t62 + t58 * t61;
t25 = t35 * t59;
t63 = cos(qJ(3));
t77 = -pkin(3) * t63 - pkin(8) * t59;
t40 = -pkin(2) + t77;
t115 = t62 * t63;
t47 = pkin(7) * t115;
t112 = t58 * t40 + t47;
t106 = qJD(4) * t58;
t102 = t63 * qJD(3);
t85 = t62 * t102;
t126 = -t59 * t106 + t85;
t53 = t62 ^ 2;
t111 = t58 ^ 2 - t53;
t81 = t111 * qJD(4);
t125 = qJD(4) + qJD(5);
t105 = qJD(4) * t62;
t76 = pkin(3) * t59 - pkin(8) * t63;
t38 = t76 * qJD(3);
t50 = t59 * qJD(3);
t104 = qJD(4) * t63;
t93 = t58 * t104;
t66 = t50 * t62 + t93;
t14 = pkin(7) * t66 - t105 * t40 - t38 * t58;
t124 = pkin(8) + pkin(9);
t123 = pkin(7) * t58;
t88 = t58 * t102;
t65 = t105 * t59 + t88;
t9 = -pkin(9) * t65 - t14;
t122 = t61 * t9;
t55 = sin(pkin(5));
t60 = sin(qJ(2));
t121 = t55 * t60;
t64 = cos(qJ(2));
t120 = t55 * t64;
t118 = t58 * t59;
t23 = -pkin(9) * t118 + t112;
t119 = t57 * t23;
t117 = t59 * t62;
t116 = t61 * t23;
t89 = t58 * t50;
t113 = pkin(7) * t89 + t38 * t62;
t52 = t59 ^ 2;
t110 = -t63 ^ 2 + t52;
t109 = qJD(2) * t60;
t56 = cos(pkin(5));
t27 = t121 * t59 - t56 * t63;
t108 = qJD(3) * t27;
t107 = qJD(3) * t62;
t103 = qJD(5) * t57;
t101 = -0.2e1 * pkin(2) * qJD(3);
t100 = -0.2e1 * pkin(3) * qJD(4);
t99 = pkin(4) * t106;
t98 = pkin(4) * t50;
t97 = pkin(4) * t103;
t96 = qJD(5) * t61 * pkin(4);
t95 = pkin(7) * t102;
t92 = t62 * t104;
t91 = t55 * t109;
t90 = qJD(2) * t120;
t87 = t58 * t105;
t86 = t59 * t102;
t8 = (pkin(4) * t59 - pkin(9) * t115) * qJD(3) + (-t47 + (pkin(9) * t59 - t40) * t58) * qJD(4) + t113;
t84 = -t57 * t9 + t61 * t8;
t33 = t62 * t40;
t16 = -pkin(9) * t117 + t33 + (-pkin(4) - t123) * t63;
t83 = pkin(4) * t63 - t16;
t82 = qJD(4) * t124;
t80 = t110 * qJD(3);
t79 = 0.2e1 * t86;
t78 = t58 * t85;
t75 = t16 * t61 - t119;
t74 = t16 * t57 + t116;
t28 = t121 * t63 + t56 * t59;
t21 = -t120 * t62 - t28 * t58;
t69 = t120 * t58 - t28 * t62;
t73 = t21 * t61 + t57 * t69;
t72 = t21 * t57 - t61 * t69;
t43 = t124 * t58;
t44 = t124 * t62;
t71 = -t43 * t61 - t44 * t57;
t70 = -t43 * t57 + t44 * t61;
t34 = t57 * t58 - t61 * t62;
t20 = qJD(3) * t28 + t59 * t90;
t68 = t105 * t27 + t20 * t58;
t67 = t106 * t27 - t20 * t62;
t49 = -pkin(4) * t62 - pkin(3);
t46 = -0.2e1 * t86;
t39 = (pkin(4) * t58 + pkin(7)) * t59;
t37 = t62 * t82;
t36 = t58 * t82;
t26 = t34 * t59;
t24 = pkin(4) * t65 + t95;
t19 = t63 * t90 - t108;
t18 = t125 * t35;
t17 = t125 * t34;
t15 = -qJD(4) * t112 + t113;
t13 = -qJD(5) * t70 + t57 * t36 - t61 * t37;
t12 = -qJD(5) * t71 + t61 * t36 + t57 * t37;
t11 = -t103 * t118 + (t117 * t125 + t88) * t61 + t126 * t57;
t10 = -t34 * t102 - t125 * t25;
t6 = qJD(4) * t21 + t19 * t62 + t58 * t91;
t5 = qJD(4) * t69 - t19 * t58 + t62 * t91;
t4 = -qJD(5) * t74 + t84;
t3 = -qJD(5) * t75 - t57 * t8 - t122;
t2 = -qJD(5) * t72 + t61 * t5 - t57 * t6;
t1 = -qJD(5) * t73 - t57 * t5 - t61 * t6;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t91, -t90, 0, 0, 0, 0, 0, (-t109 * t63 - t50 * t64) * t55, (-t102 * t64 + t109 * t59) * t55, 0, 0, 0, 0, 0, (t108 * t58 - t5) * t63 + (qJD(3) * t21 + t68) * t59, (t107 * t27 + t6) * t63 + (qJD(3) * t69 - t67) * t59, 0, 0, 0, 0, 0, t11 * t27 - t2 * t63 + t20 * t25 + t50 * t73, -t1 * t63 + t10 * t27 - t20 * t26 - t50 * t72; 0, 0, 0, 0, t79, -0.2e1 * t80, 0, 0, 0, t59 * t101, t63 * t101, -0.2e1 * t52 * t87 + 0.2e1 * t53 * t86, t127 * t78 + 0.2e1 * t52 * t81, 0.2e1 * t107 * t110 + 0.2e1 * t59 * t93, -0.2e1 * t58 * t80 + 0.2e1 * t59 * t92, t46, 0.2e1 * t33 * t50 - 0.2e1 * t15 * t63 + 0.2e1 * (t105 * t52 + t58 * t86) * pkin(7), -0.2e1 * t14 * t63 - 0.2e1 * t112 * t50 + 0.2e1 * (-t106 * t52 + t62 * t79) * pkin(7), -0.2e1 * t26 * t10, -0.2e1 * t10 * t25 + 0.2e1 * t11 * t26, -0.2e1 * t10 * t63 - 0.2e1 * t26 * t50, 0.2e1 * t11 * t63 - 0.2e1 * t25 * t50, t46, 0.2e1 * t11 * t39 + 0.2e1 * t24 * t25 - 0.2e1 * t4 * t63 + 0.2e1 * t50 * t75, 0.2e1 * t10 * t39 - 0.2e1 * t24 * t26 - 0.2e1 * t3 * t63 - 0.2e1 * t50 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t19, 0, 0, 0, 0, 0, t67, t68, 0, 0, 0, 0, 0, t18 * t27 + t20 * t34, -t17 * t27 + t20 * t35; 0, 0, 0, 0, 0, 0, t102, -t50, 0, -t95, pkin(7) * t50, -t59 * t81 + t78, -t102 * t111 + t127 * t87, t89 - t92, t66, 0, (pkin(8) * t115 + (-pkin(3) * t62 + t123) * t59) * qJD(4) + (t58 * t77 - t47) * qJD(3), (pkin(7) * t117 + t58 * t76) * qJD(4) + (t123 * t63 + t62 * t77) * qJD(3), t10 * t35 + t17 * t26, -t10 * t34 - t11 * t35 + t17 * t25 + t18 * t26, t17 * t63 + t35 * t50, t18 * t63 - t34 * t50, 0, t11 * t49 - t13 * t63 + t18 * t39 + t24 * t34 + t25 * t99 + t50 * t71, t10 * t49 - t12 * t63 - t17 * t39 + t24 * t35 - t26 * t99 - t50 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, -0.2e1 * t81, 0, 0, 0, t58 * t100, t62 * t100, -0.2e1 * t35 * t17, 0.2e1 * t17 * t34 - 0.2e1 * t18 * t35, 0, 0, 0, 0.2e1 * t18 * t49 + 0.2e1 * t34 * t99, -0.2e1 * t17 * t49 + 0.2e1 * t35 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t65, t50, t15, t14, 0, 0, t10, -t11, t50, t61 * t98 + (t57 * t83 - t116) * qJD(5) + t84, -t122 + (-t8 - t98) * t57 + (t61 * t83 + t119) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t106, 0, -pkin(8) * t105, pkin(8) * t106, 0, 0, -t17, -t18, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t97, -0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t50, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
