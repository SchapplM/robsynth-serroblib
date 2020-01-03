% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:22:01
% DurationCPUTime: 0.83s
% Computational Cost: add. (1313->146), mult. (3116->283), div. (0->0), fcn. (2966->8), ass. (0->106)
t65 = sin(pkin(9));
t66 = cos(pkin(9));
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t48 = t65 * t72 + t66 * t69;
t67 = sin(qJ(5));
t68 = sin(qJ(4));
t70 = cos(qJ(5));
t71 = cos(qJ(4));
t51 = t67 * t71 + t70 * t68;
t24 = t51 * t48;
t47 = t65 * t69 - t66 * t72;
t92 = -t72 * pkin(2) - pkin(1);
t29 = t47 * pkin(3) - t48 * pkin(7) + t92;
t109 = -qJ(3) - pkin(6);
t53 = t109 * t69;
t54 = t109 * t72;
t35 = t65 * t53 - t66 * t54;
t31 = t71 * t35;
t108 = t68 * t29 + t31;
t106 = qJD(4) * t68;
t101 = t72 * qJD(2);
t102 = t69 * qJD(2);
t43 = t66 * t101 - t65 * t102;
t110 = t71 * t43;
t74 = -t48 * t106 + t110;
t100 = qJD(4) + qJD(5);
t123 = 0.2e1 * qJD(4);
t42 = t48 * qJD(2);
t122 = t42 * pkin(4);
t121 = t47 * pkin(4);
t105 = qJD(4) * t71;
t85 = qJD(2) * t109;
t41 = t72 * qJD(3) + t69 * t85;
t73 = -t69 * qJD(3) + t72 * t85;
t21 = t66 * t41 + t65 * t73;
t62 = pkin(2) * t102;
t22 = t42 * pkin(3) - t43 * pkin(7) + t62;
t6 = -t29 * t105 + t35 * t106 - t71 * t21 - t68 * t22;
t113 = t68 * t43;
t75 = t48 * t105 + t113;
t5 = -t75 * pkin(8) - t6;
t120 = t70 * t5;
t59 = t65 * pkin(2) + pkin(7);
t119 = pkin(8) + t59;
t118 = t48 * t68;
t117 = t48 * t71;
t11 = -pkin(8) * t118 + t108;
t116 = t67 * t11;
t115 = t67 * t68;
t114 = t68 * t42;
t112 = t70 * t11;
t111 = t71 * t42;
t64 = t71 ^ 2;
t107 = t68 ^ 2 - t64;
t104 = qJD(5) * t67;
t103 = qJD(5) * t70;
t99 = -0.2e1 * pkin(1) * qJD(2);
t60 = -t66 * pkin(2) - pkin(3);
t98 = t60 * t123;
t97 = pkin(4) * t106;
t96 = pkin(4) * t104;
t95 = pkin(4) * t103;
t93 = t68 * t105;
t88 = -t68 * t21 + t71 * t22;
t4 = -pkin(8) * t110 + t122 + (-t31 + (pkin(8) * t48 - t29) * t68) * qJD(4) + t88;
t91 = t70 * t4 - t67 * t5;
t87 = t71 * t29 - t68 * t35;
t10 = -pkin(8) * t117 + t121 + t87;
t90 = -t10 - t121;
t89 = -0.4e1 * t68 * t117;
t20 = t65 * t41 - t66 * t73;
t34 = -t66 * t53 - t65 * t54;
t86 = qJD(4) * t119;
t84 = t107 * qJD(4);
t83 = t70 * t10 - t116;
t82 = t67 * t10 + t112;
t32 = t100 * t115 - t71 * t103 - t70 * t105;
t81 = t32 * t47 - t51 * t42;
t80 = -t42 * t59 + t43 * t60;
t44 = t119 * t68;
t45 = t119 * t71;
t79 = -t70 * t44 - t67 * t45;
t78 = -t67 * t44 + t70 * t45;
t77 = t47 * t59 - t48 * t60;
t50 = -t70 * t71 + t115;
t76 = t47 * t105 + t114;
t52 = -t71 * pkin(4) + t60;
t46 = t48 ^ 2;
t40 = t71 * t86;
t39 = t68 * t86;
t33 = t100 * t51;
t30 = 0.2e1 * t47 * t42;
t28 = -t47 * t106 + t111;
t25 = t50 * t48;
t19 = pkin(4) * t118 + t34;
t15 = -t33 * t47 - t50 * t42;
t14 = -t78 * qJD(5) + t67 * t39 - t70 * t40;
t13 = -t79 * qJD(5) + t70 * t39 + t67 * t40;
t12 = t75 * pkin(4) + t20;
t9 = -t104 * t118 + (t100 * t117 + t113) * t70 + t74 * t67;
t8 = -t100 * t24 - t50 * t43;
t7 = -t108 * qJD(4) + t88;
t2 = -t82 * qJD(5) + t91;
t1 = -t83 * qJD(5) - t67 * t4 - t120;
t3 = [0, 0, 0, 0.2e1 * t69 * t101, 0.2e1 * (-t69 ^ 2 + t72 ^ 2) * qJD(2), 0, 0, 0, t69 * t99, t72 * t99, 0.2e1 * t20 * t48 - 0.2e1 * t21 * t47 + 0.2e1 * t34 * t43 - 0.2e1 * t35 * t42, 0.2e1 * t34 * t20 + 0.2e1 * t35 * t21 + 0.2e1 * t92 * t62, 0.2e1 * t64 * t48 * t43 - 0.2e1 * t46 * t93, t107 * t46 * t123 + t43 * t89, 0.2e1 * t48 * t111 + 0.2e1 * t74 * t47, -0.2e1 * t48 * t114 - 0.2e1 * t75 * t47, t30, 0.2e1 * t20 * t118 + 0.2e1 * t75 * t34 + 0.2e1 * t87 * t42 + 0.2e1 * t7 * t47, -0.2e1 * t108 * t42 + 0.2e1 * t20 * t117 + 0.2e1 * t74 * t34 + 0.2e1 * t6 * t47, -0.2e1 * t25 * t8, -0.2e1 * t8 * t24 + 0.2e1 * t25 * t9, -0.2e1 * t25 * t42 + 0.2e1 * t8 * t47, -0.2e1 * t24 * t42 - 0.2e1 * t9 * t47, t30, 0.2e1 * t12 * t24 + 0.2e1 * t19 * t9 + 0.2e1 * t2 * t47 + 0.2e1 * t83 * t42, 0.2e1 * t1 * t47 - 0.2e1 * t12 * t25 + 0.2e1 * t19 * t8 - 0.2e1 * t82 * t42; 0, 0, 0, 0, 0, t101, -t102, 0, -pkin(6) * t101, pkin(6) * t102, (-t42 * t65 - t43 * t66) * pkin(2), (-t20 * t66 + t21 * t65) * pkin(2), t68 * t110 - t48 * t84, qJD(4) * t89 - t107 * t43, t76, t28, 0, -t20 * t71 + t80 * t68 + (t34 * t68 - t77 * t71) * qJD(4), t20 * t68 + t80 * t71 + (t34 * t71 + t77 * t68) * qJD(4), t25 * t32 + t8 * t51, t32 * t24 + t25 * t33 - t8 * t50 - t51 * t9, -t81, t15, 0, t12 * t50 + t14 * t47 + t19 * t33 + t24 * t97 + t79 * t42 + t52 * t9, t12 * t51 + t13 * t47 - t19 * t32 - t25 * t97 - t78 * t42 + t52 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t93, -0.2e1 * t84, 0, 0, 0, t68 * t98, t71 * t98, -0.2e1 * t51 * t32, 0.2e1 * t32 * t50 - 0.2e1 * t51 * t33, 0, 0, 0, 0.2e1 * t52 * t33 + 0.2e1 * t50 * t97, -0.2e1 * t52 * t32 + 0.2e1 * t51 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, t28, -t76, 0, 0, 0, 0, 0, t15, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t75, t42, t7, t6, 0, 0, t8, -t9, t42, t70 * t122 + (t90 * t67 - t112) * qJD(5) + t91, -t120 + (-t4 - t122) * t67 + (t90 * t70 + t116) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t106, 0, -t59 * t105, t59 * t106, 0, 0, -t32, -t33, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t96, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t42, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
