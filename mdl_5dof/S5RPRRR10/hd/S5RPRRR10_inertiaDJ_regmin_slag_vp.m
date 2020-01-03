% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR10
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
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:47
% DurationCPUTime: 0.81s
% Computational Cost: add. (1221->136), mult. (2968->256), div. (0->0), fcn. (2915->8), ass. (0->101)
t116 = cos(qJ(3));
t64 = sin(pkin(9));
t65 = cos(pkin(9));
t68 = sin(qJ(3));
t44 = t116 * t64 + t68 * t65;
t66 = sin(qJ(5));
t67 = sin(qJ(4));
t69 = cos(qJ(5));
t70 = cos(qJ(4));
t46 = t66 * t70 + t69 * t67;
t23 = t46 * t44;
t57 = -t65 * pkin(2) - pkin(1);
t89 = t116 * t65;
t71 = -t68 * t64 + t89;
t29 = -pkin(3) * t71 - t44 * pkin(7) + t57;
t105 = pkin(6) + qJ(2);
t49 = t105 * t64;
t50 = t105 * t65;
t35 = t116 * t50 - t68 * t49;
t30 = t70 * t35;
t104 = t67 * t29 + t30;
t102 = qJD(4) * t67;
t38 = t71 * qJD(3);
t107 = t70 * t38;
t73 = -t44 * t102 + t107;
t63 = t70 ^ 2;
t103 = t67 ^ 2 - t63;
t83 = t103 * qJD(4);
t98 = qJD(4) + qJD(5);
t121 = 0.2e1 * t57;
t120 = pkin(7) + pkin(8);
t39 = t44 * qJD(3);
t119 = t39 * pkin(4);
t118 = t71 * pkin(4);
t101 = qJD(4) * t70;
t90 = t116 * t49;
t17 = qJD(3) * t90 - qJD(2) * t89 + (qJD(2) * t64 + qJD(3) * t50) * t68;
t28 = t39 * pkin(3) - t38 * pkin(7);
t6 = -t29 * t101 + t35 * t102 + t70 * t17 - t67 * t28;
t110 = t67 * t38;
t74 = t44 * t101 + t110;
t5 = -t74 * pkin(8) - t6;
t117 = t69 * t5;
t115 = t44 * t38;
t114 = t44 * t67;
t113 = t44 * t70;
t11 = -pkin(8) * t114 + t104;
t112 = t66 * t11;
t111 = t66 * t67;
t109 = t67 * t39;
t108 = t69 * t11;
t106 = t70 * t39;
t100 = qJD(5) * t66;
t99 = qJD(5) * t69;
t97 = -0.2e1 * pkin(3) * qJD(4);
t96 = pkin(4) * t102;
t95 = pkin(4) * t100;
t94 = pkin(4) * t99;
t92 = t67 * t101;
t85 = t67 * t17 + t70 * t28;
t4 = -pkin(8) * t107 + t119 + (-t30 + (pkin(8) * t44 - t29) * t67) * qJD(4) + t85;
t91 = t69 * t4 - t66 * t5;
t84 = t70 * t29 - t67 * t35;
t10 = -pkin(8) * t113 - t118 + t84;
t88 = -t10 + t118;
t87 = qJD(4) * t120;
t86 = -0.4e1 * t67 * t113;
t82 = 0.2e1 * (t64 ^ 2 + t65 ^ 2) * qJD(2);
t81 = -pkin(3) * t38 - pkin(7) * t39;
t80 = pkin(3) * t44 - pkin(7) * t71;
t34 = t68 * t50 + t90;
t79 = t69 * t10 - t112;
t78 = t66 * t10 + t108;
t32 = -t69 * t101 + t98 * t111 - t70 * t99;
t77 = -t32 * t71 - t46 * t39;
t51 = t120 * t67;
t52 = t120 * t70;
t76 = -t69 * t51 - t66 * t52;
t75 = -t66 * t51 + t69 * t52;
t45 = -t69 * t70 + t111;
t72 = -t101 * t71 + t109;
t18 = t44 * qJD(2) + t35 * qJD(3);
t59 = -t70 * pkin(4) - pkin(3);
t48 = t70 * t87;
t47 = t67 * t87;
t41 = t44 ^ 2;
t33 = t98 * t46;
t31 = -0.2e1 * t71 * t39;
t27 = t102 * t71 + t106;
t24 = t45 * t44;
t19 = pkin(4) * t114 + t34;
t15 = -t75 * qJD(5) + t66 * t47 - t69 * t48;
t14 = -t76 * qJD(5) + t69 * t47 + t66 * t48;
t13 = t33 * t71 - t45 * t39;
t12 = t74 * pkin(4) + t18;
t9 = -t100 * t114 + (t98 * t113 + t110) * t69 + t73 * t66;
t8 = -t98 * t23 - t45 * t38;
t7 = -t104 * qJD(4) + t85;
t2 = -t78 * qJD(5) + t91;
t1 = -t79 * qJD(5) - t66 * t4 - t117;
t3 = [0, 0, 0, 0, 0, t82, qJ(2) * t82, 0.2e1 * t115, 0.2e1 * t38 * t71 - 0.2e1 * t44 * t39, 0, 0, 0, t39 * t121, t38 * t121, 0.2e1 * t63 * t115 - 0.2e1 * t41 * t92, t38 * t86 + 0.2e1 * t41 * t83, 0.2e1 * t44 * t106 - 0.2e1 * t71 * t73, -0.2e1 * t44 * t109 + 0.2e1 * t71 * t74, t31, 0.2e1 * t18 * t114 + 0.2e1 * t74 * t34 + 0.2e1 * t84 * t39 - 0.2e1 * t7 * t71, -0.2e1 * t104 * t39 + 0.2e1 * t18 * t113 + 0.2e1 * t73 * t34 - 0.2e1 * t6 * t71, -0.2e1 * t24 * t8, -0.2e1 * t8 * t23 + 0.2e1 * t24 * t9, -0.2e1 * t24 * t39 - 0.2e1 * t71 * t8, -0.2e1 * t23 * t39 + 0.2e1 * t71 * t9, t31, 0.2e1 * t12 * t23 + 0.2e1 * t19 * t9 - 0.2e1 * t2 * t71 + 0.2e1 * t79 * t39, -0.2e1 * t1 * t71 - 0.2e1 * t12 * t24 + 0.2e1 * t19 * t8 - 0.2e1 * t78 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t38, 0, 0, 0, 0, 0, t27, -t72, 0, 0, 0, 0, 0, t13, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, -t18, t17, t67 * t107 - t44 * t83, qJD(4) * t86 - t103 * t38, t72, t27, 0, -t18 * t70 + t81 * t67 + (t34 * t67 - t80 * t70) * qJD(4), t18 * t67 + t81 * t70 + (t34 * t70 + t80 * t67) * qJD(4), t24 * t32 + t8 * t46, t32 * t23 + t24 * t33 - t8 * t45 - t46 * t9, -t77, t13, 0, t12 * t45 - t15 * t71 + t19 * t33 + t23 * t96 + t76 * t39 + t59 * t9, t12 * t46 - t14 * t71 - t19 * t32 - t24 * t96 - t75 * t39 + t59 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t92, -0.2e1 * t83, 0, 0, 0, t67 * t97, t70 * t97, -0.2e1 * t46 * t32, 0.2e1 * t32 * t45 - 0.2e1 * t46 * t33, 0, 0, 0, 0.2e1 * t59 * t33 + 0.2e1 * t45 * t96, -0.2e1 * t59 * t32 + 0.2e1 * t46 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t74, t39, t7, t6, 0, 0, t8, -t9, t39, t69 * t119 + (t88 * t66 - t108) * qJD(5) + t91, -t117 + (-t4 - t119) * t66 + (t88 * t69 + t112) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t101, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t102, 0, -pkin(7) * t101, pkin(7) * t102, 0, 0, -t32, -t33, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t39, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
