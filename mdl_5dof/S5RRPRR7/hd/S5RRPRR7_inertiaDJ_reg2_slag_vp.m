% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:43
% DurationCPUTime: 0.82s
% Computational Cost: add. (988->110), mult. (1913->172), div. (0->0), fcn. (1512->6), ass. (0->92)
t119 = sin(qJ(5));
t123 = qJD(4) + qJD(5);
t67 = sin(qJ(4));
t69 = cos(qJ(5));
t70 = cos(qJ(4));
t42 = t119 * t70 + t69 * t67;
t24 = t123 * t42;
t115 = t69 * t70;
t93 = qJD(5) * t119;
t97 = t119 * t67;
t25 = -qJD(4) * t97 + t123 * t115 - t67 * t93;
t43 = -t97 + t115;
t4 = ((t119 * t43 - t42 * t69) * qJD(5) - t119 * t25 + t24 * t69) * pkin(4);
t65 = t67 ^ 2;
t66 = t70 ^ 2;
t108 = t65 + t66;
t71 = cos(qJ(2));
t120 = t71 * pkin(1);
t98 = -pkin(2) - t120;
t89 = -pkin(7) + t98;
t87 = pkin(8) - t89;
t39 = t87 * t67;
t81 = t87 * t70;
t77 = t69 * t81;
t20 = t119 * t39 - t77;
t75 = t119 * t81;
t21 = -t69 * t39 - t75;
t105 = t67 * qJD(4);
t107 = pkin(1) * qJD(2);
t68 = sin(qJ(2));
t60 = t68 * t107;
t90 = t70 * t60;
t73 = t87 * t105 + t90;
t91 = t67 * t60;
t74 = -qJD(4) * t81 + t91;
t5 = qJD(5) * t77 - t119 * t73 - t39 * t93 - t69 * t74;
t106 = qJD(5) * t69;
t6 = qJD(5) * t75 + t39 * t106 - t119 * t74 + t69 * t73;
t95 = t20 * t24 - t21 * t25 + t5 * t42 - t6 * t43;
t121 = pkin(2) + pkin(7);
t102 = pkin(8) + t121;
t47 = t102 * t67;
t84 = t102 * t105;
t92 = t102 * t70;
t86 = t69 * t92;
t12 = -t119 * t84 + t123 * t86 - t47 * t93;
t80 = t119 * t92;
t13 = t47 * t106 + t123 * t80 + t69 * t84;
t28 = t119 * t47 - t86;
t29 = -t69 * t47 - t80;
t94 = t12 * t42 - t13 * t43 + t28 * t24 - t29 * t25;
t72 = 2 * qJD(3);
t118 = t42 * t25;
t117 = t43 * t24;
t101 = t71 * t107;
t51 = qJD(3) + t101;
t55 = t68 * pkin(1) + qJ(3);
t116 = t55 * t51;
t104 = t70 * qJD(4);
t59 = pkin(4) * t104;
t40 = t51 + t59;
t63 = t67 * pkin(4);
t48 = t55 + t63;
t114 = t48 * t25 + t40 * t42;
t113 = -t48 * t24 + t40 * t43;
t52 = qJD(3) + t59;
t56 = qJ(3) + t63;
t112 = t56 * t25 + t52 * t42;
t111 = -t56 * t24 + t52 * t43;
t110 = t55 * t104 + t51 * t67;
t103 = qJ(3) * qJD(4);
t109 = qJD(3) * t67 + t70 * t103;
t100 = pkin(4) * t106;
t99 = t67 * t104;
t96 = t121 * qJD(4);
t88 = pkin(4) * t93;
t85 = t117 - t118;
t83 = t51 * qJ(3) + t55 * qJD(3);
t82 = t89 * qJD(4);
t76 = t108 * t121;
t64 = qJ(3) * t72;
t62 = qJD(3) * t70;
t50 = -0.2e1 * t99;
t49 = 0.2e1 * t99;
t45 = t51 * t70;
t41 = 0.2e1 * (t65 - t66) * qJD(4);
t35 = t108 * t60;
t17 = -0.2e1 * t117;
t16 = 0.2e1 * t118;
t7 = 0.2e1 * t24 * t42 - 0.2e1 * t43 * t25;
t1 = 0.2e1 * t85;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t60, -0.2e1 * t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t60, 0.2e1 * t51, 0.2e1 * t98 * t60 + 0.2e1 * t116, t50, t41, 0, t49, 0, 0, 0.2e1 * t110, -0.2e1 * t55 * t105 + 0.2e1 * t45, -0.2e1 * t35, 0.2e1 * t116 + 0.2e1 * (-t108 * t120 - t76) * t60, t17, t7, 0, t16, 0, 0, 0.2e1 * t114, 0.2e1 * t113, 0.2e1 * t95, 0.2e1 * t20 * t6 - 0.2e1 * t21 * t5 + 0.2e1 * t48 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t72 + t101, -pkin(2) * t60 + t83, t50, t41, 0, t49, 0, 0, t109 + t110, t45 + t62 + (-qJ(3) - t55) * t105, -t35, -t76 * t60 + t83, t17, t7, 0, t16, 0, 0, t112 + t114, t111 + t113, t94 + t95, -t21 * t12 + t20 * t13 + t6 * t28 - t5 * t29 + t40 * t56 + t48 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t64, t50, t41, 0, t49, 0, 0, 0.2e1 * t109, -0.2e1 * t67 * t103 + 0.2e1 * t62, 0, t64, t17, t7, 0, t16, 0, 0, 0.2e1 * t112, 0.2e1 * t111, 0.2e1 * t94, -0.2e1 * t29 * t12 + 0.2e1 * t28 * t13 + 0.2e1 * t56 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, 0, -t104, 0, -t67 * t82 + t90, -t70 * t82 - t91, 0, 0, 0, 0, -t24, 0, -t25, 0, t6, t5, t4, (-t119 * t5 + t6 * t69 + (-t119 * t20 + t21 * t69) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, 0, -t104, 0, t67 * t96, t70 * t96, 0, 0, 0, 0, -t24, 0, -t25, 0, t13, t12, t4, (-t119 * t12 + t13 * t69 + (-t119 * t28 + t29 * t69) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, -0.2e1 * t100, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, 0, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t100, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
