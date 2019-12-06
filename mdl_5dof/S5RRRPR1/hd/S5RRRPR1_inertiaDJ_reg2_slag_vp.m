% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:38:56
% DurationCPUTime: 1.36s
% Computational Cost: add. (3337->155), mult. (7486->292), div. (0->0), fcn. (7072->8), ass. (0->91)
t80 = sin(qJ(3));
t81 = sin(qJ(2));
t82 = cos(qJ(3));
t83 = cos(qJ(2));
t61 = t80 * t83 + t82 * t81;
t118 = -pkin(7) - pkin(6);
t67 = t118 * t81;
t62 = t80 * t67;
t68 = t118 * t83;
t42 = -t82 * t68 + t62;
t109 = qJD(3) * t82;
t110 = qJD(3) * t80;
t106 = t83 * qJD(2);
t114 = t80 * t81;
t119 = qJD(2) + qJD(3);
t39 = -t82 * t106 - t83 * t109 + t119 * t114;
t111 = t82 * t83;
t88 = (t118 * t111 - t62) * qJD(2);
t120 = t39 * qJ(4) - t61 * qJD(4) + t68 * t109 - t67 * t110 + t88;
t77 = sin(pkin(9));
t117 = t77 * pkin(3);
t116 = cos(qJ(5));
t78 = cos(pkin(9));
t115 = t78 * t80;
t41 = t82 * t67 + t80 * t68;
t32 = -t61 * qJ(4) + t41;
t60 = -t111 + t114;
t33 = -t60 * qJ(4) + t42;
t15 = t77 * t32 + t78 * t33;
t79 = sin(qJ(5));
t108 = qJD(5) * t79;
t107 = t81 * qJD(2);
t105 = t79 * t117;
t104 = -0.2e1 * pkin(1) * qJD(2);
t74 = t82 * pkin(2) + pkin(3);
t53 = -t77 * t80 * pkin(2) + t78 * t74;
t50 = pkin(4) + t53;
t51 = (-t77 * t82 - t115) * qJD(3) * pkin(2);
t101 = pkin(2) * t109;
t102 = pkin(2) * t110;
t52 = t78 * t101 - t77 * t102;
t95 = qJD(5) * t116;
t103 = -t116 * t52 - t50 * t95 - t79 * t51;
t76 = pkin(2) * t107;
t100 = t81 * t106;
t75 = -t83 * pkin(2) - pkin(1);
t55 = pkin(2) * t115 + t77 * t74;
t99 = t116 * t55;
t40 = t119 * t61;
t34 = t40 * pkin(3) + t76;
t14 = t78 * t32 - t77 * t33;
t97 = -t77 * t39 + t78 * t40;
t96 = t116 * t51 - t79 * t52;
t94 = t116 * t117;
t93 = t78 * t60 + t77 * t61;
t49 = t60 * pkin(3) + t75;
t38 = -t77 * t60 + t78 * t61;
t92 = -t38 * pkin(8) + t14;
t31 = t79 * t50 + t99;
t24 = -t61 * qJD(2) * t118 - t67 * t109 - t68 * t110;
t72 = t78 * pkin(3) + pkin(4);
t56 = t79 * t72 + t94;
t91 = t79 * t92;
t90 = t116 * t93;
t89 = t116 * t92;
t19 = t116 * t38 - t79 * t93;
t87 = -t40 * qJ(4) - t60 * qJD(4) - t24;
t9 = t120 * t78 - t77 * t87;
t10 = t120 * t77 + t78 * t87;
t23 = -t78 * t39 - t77 * t40;
t85 = -t23 * pkin(8) + t9;
t84 = -pkin(8) * t97 + t10;
t65 = t72 * t95;
t54 = t116 * t72 - t105;
t48 = t56 * qJD(5);
t47 = qJD(5) * t105 - t65;
t30 = t116 * t50 - t79 * t55;
t26 = pkin(4) * t93 + t49;
t25 = -t42 * qJD(3) + t88;
t18 = t79 * t38 + t90;
t17 = -qJD(5) * t31 + t96;
t16 = t55 * t108 + t103;
t13 = pkin(4) * t97 + t34;
t12 = -pkin(8) * t93 + t15;
t8 = t19 * qJD(5) + t116 * t97 + t79 * t23;
t7 = qJD(5) * t90 + t38 * t108 - t116 * t23 + t79 * t97;
t6 = t116 * t12 + t91;
t5 = -t79 * t12 + t89;
t2 = -qJD(5) * t91 + t116 * t85 - t12 * t95 - t79 * t84;
t1 = -qJD(5) * t89 + t12 * t108 - t116 * t84 - t79 * t85;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100, 0.2e1 * (-t81 ^ 2 + t83 ^ 2) * qJD(2), 0, -0.2e1 * t100, 0, 0, t81 * t104, t83 * t104, 0, 0, -0.2e1 * t61 * t39, 0.2e1 * t39 * t60 - 0.2e1 * t61 * t40, 0, 0.2e1 * t60 * t40, 0, 0, 0.2e1 * t75 * t40 + 0.2e1 * t60 * t76, -0.2e1 * t75 * t39 + 0.2e1 * t61 * t76, 0.2e1 * t24 * t60 - 0.2e1 * t25 * t61 + 0.2e1 * t41 * t39 - 0.2e1 * t42 * t40, -0.2e1 * t42 * t24 + 0.2e1 * t41 * t25 + 0.2e1 * t75 * t76, 0.2e1 * t38 * t23, -0.2e1 * t23 * t93 - 0.2e1 * t38 * t97, 0, 0.2e1 * t93 * t97, 0, 0, 0.2e1 * t34 * t93 + 0.2e1 * t49 * t97, 0.2e1 * t49 * t23 + 0.2e1 * t34 * t38, -0.2e1 * t10 * t93 - 0.2e1 * t14 * t23 - 0.2e1 * t15 * t97 - 0.2e1 * t9 * t38, 0.2e1 * t15 * t10 + 0.2e1 * t14 * t9 + 0.2e1 * t49 * t34, -0.2e1 * t19 * t7, 0.2e1 * t7 * t18 - 0.2e1 * t19 * t8, 0, 0.2e1 * t18 * t8, 0, 0, 0.2e1 * t13 * t18 + 0.2e1 * t26 * t8, 0.2e1 * t13 * t19 - 0.2e1 * t26 * t7, 0.2e1 * t1 * t18 - 0.2e1 * t2 * t19 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t26 * t13 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, -t107, 0, -pkin(6) * t106, pkin(6) * t107, 0, 0, 0, 0, -t39, 0, -t40, 0, t25, t24, (t39 * t82 - t40 * t80 + (-t60 * t82 + t61 * t80) * qJD(3)) * pkin(2), (-t24 * t80 + t25 * t82 + (-t41 * t80 + t42 * t82) * qJD(3)) * pkin(2), 0, 0, t23, 0, -t97, 0, t9, -t10, -t53 * t23 - t51 * t38 - t52 * t93 - t55 * t97, t10 * t55 + t14 * t51 + t15 * t52 + t9 * t53, 0, 0, -t7, 0, -t8, 0, t2, t1, t16 * t18 - t17 * t19 + t30 * t7 - t31 * t8, -t1 * t31 - t6 * t16 + t5 * t17 + t2 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t102, -0.2e1 * t101, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t51, -0.2e1 * t52, 0, 0.2e1 * t53 * t51 + 0.2e1 * t55 * t52, 0, 0, 0, 0, 0, 0, 0.2e1 * t17, 0.2e1 * t16, 0, -0.2e1 * t31 * t16 + 0.2e1 * t30 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t40, 0, t25, t24, 0, 0, 0, 0, t23, 0, -t97, 0, t9, -t10, (-t78 * t23 - t77 * t97) * pkin(3), (t10 * t77 + t78 * t9) * pkin(3), 0, 0, -t7, 0, -t8, 0, t2, t1, t47 * t18 + t48 * t19 + t54 * t7 - t56 * t8, -t1 * t56 + t2 * t54 - t6 * t47 - t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t101, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t52, 0, (t51 * t78 + t52 * t77) * pkin(3), 0, 0, 0, 0, 0, 0, (-t94 - t99 + (-t50 - t72) * t79) * qJD(5) + t96, -t65 + (t55 + t117) * t108 + t103, 0, -t16 * t56 + t17 * t54 - t30 * t48 - t31 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t48, 0.2e1 * t47, 0, -0.2e1 * t56 * t47 - 0.2e1 * t54 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t23, 0, t34, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
