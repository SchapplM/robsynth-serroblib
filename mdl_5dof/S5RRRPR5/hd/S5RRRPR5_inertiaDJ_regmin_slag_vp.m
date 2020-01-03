% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR5
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
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:17
% EndTime: 2019-12-31 21:14:19
% DurationCPUTime: 0.71s
% Computational Cost: add. (1549->119), mult. (3597->235), div. (0->0), fcn. (3378->8), ass. (0->95)
t62 = cos(qJ(5));
t57 = t62 ^ 2;
t59 = sin(qJ(5));
t97 = t59 ^ 2 - t57;
t81 = t97 * qJD(5);
t106 = qJD(2) + qJD(3);
t61 = sin(qJ(2));
t105 = pkin(6) + pkin(7);
t85 = qJD(2) * t105;
t43 = t61 * t85;
t64 = cos(qJ(2));
t44 = t64 * t85;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t46 = t105 * t61;
t47 = t105 * t64;
t73 = -t63 * t46 - t60 * t47;
t20 = -t73 * qJD(3) + t63 * t43 + t60 * t44;
t42 = t60 * t64 + t63 * t61;
t30 = t106 * t42;
t71 = t60 * t61 - t63 * t64;
t12 = -t30 * qJ(4) - t71 * qJD(4) - t20;
t58 = sin(pkin(9));
t72 = t60 * t46 - t63 * t47;
t21 = t72 * qJD(3) + t60 * t43 - t63 * t44;
t29 = t106 * t71;
t65 = t29 * qJ(4) - t42 * qJD(4) + t21;
t95 = cos(pkin(9));
t4 = t58 * t12 - t95 * t65;
t3 = t4 * t59;
t23 = -t71 * qJ(4) - t72;
t68 = -t42 * qJ(4) + t73;
t16 = t58 * t23 - t95 * t68;
t55 = qJD(5) * t62;
t104 = t16 * t55 + t3;
t28 = t95 * t42 - t58 * t71;
t103 = t28 * t62;
t102 = t58 * t60;
t18 = -t58 * t29 + t95 * t30;
t101 = t59 * t18;
t100 = t62 * t18;
t19 = -t95 * t29 - t58 * t30;
t99 = t62 * t19;
t52 = t63 * pkin(2) + pkin(3);
t38 = -pkin(2) * t102 + t95 * t52;
t34 = -pkin(4) - t38;
t82 = t95 * t60;
t96 = pkin(2) * qJD(3);
t36 = (t58 * t63 + t82) * t96;
t98 = t34 * t55 + t36 * t59;
t39 = pkin(2) * t82 + t58 * t52;
t94 = qJD(5) * t59;
t93 = t61 * qJD(2);
t92 = t64 * qJD(2);
t91 = -0.2e1 * pkin(1) * qJD(2);
t54 = pkin(2) * t93;
t90 = t60 * t96;
t89 = t63 * t96;
t51 = -t95 * pkin(3) - pkin(4);
t88 = t51 * t94;
t87 = t51 * t55;
t86 = t59 * t55;
t53 = -t64 * pkin(2) - pkin(1);
t24 = t30 * pkin(3) + t54;
t84 = -0.4e1 * t59 * t103;
t83 = t34 * t94 - t36 * t62;
t27 = t58 * t42 + t95 * t71;
t66 = t71 * pkin(3) + t53;
t15 = t27 * pkin(4) - t28 * pkin(8) + t66;
t17 = t95 * t23 + t58 * t68;
t80 = t62 * t15 - t59 * t17;
t79 = t59 * t15 + t62 * t17;
t50 = t58 * pkin(3) + pkin(8);
t78 = -t18 * t50 + t19 * t51;
t35 = pkin(8) + t39;
t77 = t27 * t35 - t28 * t34;
t37 = (t95 * t63 - t102) * t96;
t76 = -t37 * t27 + t36 * t28;
t75 = t27 * t50 - t28 * t51;
t10 = t27 * t55 + t101;
t70 = t59 * t19 + t28 * t55;
t69 = -t28 * t94 + t99;
t67 = -t18 * t35 + t19 * t34 + t76;
t48 = 0.2e1 * t86;
t41 = -0.2e1 * t81;
t25 = t28 ^ 2;
t13 = t16 * t94;
t9 = -t27 * t94 + t100;
t8 = -t28 * t81 + t59 * t99;
t7 = t18 * pkin(4) - t19 * pkin(8) + t24;
t6 = qJD(5) * t84 - t97 * t19;
t5 = t95 * t12 + t58 * t65;
t2 = -t79 * qJD(5) - t59 * t5 + t62 * t7;
t1 = -t80 * qJD(5) - t62 * t5 - t59 * t7;
t11 = [0, 0, 0, 0.2e1 * t61 * t92, 0.2e1 * (-t61 ^ 2 + t64 ^ 2) * qJD(2), 0, 0, 0, t61 * t91, t64 * t91, -0.2e1 * t42 * t29, 0.2e1 * t29 * t71 - 0.2e1 * t42 * t30, 0, 0, 0, 0.2e1 * t53 * t30 + 0.2e1 * t71 * t54, -0.2e1 * t53 * t29 + 0.2e1 * t42 * t54, 0.2e1 * t16 * t19 - 0.2e1 * t17 * t18 - 0.2e1 * t5 * t27 + 0.2e1 * t4 * t28, 0.2e1 * t16 * t4 + 0.2e1 * t17 * t5 + 0.2e1 * t66 * t24, 0.2e1 * t57 * t28 * t19 - 0.2e1 * t25 * t86, t19 * t84 + 0.2e1 * t25 * t81, 0.2e1 * t28 * t100 + 0.2e1 * t69 * t27, -0.2e1 * t28 * t101 - 0.2e1 * t70 * t27, 0.2e1 * t27 * t18, 0.2e1 * t70 * t16 + 0.2e1 * t80 * t18 + 0.2e1 * t2 * t27 + 0.2e1 * t28 * t3, 0.2e1 * t1 * t27 + 0.2e1 * t4 * t103 + 0.2e1 * t69 * t16 - 0.2e1 * t79 * t18; 0, 0, 0, 0, 0, t92, -t93, 0, -pkin(6) * t92, pkin(6) * t93, 0, 0, -t29, -t30, 0, t21, t20, -t39 * t18 - t38 * t19 + t76, t16 * t36 + t17 * t37 - t4 * t38 + t5 * t39, t8, t6, t10, t9, 0, t13 + (-t77 * qJD(5) - t4) * t62 + t67 * t59, t67 * t62 + t77 * t94 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t90, -0.2e1 * t89, 0, -0.2e1 * t38 * t36 + 0.2e1 * t39 * t37, t48, t41, 0, 0, 0, 0.2e1 * t83, 0.2e1 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, t21, t20, (-t18 * t58 - t95 * t19) * pkin(3), (-t95 * t4 + t5 * t58) * pkin(3), t8, t6, t10, t9, 0, t13 + t78 * t59 + (-t75 * qJD(5) - t4) * t62, t78 * t62 + t75 * t94 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89, 0, (-t95 * t36 + t37 * t58) * pkin(3), t48, t41, 0, 0, 0, t83 + t88, t87 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t41, 0, 0, 0, 0.2e1 * t88, 0.2e1 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t70, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t94, 0, -t35 * t55 - t59 * t37, t35 * t94 - t62 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t94, 0, -t50 * t55, t50 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
