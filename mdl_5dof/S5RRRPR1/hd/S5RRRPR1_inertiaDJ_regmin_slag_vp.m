% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:57
% EndTime: 2021-01-15 22:49:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (1681->105), mult. (3891->202), div. (0->0), fcn. (3601->8), ass. (0->76)
t88 = qJD(2) + qJD(3);
t66 = sin(qJ(2));
t87 = pkin(6) + pkin(7);
t76 = qJD(2) * t87;
t51 = t66 * t76;
t69 = cos(qJ(2));
t52 = t69 * t76;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t54 = t87 * t66;
t55 = t87 * t69;
t72 = -t68 * t54 - t65 * t55;
t24 = -qJD(3) * t72 + t68 * t51 + t65 * t52;
t62 = sin(pkin(9));
t86 = pkin(3) * t62;
t63 = cos(pkin(9));
t85 = t63 * t65;
t50 = t65 * t69 + t68 * t66;
t30 = -t50 * qJ(4) + t72;
t70 = t65 * t66 - t68 * t69;
t71 = t65 * t54 - t68 * t55;
t31 = -qJ(4) * t70 - t71;
t16 = t62 * t30 + t63 * t31;
t59 = t68 * pkin(2) + pkin(3);
t46 = -t62 * t65 * pkin(2) + t63 * t59;
t43 = pkin(4) + t46;
t58 = t63 * pkin(3) + pkin(4);
t84 = -t43 - t58;
t83 = pkin(2) * qJD(3);
t64 = sin(qJ(5));
t82 = qJD(5) * t64;
t81 = t66 * qJD(2);
t80 = t69 * qJD(2);
t79 = -0.2e1 * pkin(1) * qJD(2);
t61 = pkin(2) * t81;
t78 = t65 * t83;
t77 = t68 * t83;
t60 = -t69 * pkin(2) - pkin(1);
t37 = t88 * t50;
t32 = t37 * pkin(3) + t61;
t13 = -t37 * qJ(4) - qJD(4) * t70 - t24;
t25 = qJD(3) * t71 + t65 * t51 - t68 * t52;
t36 = t88 * t70;
t14 = t36 * qJ(4) - t50 * qJD(4) + t25;
t7 = -t62 * t13 + t63 * t14;
t15 = t63 * t30 - t62 * t31;
t44 = (-t62 * t68 - t85) * t83;
t47 = pkin(2) * t85 + t62 * t59;
t75 = -t64 * t44 + t47 * t82;
t45 = -t62 * t78 + t63 * t77;
t67 = cos(qJ(5));
t74 = t67 * t44 - t64 * t45;
t8 = t63 * t13 + t62 * t14;
t34 = t62 * t50 + t63 * t70;
t35 = t63 * t50 - t62 * t70;
t19 = t67 * t34 + t64 * t35;
t20 = -t64 * t34 + t67 * t35;
t42 = pkin(3) * t70 + t60;
t57 = t82 * t86;
t41 = (-t58 * t64 - t67 * t86) * qJD(5);
t40 = -qJD(5) * t67 * t58 + t57;
t26 = t34 * pkin(4) + t42;
t23 = -t63 * t36 - t62 * t37;
t22 = -t62 * t36 + t63 * t37;
t18 = (-t43 * t64 - t47 * t67) * qJD(5) + t74;
t17 = (-qJD(5) * t43 - t45) * t67 + t75;
t12 = t22 * pkin(4) + t32;
t10 = -t34 * pkin(8) + t16;
t9 = -t35 * pkin(8) + t15;
t6 = qJD(5) * t20 + t67 * t22 + t64 * t23;
t5 = -qJD(5) * t19 - t64 * t22 + t67 * t23;
t4 = -t22 * pkin(8) + t8;
t3 = -t23 * pkin(8) + t7;
t2 = t67 * t3 - t64 * t4 + (-t10 * t67 - t64 * t9) * qJD(5);
t1 = -t64 * t3 - t67 * t4 + (t10 * t64 - t67 * t9) * qJD(5);
t11 = [0, 0, 0, 0.2e1 * t66 * t80, 0.2e1 * (-t66 ^ 2 + t69 ^ 2) * qJD(2), 0, 0, 0, t66 * t79, t69 * t79, -0.2e1 * t50 * t36, 0.2e1 * t36 * t70 - 0.2e1 * t50 * t37, 0, 0, 0, 0.2e1 * t60 * t37 + 0.2e1 * t61 * t70, -0.2e1 * t60 * t36 + 0.2e1 * t50 * t61, 0.2e1 * t42 * t22 + 0.2e1 * t32 * t34, 0.2e1 * t42 * t23 + 0.2e1 * t32 * t35, -0.2e1 * t15 * t23 - 0.2e1 * t16 * t22 - 0.2e1 * t8 * t34 - 0.2e1 * t7 * t35, 0.2e1 * t15 * t7 + 0.2e1 * t16 * t8 + 0.2e1 * t42 * t32, 0.2e1 * t20 * t5, -0.2e1 * t5 * t19 - 0.2e1 * t20 * t6, 0, 0, 0, 0.2e1 * t12 * t19 + 0.2e1 * t26 * t6, 0.2e1 * t12 * t20 + 0.2e1 * t26 * t5; 0, 0, 0, 0, 0, t80, -t81, 0, -pkin(6) * t80, pkin(6) * t81, 0, 0, -t36, -t37, 0, t25, t24, t7, -t8, -t47 * t22 - t46 * t23 - t45 * t34 - t44 * t35, t15 * t44 + t16 * t45 + t7 * t46 + t8 * t47, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t78, -0.2e1 * t77, 0.2e1 * t44, -0.2e1 * t45, 0, 0.2e1 * t46 * t44 + 0.2e1 * t47 * t45, 0, 0, 0, 0, 0, 0.2e1 * t18, 0.2e1 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t25, t24, t7, -t8, (-t22 * t62 - t23 * t63) * pkin(3), (t62 * t8 + t63 * t7) * pkin(3), 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, t44, -t45, 0, (t44 * t63 + t45 * t62) * pkin(3), 0, 0, 0, 0, 0, ((-t47 - t86) * t67 + t84 * t64) * qJD(5) + t74, t57 + (qJD(5) * t84 - t45) * t67 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t41, 0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t32, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
