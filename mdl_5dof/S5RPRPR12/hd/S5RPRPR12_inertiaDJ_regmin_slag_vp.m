% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:17
% DurationCPUTime: 0.66s
% Computational Cost: add. (896->100), mult. (2211->202), div. (0->0), fcn. (2156->8), ass. (0->69)
t57 = sin(pkin(9));
t59 = cos(pkin(9));
t61 = sin(qJ(5));
t63 = cos(qJ(5));
t91 = -t61 * t57 + t63 * t59;
t90 = t91 * qJD(5);
t41 = t63 * t57 + t61 * t59;
t34 = t41 * qJD(5);
t89 = 0.2e1 * t90;
t60 = cos(pkin(8));
t52 = -t60 * pkin(2) - pkin(1);
t88 = 0.2e1 * t52;
t58 = sin(pkin(8));
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t42 = t64 * t58 + t62 * t60;
t87 = t42 * t57;
t81 = t64 * t60;
t40 = t62 * t58 - t81;
t35 = t40 * qJD(3);
t86 = t57 * t35;
t85 = t59 * t35;
t80 = pkin(6) + qJ(2);
t44 = t80 * t58;
t82 = t64 * t44;
t79 = pkin(7) + qJ(4);
t36 = t42 * qJD(3);
t17 = t36 * pkin(3) + t35 * qJ(4) - t42 * qJD(4);
t46 = t80 * t60;
t21 = qJD(3) * t82 - qJD(2) * t81 + (qJD(2) * t58 + qJD(3) * t46) * t62;
t6 = t57 * t17 - t59 * t21;
t28 = t40 * pkin(3) - t42 * qJ(4) + t52;
t31 = -t62 * t44 + t64 * t46;
t12 = t57 * t28 + t59 * t31;
t78 = t57 ^ 2 + t59 ^ 2;
t5 = t59 * t17 + t57 * t21;
t11 = t59 * t28 - t57 * t31;
t30 = t62 * t46 + t82;
t75 = 0.2e1 * t78 * qJD(4);
t74 = 0.2e1 * (t58 ^ 2 + t60 ^ 2) * qJD(2);
t73 = t5 * t59 + t6 * t57;
t72 = -t5 * t57 + t6 * t59;
t7 = -t59 * t42 * pkin(7) + t40 * pkin(4) + t11;
t8 = -pkin(7) * t87 + t12;
t71 = t61 * t8 - t63 * t7;
t70 = t61 * t7 + t63 * t8;
t22 = t42 * qJD(2) + t31 * qJD(3);
t69 = t22 * t42 - t30 * t35;
t68 = -t41 * t36 - t40 * t90;
t43 = t79 * t57;
t45 = t79 * t59;
t67 = -t63 * t43 - t61 * t45;
t66 = -t61 * t43 + t63 * t45;
t65 = pkin(3) * t35 - qJ(4) * t36 - qJD(4) * t40;
t51 = -t59 * pkin(4) - pkin(3);
t25 = t91 * t42;
t24 = t41 * t42;
t23 = pkin(4) * t87 + t30;
t20 = -t41 * qJD(4) - t66 * qJD(5);
t19 = -qJD(4) * t91 - t67 * qJD(5);
t14 = -t34 * t40 + t36 * t91;
t13 = -pkin(4) * t86 + t22;
t10 = -t41 * t35 + t42 * t90;
t9 = -t34 * t42 - t35 * t91;
t4 = pkin(7) * t86 + t6;
t3 = t36 * pkin(4) + pkin(7) * t85 + t5;
t2 = -t70 * qJD(5) + t63 * t3 - t61 * t4;
t1 = t71 * qJD(5) - t61 * t3 - t63 * t4;
t15 = [0, 0, 0, 0, 0, t74, qJ(2) * t74, -0.2e1 * t42 * t35, 0.2e1 * t35 * t40 - 0.2e1 * t42 * t36, 0, 0, 0, t36 * t88, -t35 * t88, 0.2e1 * t11 * t36 + 0.2e1 * t5 * t40 + 0.2e1 * t69 * t57, -0.2e1 * t12 * t36 - 0.2e1 * t6 * t40 + 0.2e1 * t69 * t59, -0.2e1 * t73 * t42 - 0.2e1 * (-t11 * t59 - t12 * t57) * t35, 0.2e1 * t11 * t5 + 0.2e1 * t12 * t6 + 0.2e1 * t30 * t22, 0.2e1 * t25 * t9, -0.2e1 * t25 * t10 - 0.2e1 * t9 * t24, 0.2e1 * t25 * t36 + 0.2e1 * t9 * t40, -0.2e1 * t10 * t40 - 0.2e1 * t24 * t36, 0.2e1 * t40 * t36, 0.2e1 * t23 * t10 + 0.2e1 * t13 * t24 + 0.2e1 * t2 * t40 - 0.2e1 * t71 * t36, 0.2e1 * t1 * t40 + 0.2e1 * t13 * t25 + 0.2e1 * t23 * t9 - 0.2e1 * t70 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, t59 * t36, -t57 * t36, t78 * t35, t73, 0, 0, 0, 0, 0, t14, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t22, t21, -t22 * t59 + t65 * t57, t22 * t57 + t65 * t59, t72, -t22 * pkin(3) + (-t11 * t57 + t12 * t59) * qJD(4) + t72 * qJ(4), t25 * t90 + t9 * t41, -t41 * t10 - t24 * t90 - t25 * t34 + t9 * t91, -t68, t14, 0, t51 * t10 - t13 * t91 + t20 * t40 + t23 * t34 + t67 * t36, t13 * t41 + t19 * t40 + t23 * t90 - t66 * t36 + t51 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, qJ(4) * t75, t41 * t89, -0.2e1 * t41 * t34 + 0.2e1 * t90 * t91, 0, 0, 0, 0.2e1 * t51 * t34, t51 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, t22, 0, 0, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, t36, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t34, 0, t20, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
