% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:14
% DurationCPUTime: 0.73s
% Computational Cost: add. (1323->96), mult. (2983->164), div. (0->0), fcn. (2908->6), ass. (0->58)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t33 = t83 * t44 + t81 * t45;
t80 = pkin(6) + qJ(2);
t86 = t83 * t80;
t84 = t81 * qJD(2) + qJD(3) * t86;
t87 = t80 * t81;
t89 = -t83 * qJD(2) + qJD(3) * t87;
t18 = t44 * t89 - t84 * t45;
t79 = t33 * qJD(3);
t88 = 0.2e1 * t79;
t17 = t84 * t44 + t89 * t45;
t47 = 2 * qJD(5);
t82 = cos(qJ(4));
t24 = -t44 * t87 + t45 * t86;
t46 = sin(qJ(4));
t78 = qJD(4) * t46;
t60 = t81 * t44 - t83 * t45;
t57 = t60 * qJD(3);
t58 = t46 * t60;
t71 = qJD(4) * t82;
t12 = -qJD(4) * t58 + t33 * t71 - t46 * t57 + t82 * t79;
t55 = t82 * t60;
t21 = t46 * t33 + t55;
t77 = 0.2e1 * t21 * t12;
t76 = pkin(3) * t78;
t20 = -t60 * pkin(7) + t24;
t23 = -t44 * t86 - t45 * t87;
t52 = -t33 * pkin(7) + t23;
t51 = t46 * t52;
t10 = t82 * t20 + t51;
t48 = -pkin(7) * t57 - t18;
t49 = -t79 * pkin(7) - t17;
t50 = t82 * t52;
t3 = -qJD(4) * t50 + t20 * t78 + t46 * t48 - t82 * t49;
t4 = qJD(4) * t51 + t20 * t71 + t46 * t49 + t82 * t48;
t9 = t46 * t20 - t50;
t75 = -t10 * t3 + t9 * t4;
t39 = -t45 * pkin(2) - pkin(1);
t74 = t79 * pkin(3);
t68 = 0.2e1 * (t44 ^ 2 + t45 ^ 2) * qJD(2);
t11 = t33 * t78 + t46 * t79 - (-qJD(3) - qJD(4)) * t55;
t22 = t82 * t33 - t58;
t63 = t11 * t21 - t22 * t12;
t59 = -0.2e1 * t10 * t12 - 0.2e1 * t9 * t11 + 0.2e1 * t3 * t21 + 0.2e1 * t4 * t22;
t56 = -0.2e1 * t57;
t25 = t60 * pkin(3) + t39;
t41 = pkin(3) * t71;
t40 = -t82 * pkin(3) - pkin(4);
t38 = t46 * pkin(3) + qJ(5);
t37 = -0.2e1 * t76;
t34 = t41 + qJD(5);
t8 = t21 * pkin(4) - t22 * qJ(5) + t25;
t7 = -0.2e1 * t22 * t11;
t5 = t12 * pkin(4) + t11 * qJ(5) - t22 * qJD(5) + t74;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, qJ(2) * t68, t33 * t56, 0.2e1 * t60 ^ 2 * qJD(3) - 0.2e1 * t33 * t79, 0, t60 * t88, 0, 0, t39 * t88, t39 * t56, -0.2e1 * t18 * t33 - 0.2e1 * t24 * t79 + 0.2e1 * (qJD(3) * t23 + t17) * t60, -0.2e1 * t24 * t17 + 0.2e1 * t23 * t18, t7, 0.2e1 * t63, 0, t77, 0, 0, 0.2e1 * t25 * t12 + 0.2e1 * t21 * t74, -0.2e1 * t25 * t11 + 0.2e1 * t22 * t74, t59, 0.2e1 * t25 * t74 + 0.2e1 * t75, t7, 0, -0.2e1 * t63, 0, 0, t77, 0.2e1 * t8 * t12 + 0.2e1 * t5 * t21, t59, 0.2e1 * t8 * t11 - 0.2e1 * t5 * t22, 0.2e1 * t8 * t5 + 0.2e1 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t57, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t74, 0, 0, 0, 0, 0, 0, t12, 0, t11, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, -t79, 0, t18, t17, 0, 0, 0, 0, -t11, 0, -t12, 0, -t4, t3, (t82 * t11 - t12 * t46 + (-t82 * t21 + t22 * t46) * qJD(4)) * pkin(3), (-t82 * t4 - t3 * t46 + (t82 * t10 + t46 * t9) * qJD(4)) * pkin(3), 0, -t11, 0, 0, t12, 0, -t4, -t40 * t11 - t38 * t12 - t34 * t21 + t22 * t76, -t3, t10 * t34 - t3 * t38 + t4 * t40 + t9 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -0.2e1 * t41, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0.2e1 * t34, 0.2e1 * t38 * t34 + 0.2e1 * t40 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, -t4, t3, 0, 0, 0, -t11, 0, 0, t12, 0, -t4, pkin(4) * t11 - t12 * qJ(5) - t21 * qJD(5), -t3, -t4 * pkin(4) - t3 * qJ(5) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t41, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, t47 + t41, -pkin(4) * t76 + t34 * qJ(5) + t38 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, qJ(5) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
