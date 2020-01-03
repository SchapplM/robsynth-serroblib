% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP2
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:32
% DurationCPUTime: 0.66s
% Computational Cost: add. (377->78), mult. (896->128), div. (0->0), fcn. (604->6), ass. (0->51)
t45 = sin(qJ(3));
t49 = cos(pkin(8)) * pkin(1) + pkin(2);
t60 = cos(qJ(3));
t62 = sin(pkin(8)) * pkin(1);
t64 = t45 * t62 - t60 * t49;
t57 = t45 * t49 + t60 * t62;
t18 = t57 * qJD(3);
t44 = sin(qJ(4));
t37 = t44 * qJD(4);
t36 = pkin(4) * t37;
t10 = t36 + t18;
t19 = -pkin(3) + t64;
t46 = cos(qJ(4));
t61 = t46 * pkin(4);
t15 = t19 - t61;
t39 = t46 * qJD(4);
t63 = t10 * t44 + t15 * t39;
t59 = -qJ(5) - pkin(7);
t58 = t18 * t44 + t19 * t39;
t34 = -pkin(3) - t61;
t41 = t44 ^ 2;
t56 = qJD(4) * t41 * pkin(4) + t34 * t39;
t20 = pkin(7) + t57;
t55 = qJ(5) + t20;
t54 = pkin(3) * t37;
t53 = pkin(3) * t39;
t52 = pkin(4) * t39;
t6 = t15 * t37;
t25 = t34 * t37;
t51 = t44 * t39;
t17 = t64 * qJD(3);
t42 = t46 ^ 2;
t4 = (-t41 - t42) * t17;
t50 = -t18 * t46 + t19 * t37;
t40 = t46 * qJ(5);
t38 = t46 * qJD(5);
t31 = -0.2e1 * t51;
t30 = 0.2e1 * t51;
t28 = pkin(7) * t46 + t40;
t27 = t59 * t44;
t23 = 0.2e1 * (-t41 + t42) * qJD(4);
t22 = -t44 * qJD(5) + t39 * t59;
t21 = -t37 * t59 - t38;
t16 = t21 * t46;
t11 = t46 * t17;
t9 = t20 * t46 + t40;
t8 = t55 * t44;
t3 = (-qJD(5) + t17) * t44 - t55 * t39;
t2 = t37 * t55 + t11 - t38;
t1 = t2 * t46;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t18, 0.2e1 * t17, 0, -0.2e1 * t17 * t57 + 0.2e1 * t18 * t64, t30, t23, 0, t31, 0, 0, 0.2e1 * t50, 0.2e1 * t58, 0.2e1 * t4, 0.2e1 * t19 * t18 + 0.2e1 * t20 * t4, t30, t23, 0, t31, 0, 0, -0.2e1 * t10 * t46 + 0.2e1 * t6, 0.2e1 * t63, -0.2e1 * t3 * t44 - 0.2e1 * t1 + 0.2e1 * (-t44 * t9 + t46 * t8) * qJD(4), 0.2e1 * t10 * t15 - 0.2e1 * t2 * t9 - 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t44 + t3 * t46 + (t44 * t8 + t46 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, 0, t30, t23, 0, t31, 0, 0, t50 - t54, -t53 + t58, t4, -t18 * pkin(3) + pkin(7) * t4, t30, t23, 0, t31, 0, 0, t25 + t6 + (-t10 - t36) * t46, t56 + t63, -t1 - t16 + (-t22 - t3) * t44 + ((-t27 + t8) * t46 + (-t28 - t9) * t44) * qJD(4), pkin(4) * t6 + t10 * t34 - t2 * t28 - t21 * t9 - t22 * t8 + t27 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t21 + t46 * t22 + (-t27 * t44 + t28 * t46) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t23, 0, t31, 0, 0, -0.2e1 * t54, -0.2e1 * t53, 0, 0, t30, t23, 0, t31, 0, 0, -0.2e1 * pkin(4) * t51 + 0.2e1 * t25, 0.2e1 * t56, -0.2e1 * t22 * t44 - 0.2e1 * t16 + 0.2e1 * (-t27 * t46 - t28 * t44) * qJD(4), 0.2e1 * pkin(4) * t25 - 0.2e1 * t21 * t28 + 0.2e1 * t22 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, 0, t44 * t17 - t20 * t39, t20 * t37 + t11, 0, 0, 0, 0, t39, 0, -t37, 0, t3, t2, -t52, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, 0, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, 0, -pkin(7) * t39, pkin(7) * t37, 0, 0, 0, 0, t39, 0, -t37, 0, t22, t21, -t52, t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
