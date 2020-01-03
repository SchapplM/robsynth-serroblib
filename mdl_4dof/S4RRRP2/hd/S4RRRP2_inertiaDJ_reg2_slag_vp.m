% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:10
% DurationCPUTime: 0.32s
% Computational Cost: add. (185->64), mult. (540->113), div. (0->0), fcn. (313->4), ass. (0->49)
t40 = cos(qJ(2));
t55 = t40 * pkin(1);
t37 = sin(qJ(3));
t31 = t37 * qJD(3);
t30 = pkin(3) * t31;
t38 = sin(qJ(2));
t50 = pkin(1) * qJD(2);
t46 = t38 * t50;
t14 = t30 + t46;
t39 = cos(qJ(3));
t28 = -t39 * pkin(3) - pkin(2);
t19 = t28 - t55;
t33 = t39 * qJD(3);
t54 = t14 * t37 + t19 * t33;
t53 = -qJ(4) - pkin(6);
t27 = -pkin(2) - t55;
t52 = t27 * t33 + t37 * t46;
t35 = t37 ^ 2;
t51 = qJD(3) * t35 * pkin(3) + t28 * t33;
t26 = t38 * pkin(1) + pkin(6);
t49 = qJ(4) + t26;
t48 = pkin(2) * t31;
t47 = pkin(2) * t33;
t45 = t40 * t50;
t44 = pkin(3) * t33;
t9 = t19 * t31;
t17 = t28 * t31;
t43 = t37 * t33;
t36 = t39 ^ 2;
t42 = (t35 + t36) * t40;
t41 = t27 * t31 - t39 * t46;
t34 = t39 * qJ(4);
t32 = t39 * qJD(4);
t25 = -0.2e1 * t43;
t24 = 0.2e1 * t43;
t22 = t39 * t45;
t21 = t39 * pkin(6) + t34;
t20 = t53 * t37;
t13 = 0.2e1 * (-t35 + t36) * qJD(3);
t12 = t39 * t26 + t34;
t11 = t49 * t37;
t7 = -t37 * qJD(4) + t53 * t33;
t6 = -t53 * t31 - t32;
t5 = t42 * t50;
t4 = t6 * t39;
t3 = (-qJD(4) - t45) * t37 - t49 * t33;
t2 = t49 * t31 - t22 - t32;
t1 = t2 * t39;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46, -0.2e1 * t45, 0, 0, t24, t13, 0, t25, 0, 0, 0.2e1 * t41, 0.2e1 * t52, 0.2e1 * t5, 0.2e1 * (t26 * t42 + t27 * t38) * t50, t24, t13, 0, t25, 0, 0, -0.2e1 * t14 * t39 + 0.2e1 * t9, 0.2e1 * t54, -0.2e1 * t3 * t37 - 0.2e1 * t1 + 0.2e1 * (t11 * t39 - t12 * t37) * qJD(3), -0.2e1 * t11 * t3 - 0.2e1 * t12 * t2 + 0.2e1 * t19 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, 0, t24, t13, 0, t25, 0, 0, t41 - t48, -t47 + t52, t5, (-pkin(2) * t38 + pkin(6) * t42) * t50, t24, t13, 0, t25, 0, 0, t17 + t9 + (-t14 - t30) * t39, t51 + t54, -t1 - t4 + (-t3 - t7) * t37 + ((t11 - t20) * t39 + (-t12 - t21) * t37) * qJD(3), pkin(3) * t9 - t11 * t7 - t12 * t6 + t14 * t28 - t2 * t21 + t3 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t13, 0, t25, 0, 0, -0.2e1 * t48, -0.2e1 * t47, 0, 0, t24, t13, 0, t25, 0, 0, -0.2e1 * pkin(3) * t43 + 0.2e1 * t17, 0.2e1 * t51, -0.2e1 * t7 * t37 - 0.2e1 * t4 + 0.2e1 * (-t20 * t39 - t21 * t37) * qJD(3), 0.2e1 * pkin(3) * t17 + 0.2e1 * t20 * t7 - 0.2e1 * t21 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -t26 * t33 - t37 * t45, t26 * t31 - t22, 0, 0, 0, 0, t33, 0, -t31, 0, t3, t2, -t44, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, 0, -pkin(6) * t33, pkin(6) * t31, 0, 0, 0, 0, t33, 0, -t31, 0, t7, t6, -t44, t7 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
