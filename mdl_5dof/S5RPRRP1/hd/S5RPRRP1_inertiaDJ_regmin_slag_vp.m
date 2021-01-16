% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:19
% EndTime: 2021-01-15 12:27:22
% DurationCPUTime: 0.34s
% Computational Cost: add. (533->73), mult. (1078->114), div. (0->0), fcn. (914->4), ass. (0->45)
t31 = sin(qJ(4));
t32 = sin(qJ(3));
t46 = t32 * qJD(3);
t48 = qJD(4) * t31;
t33 = cos(qJ(4));
t34 = cos(qJ(3));
t49 = t33 * t34;
t58 = qJD(3) + qJD(4);
t10 = -t31 * t46 - t32 * t48 + t58 * t49;
t17 = t31 * t34 + t33 * t32;
t18 = -t31 * t32 + t49;
t52 = t33 * pkin(3);
t28 = pkin(4) + t52;
t9 = t58 * t17;
t59 = -t28 * t9 + (t10 * t31 + (t17 * t33 - t18 * t31) * qJD(4)) * pkin(3);
t56 = 2 * qJD(2);
t55 = t9 * pkin(4);
t54 = t18 * t9;
t35 = -pkin(1) - pkin(6);
t51 = pkin(7) - t35;
t26 = t32 * pkin(3) + qJ(2);
t47 = qJD(4) * t33;
t45 = t34 * qJD(3);
t21 = pkin(3) * t45 + qJD(2);
t44 = qJ(2) * qJD(3);
t43 = pkin(3) * t48;
t42 = pkin(3) * t47;
t20 = t51 * t34;
t41 = qJD(3) * t20;
t40 = t51 * t46;
t39 = t17 * t10 - t54;
t19 = t51 * t32;
t37 = t33 * t19 + t31 * t20;
t3 = -t19 * t48 + t20 * t47 - t31 * t40 + t33 * t41;
t1 = t10 * qJ(5) + t17 * qJD(5) + t3;
t4 = t37 * qJD(4) + t31 * t41 + t33 * t40;
t2 = t9 * qJ(5) - t18 * qJD(5) + t4;
t5 = -t18 * qJ(5) + t31 * t19 - t33 * t20;
t6 = -t17 * qJ(5) - t37;
t36 = t1 * t17 - t6 * t10 - t2 * t18 + t5 * t9;
t25 = -0.2e1 * t42;
t24 = -0.2e1 * t43;
t11 = t17 * pkin(4) + t26;
t7 = t10 * pkin(4) + t21;
t8 = [0, 0, 0, 0, t56, qJ(2) * t56, -0.2e1 * t32 * t45, 0.2e1 * (t32 ^ 2 - t34 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t32 + 0.2e1 * t34 * t44, 0.2e1 * qJD(2) * t34 - 0.2e1 * t32 * t44, -0.2e1 * t54, -0.2e1 * t10 * t18 + 0.2e1 * t17 * t9, 0, 0, 0, 0.2e1 * t26 * t10 + 0.2e1 * t21 * t17, 0.2e1 * t21 * t18 - 0.2e1 * t26 * t9, 0.2e1 * t11 * t10 + 0.2e1 * t7 * t17, -0.2e1 * t11 * t9 + 0.2e1 * t7 * t18, 0.2e1 * t36, -0.2e1 * t6 * t1 + 0.2e1 * t11 * t7 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t39, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, -t35 * t46, -t35 * t45, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, -t59, t2 * t28 + (-t1 * t31 + (-t31 * t5 + t33 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, 0, 0, 0, 0, -t9, -t10, -t9, -t10, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, t24, t25, 0, 0.2e1 * (-t28 + t52) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, t55, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -t9, -t10, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, -t43, -t42, 0, -pkin(4) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
