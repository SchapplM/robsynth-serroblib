% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP8
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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:31
% DurationCPUTime: 0.34s
% Computational Cost: add. (306->70), mult. (615->113), div. (0->0), fcn. (393->4), ass. (0->50)
t28 = sin(qJ(4));
t26 = t28 ^ 2;
t30 = cos(qJ(4));
t27 = t30 ^ 2;
t59 = t26 + t27;
t58 = -t28 * pkin(4) + t30 * qJ(5);
t49 = cos(qJ(3));
t39 = qJD(3) * t49;
t9 = -qJD(4) * t58 - t28 * qJD(5);
t35 = t30 * pkin(4) + t28 * qJ(5);
t8 = t35 * qJD(4) - t30 * qJD(5);
t57 = 2 * qJD(2);
t56 = 0.2e1 * qJD(5);
t55 = -pkin(1) - pkin(2);
t29 = sin(qJ(3));
t33 = t49 * qJ(2) + t29 * t55;
t7 = t29 * qJD(2) + t33 * qJD(3);
t4 = t7 - t9;
t54 = -t4 + t9;
t52 = t7 * t28;
t51 = t7 * t30;
t18 = -pkin(3) - t35;
t36 = t49 * t55;
t16 = t29 * qJ(2) + pkin(3) - t36;
t5 = t16 + t35;
t50 = t18 - t5;
t47 = qJD(3) * t29;
t46 = t28 * qJD(4);
t23 = t30 * qJD(4);
t43 = -0.2e1 * pkin(3) * qJD(4);
t42 = pkin(7) * t46;
t41 = pkin(7) * t23;
t40 = t28 * t23;
t6 = qJ(2) * t47 - t49 * qJD(2) - qJD(3) * t36;
t1 = t59 * t6;
t38 = qJD(4) * (pkin(3) + t16);
t37 = qJD(4) * t49;
t32 = t59 * t49;
t20 = 0.2e1 * t40;
t19 = (-t26 + t27) * qJD(4);
t17 = -pkin(7) + t33;
t15 = 0.2e1 * t19;
t14 = t28 * t37 + t30 * t47;
t13 = t28 * t47 - t30 * t37;
t12 = t29 * t23 + t28 * t39;
t11 = t32 * qJD(3);
t10 = t29 * t46 - t30 * t39;
t3 = t17 * t23 - t28 * t6;
t2 = t17 * t46 + t30 * t6;
t21 = [0, 0, 0, 0, t57, qJ(2) * t57, 0, 0.2e1 * t7, -0.2e1 * t6, t20, t15, 0, 0, 0, -0.2e1 * t16 * t46 + 0.2e1 * t51, -0.2e1 * t16 * t23 - 0.2e1 * t52, 0.2e1 * t4 * t30 - 0.2e1 * t5 * t46, 0.2e1 * t1, 0.2e1 * t5 * t23 + 0.2e1 * t4 * t28, -0.2e1 * t17 * t1 + 0.2e1 * t5 * t4; 0, 0, 0, 0, 0, 0, 0, t47, t39, 0, 0, 0, 0, 0, t14, -t13, t14, -t11, t13, -t4 * t49 - t29 * t1 + (t32 * t17 + t29 * t5) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-t49 + t32) * t47; 0, 0, 0, 0, 0, 0, 0, -t7, t6, -0.2e1 * t40, -0.2e1 * t19, 0, 0, 0, t28 * t38 - t51, t30 * t38 + t52, t54 * t30 - t50 * t46, -t1, t50 * t23 + t54 * t28, -pkin(7) * t1 + t4 * t18 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, -t47, -t39, 0, 0, 0, 0, 0, -t14, t13, -t14, t11, -t13, -t49 * t9 + (t32 * pkin(7) + t18 * t29) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t15, 0, 0, 0, t28 * t43, t30 * t43, 0.2e1 * t18 * t46 - 0.2e1 * t9 * t30, 0, -0.2e1 * t18 * t23 - 0.2e1 * t9 * t28, 0.2e1 * t18 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t46, 0, -t3, t2, -t3, t8, -t2, -t17 * t8 - t58 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t10, -t12, 0, -t10, -t29 * t8 + t58 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t46, 0, -t41, t42, -t41, -t8, -t42, -t8 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, qJ(5) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t21;
