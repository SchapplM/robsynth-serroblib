% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:01
% EndTime: 2021-01-15 11:24:03
% DurationCPUTime: 0.29s
% Computational Cost: add. (385->61), mult. (771->102), div. (0->0), fcn. (631->4), ass. (0->39)
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t50 = sin(pkin(7));
t42 = qJD(3) * t50;
t51 = cos(pkin(7));
t43 = qJD(3) * t51;
t16 = t32 * t42 - t33 * t43;
t17 = -t32 * t43 - t33 * t42;
t19 = -t50 * t32 + t51 * t33;
t20 = t51 * t32 + t50 * t33;
t60 = 0.2e1 * t20 * t16 - 0.2e1 * t19 * t17;
t59 = (-t50 * t16 + t51 * t17) * pkin(3);
t55 = -pkin(1) - pkin(6);
t46 = qJ(4) - t55;
t21 = t46 * t32;
t41 = t46 * t33;
t11 = -t50 * t21 + t51 * t41;
t12 = -t51 * t21 - t50 * t41;
t15 = -qJD(3) * t41 - t32 * qJD(4);
t49 = t32 * qJD(3);
t34 = -t33 * qJD(4) + t46 * t49;
t8 = t50 * t15 - t51 * t34;
t9 = t51 * t15 + t50 * t34;
t38 = -t11 * t17 - t12 * t16 - t8 * t19 + t9 * t20;
t57 = 2 * qJD(2);
t56 = 2 * qJD(5);
t28 = t32 * pkin(3) + qJ(2);
t48 = t33 * qJD(3);
t22 = pkin(3) * t48 + qJD(2);
t47 = qJ(2) * qJD(3);
t45 = t11 * t8 + t12 * t9;
t44 = t55 * qJD(3);
t26 = t50 * pkin(3) + qJ(5);
t29 = -t51 * pkin(3) - pkin(4);
t37 = t20 * qJD(5) - t26 * t16 - t29 * t17;
t35 = 0.2e1 * t38;
t10 = t20 * pkin(4) - t19 * qJ(5) + t28;
t5 = -t16 * pkin(4) - t17 * qJ(5) - t19 * qJD(5) + t22;
t1 = [0, 0, 0, 0, t57, qJ(2) * t57, -0.2e1 * t32 * t48, 0.2e1 * (t32 ^ 2 - t33 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t32 + 0.2e1 * t33 * t47, 0.2e1 * qJD(2) * t33 - 0.2e1 * t32 * t47, -0.2e1 * t28 * t16 + 0.2e1 * t22 * t20, 0.2e1 * t28 * t17 + 0.2e1 * t22 * t19, -t35, 0.2e1 * t28 * t22 + 0.2e1 * t45, -0.2e1 * t10 * t16 + 0.2e1 * t5 * t20, -t35, -0.2e1 * t10 * t17 - 0.2e1 * t5 * t19, 0.2e1 * t10 * t5 + 0.2e1 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t38, 0, t60, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, 0, -t32 * t44, -t33 * t44, -t8, -t9, -t59, (t50 * t9 - t51 * t8) * pkin(3), -t8, -t37, t9, t12 * qJD(5) + t9 * t26 + t8 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t48, t17, t16, 0, t59, t17, 0, -t16, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t26 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t17, 0, t22, -t16, 0, -t17, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
