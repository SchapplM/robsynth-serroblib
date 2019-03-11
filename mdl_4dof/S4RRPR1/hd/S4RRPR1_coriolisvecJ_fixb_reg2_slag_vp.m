% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:00
% EndTime: 2019-03-08 18:35:01
% DurationCPUTime: 0.23s
% Computational Cost: add. (401->57), mult. (992->94), div. (0->0), fcn. (582->6), ass. (0->48)
t32 = sin(pkin(7));
t37 = cos(qJ(2));
t33 = cos(pkin(7));
t35 = sin(qJ(2));
t55 = t33 * t35;
t40 = pkin(1) * (-t32 * t37 - t55);
t20 = qJD(2) * t40;
t14 = qJD(1) * t20;
t53 = pkin(1) * qJD(1);
t52 = t35 * t53;
t48 = t32 * t52;
t51 = t37 * t53;
t15 = (t33 * t51 - t48) * qJD(2);
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t31 = qJD(1) + qJD(2);
t25 = t31 * pkin(2) + t51;
t12 = t32 * t25 + t33 * t52;
t54 = t34 * t12;
t11 = t33 * t25 - t48;
t9 = t31 * pkin(3) + t11;
t1 = (qJD(4) * t9 + t15) * t36 - qJD(4) * t54 + t34 * t14;
t59 = pkin(2) * t32;
t19 = qJD(1) * t40;
t56 = t32 * t35;
t39 = pkin(1) * (t33 * t37 - t56);
t21 = qJD(1) * t39;
t28 = t33 * pkin(2) + pkin(3);
t41 = t36 * t28 - t34 * t59;
t58 = t41 * qJD(4) - t34 * t19 - t36 * t21;
t42 = t34 * t28 + t36 * t59;
t57 = -t42 * qJD(4) - t36 * t19 + t34 * t21;
t29 = t37 * pkin(1) + pkin(2);
t47 = -pkin(1) * t56 + t33 * t29;
t46 = (-qJD(2) + t31) * t53;
t45 = pkin(1) * qJD(2) * (-qJD(1) - t31);
t6 = t36 * t12 + t34 * t9;
t18 = pkin(3) + t47;
t23 = pkin(1) * t55 + t32 * t29;
t44 = t36 * t18 - t34 * t23;
t43 = t34 * t18 + t36 * t23;
t2 = -t6 * qJD(4) + t36 * t14 - t34 * t15;
t30 = qJD(4) + t31;
t22 = qJD(2) * t39;
t5 = t36 * t9 - t54;
t4 = -t43 * qJD(4) + t36 * t20 - t34 * t22;
t3 = t44 * qJD(4) + t34 * t20 + t36 * t22;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t45, t37 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t31 + t14, -t22 * t31 - t15, 0, t11 * t20 + t12 * t22 + t14 * t47 + t15 * t23, 0, 0, 0, 0, 0, 0, t4 * t30 + t2, -t3 * t30 - t1, 0, t1 * t43 + t2 * t44 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t46, t37 * t46, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t31 + t14, t21 * t31 - t15, 0, -t11 * t19 - t12 * t21 + (t14 * t33 + t15 * t32) * pkin(2), 0, 0, 0, 0, 0, 0, t57 * t30 + t2, -t58 * t30 - t1, 0, t1 * t42 + t2 * t41 + t57 * t5 + t58 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t30 + t2, t5 * t30 - t1, 0, 0;];
tauc_reg  = t7;
