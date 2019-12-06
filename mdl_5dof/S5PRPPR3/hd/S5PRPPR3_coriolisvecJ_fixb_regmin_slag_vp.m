% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:55
% DurationCPUTime: 0.16s
% Computational Cost: add. (171->55), mult. (442->85), div. (0->0), fcn. (304->6), ass. (0->47)
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t14 = t23 * t28 + t24 * t26;
t10 = t14 * qJD(1);
t39 = t28 * qJD(1);
t16 = qJD(2) * pkin(2) + t39;
t42 = qJD(1) * t26;
t6 = t23 * t16 + t24 * t42;
t3 = qJD(2) * qJ(4) + t6;
t51 = -t3 + t10;
t48 = t24 * t28;
t49 = t23 * t26;
t13 = -t48 + t49;
t9 = t14 * qJD(2);
t7 = qJD(1) * t9;
t50 = t7 * t13;
t25 = sin(qJ(5));
t27 = cos(qJ(5));
t47 = t25 * t27;
t29 = qJD(5) ^ 2;
t46 = t29 * t25;
t45 = t29 * t27;
t44 = t25 ^ 2 - t27 ^ 2;
t30 = qJD(2) ^ 2;
t43 = -t29 - t30;
t41 = t25 * qJD(5);
t40 = t27 * qJD(5);
t17 = t23 * t42;
t36 = t24 * t39;
t12 = -t17 + t36;
t38 = qJD(4) - t12;
t37 = qJD(2) * qJD(5);
t35 = qJD(2) * t48;
t34 = -t24 * pkin(2) - pkin(3);
t33 = t3 * qJD(2) - t7;
t5 = t24 * t16 - t17;
t20 = t23 * pkin(2) + qJ(4);
t32 = qJD(2) * t20 - t51;
t15 = qJD(2) * t17;
t4 = -t15 + (qJD(4) + t36) * qJD(2);
t31 = t38 * qJD(2) - (-pkin(6) + t34) * t29 + t4;
t11 = -qJD(2) * t49 + t35;
t8 = qJD(1) * t35 - t15;
t2 = -qJD(2) * pkin(3) + qJD(4) - t5;
t1 = [0, 0, -t30 * t26, -t30 * t28, t6 * t11 + t8 * t14 - t5 * t9 + t50, t9 * qJD(2), t11 * qJD(2), t3 * t11 + t4 * t14 + t2 * t9 + t50, 0, 0, 0, 0, 0, t9 * t40 - t13 * t46 + (t11 * t25 + t14 * t40) * qJD(2), -t9 * t41 - t13 * t45 + (t11 * t27 - t14 * t41) * qJD(2); 0, 0, 0, 0, t5 * t10 - t6 * t12 + (t23 * t8 - t24 * t7) * pkin(2), 0, -t15 + (0.2e1 * qJD(4) - t12 + t36) * qJD(2), -t2 * t10 + t4 * t20 + t38 * t3 + t7 * t34, -0.2e1 * t37 * t47, 0.2e1 * t44 * t37, -t46, -t45, 0, t31 * t25 + t32 * t40, t31 * t27 - t32 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t46; 0, 0, 0, 0, 0, 0, -t30, t51 * qJD(2), 0, 0, 0, 0, 0, t43 * t25, t43 * t27; 0, 0, 0, 0, 0, 0, 0, 0, t30 * t47, -t44 * t30, 0, 0, 0, -t33 * t27, t33 * t25;];
tauc_reg = t1;
