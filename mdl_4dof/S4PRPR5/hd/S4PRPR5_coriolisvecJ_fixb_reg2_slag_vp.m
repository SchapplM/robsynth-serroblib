% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.24s
% Computational Cost: add. (284->55), mult. (749->100), div. (0->0), fcn. (535->6), ass. (0->48)
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t15 = t23 * t26 - t24 * t28;
t14 = t15 * qJD(1);
t16 = t23 * t28 + t24 * t26;
t12 = t16 * qJD(1);
t11 = t16 * qJD(2);
t9 = qJD(1) * t11;
t48 = t9 * t15;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t47 = t25 * t27;
t29 = qJD(4) ^ 2;
t46 = t29 * t25;
t45 = t29 * t27;
t21 = t25 ^ 2;
t22 = t27 ^ 2;
t44 = t21 - t22;
t43 = t21 + t22;
t42 = qJD(1) * t26;
t41 = qJD(4) * t25;
t40 = qJD(4) * t27;
t13 = t15 * qJD(2);
t39 = t13 * qJD(2);
t38 = qJD(2) * qJD(4);
t30 = qJD(2) ^ 2;
t37 = t30 * t47;
t10 = qJD(1) * t13;
t18 = qJD(2) * pkin(2) + t28 * qJD(1);
t7 = t24 * t18 - t23 * t42;
t5 = -qJD(2) * pkin(3) - t7;
t36 = -qJD(2) * t5 + t10;
t35 = t38 * t47;
t8 = t23 * t18 + t24 * t42;
t6 = qJD(2) * pkin(5) + t8;
t3 = t27 * qJD(3) - t25 * t6;
t4 = t25 * qJD(3) + t27 * t6;
t34 = t25 * t3 - t27 * t4;
t19 = t23 * pkin(2) + pkin(5);
t33 = t12 * qJD(2) - t19 * t29 - t9;
t20 = -t24 * pkin(2) - pkin(3);
t32 = qJD(4) * (qJD(2) * t20 - t14 + t5);
t1 = t3 * qJD(4) - t27 * t10;
t2 = -t4 * qJD(4) + t25 * t10;
t31 = t1 * t27 - t2 * t25 + (-t25 * t4 - t27 * t3) * qJD(4);
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t26, -t30 * t28, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * qJD(2), t39, 0, -t10 * t16 - t7 * t11 - t8 * t13 + t48, 0, 0, 0, 0, 0, 0, t13 * t41 - t16 * t45 + (-t11 * t27 + t15 * t41) * qJD(2), t13 * t40 + t16 * t46 + (t11 * t25 + t15 * t40) * qJD(2), -t43 * t39, t5 * t11 + t34 * t13 + t31 * t16 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t12 + t8 * t14 + (-t10 * t23 - t24 * t9) * pkin(2), 0.2e1 * t35, -0.2e1 * t44 * t38, t45, -0.2e1 * t35, -t46, 0, t25 * t32 + t33 * t27, -t33 * t25 + t27 * t32, t43 * t14 * qJD(2) + t31, -t5 * t12 - t34 * t14 + t31 * t19 + t9 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, -t34 * qJD(4) + t1 * t25 + t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t44 * t30, 0, t37, 0, 0, t36 * t25, t36 * t27, 0, 0;];
tauc_reg = t17;
