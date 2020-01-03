% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.13s
% Computational Cost: add. (102->36), mult. (289->66), div. (0->0), fcn. (204->6), ass. (0->34)
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t39 = t22 * t24;
t26 = qJD(4) ^ 2;
t38 = t26 * t22;
t37 = t26 * t24;
t36 = t22 ^ 2 - t24 ^ 2;
t23 = sin(qJ(2));
t35 = qJD(1) * t23;
t34 = t22 * qJD(4);
t33 = t24 * qJD(4);
t25 = cos(qJ(2));
t32 = t25 * qJD(1);
t31 = 0.2e1 * qJD(2) * qJD(4);
t14 = qJD(2) * pkin(2) + t32;
t20 = sin(pkin(7));
t21 = cos(pkin(7));
t3 = t21 * t14 - t20 * t35;
t1 = -qJD(2) * pkin(3) - t3;
t11 = t20 * t23 - t21 * t25;
t9 = t11 * qJD(2);
t6 = qJD(1) * t9;
t30 = -t1 * qJD(2) + t6;
t12 = t20 * t25 + t21 * t23;
t7 = t12 * qJD(2);
t5 = qJD(1) * t7;
t15 = t21 * t35;
t8 = t20 * t32 + t15;
t29 = qJD(2) * t8 - (t20 * pkin(2) + pkin(5)) * t26 - t5;
t10 = t11 * qJD(1);
t28 = qJD(4) * (qJD(2) * (-t21 * pkin(2) - pkin(3)) + t1 - t10);
t27 = qJD(2) ^ 2;
t4 = t20 * t14 + t15;
t2 = [0, 0, -t27 * t23, -t27 * t25, t5 * t11 - t6 * t12 - t3 * t7 - t4 * t9, 0, 0, 0, 0, 0, t9 * t34 - t12 * t37 + (t11 * t34 - t24 * t7) * qJD(2), t9 * t33 + t12 * t38 + (t11 * t33 + t22 * t7) * qJD(2); 0, 0, 0, 0, t4 * t10 + t3 * t8 + (-t20 * t6 - t21 * t5) * pkin(2), t31 * t39, -t36 * t31, t37, -t38, 0, t22 * t28 + t29 * t24, -t29 * t22 + t24 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37; 0, 0, 0, 0, 0, -t27 * t39, t36 * t27, 0, 0, 0, t30 * t22, t30 * t24;];
tauc_reg = t2;
