% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR6
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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:42
% DurationCPUTime: 0.19s
% Computational Cost: add. (171->57), mult. (503->94), div. (0->0), fcn. (342->6), ass. (0->43)
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t12 = t27 * t23 + t25 * t24;
t9 = t12 * qJD(4);
t26 = sin(qJ(2));
t41 = t26 * qJD(1);
t18 = qJD(2) * qJ(3) + t41;
t43 = t23 ^ 2 + t24 ^ 2;
t49 = t43 * t18 - t41;
t48 = t25 * t23;
t47 = t27 * t24;
t29 = qJD(2) ^ 2;
t46 = t29 * t26;
t28 = cos(qJ(2));
t45 = t29 * t28;
t44 = pkin(5) + qJ(3);
t42 = qJD(2) * pkin(2);
t7 = qJD(2) * t12;
t40 = t26 * qJD(4) ^ 2;
t39 = t28 * qJD(1);
t38 = qJD(2) * t48;
t37 = qJD(2) * t47;
t20 = -t24 * pkin(3) - pkin(2);
t36 = t43 * t28;
t13 = (qJD(3) + t39) * qJD(2);
t35 = t43 * t13;
t33 = t43 * qJD(3);
t32 = qJD(3) - t39;
t11 = -t47 + t48;
t8 = t11 * qJD(4);
t31 = t28 * t9;
t30 = t28 * t8;
t17 = qJD(4) * t37;
t16 = t32 - t42;
t15 = t44 * t24;
t14 = t44 * t23;
t10 = t20 * qJD(2) + t32;
t5 = -t37 + t38;
t4 = qJD(2) * t9;
t3 = -qJD(4) * t38 + t17;
t1 = [0, 0, -t46, -t45, -t24 * t46, t23 * t46, t43 * t45, t26 * t35 + (t16 * t26 + t49 * t28) * qJD(2), 0, 0, 0, 0, 0, -t28 * t4 + t11 * t40 + (t26 * t5 - t31) * qJD(2), -t28 * t3 + t12 * t40 + (t26 * t7 + t30) * qJD(2); 0, 0, 0, 0, 0, 0, t35 + (-qJD(1) * t36 + t33) * qJD(2), t18 * t33 + qJ(3) * t35 + ((-t16 - t42) * t26 - t18 * t36) * qJD(1), t3 * t12 - t7 * t8, -t3 * t11 - t12 * t4 + t8 * t5 - t7 * t9, -t8 * qJD(4), -t9 * qJD(4), 0, t10 * t9 + t20 * t4 + ((t14 * t25 - t15 * t27) * qJD(4) - t12 * qJD(3)) * qJD(4) + (t31 + (qJD(2) * t11 - t5) * t26) * qJD(1), -t10 * t8 + t20 * t3 + ((t14 * t27 + t15 * t25) * qJD(4) + t11 * qJD(3)) * qJD(4) - t30 * qJD(1); 0, 0, 0, 0, 0, 0, -t43 * t29, -t49 * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t7 * qJD(4), t17 + (-t5 - t38) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, t7 * t5, -t5 ^ 2 + t7 ^ 2, t17 + (t5 - t38) * qJD(4), 0, 0, -t10 * t7 - t12 * t13, t10 * t5 + t11 * t13;];
tauc_reg = t1;
