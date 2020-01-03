% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauc_reg [4x13]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (117->40), mult. (321->72), div. (0->0), fcn. (150->2), ass. (0->31)
t24 = qJD(3) * pkin(3);
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t26 = -qJ(4) - pkin(5);
t9 = t26 * t14;
t4 = t15 * qJD(1) + qJD(2) * t9;
t3 = t4 + t24;
t31 = t3 - t4;
t10 = t26 * t15;
t6 = t14 * qJD(1) - qJD(2) * t10;
t30 = t15 * t6;
t17 = qJD(2) ^ 2;
t29 = t15 * t17;
t16 = qJD(3) ^ 2;
t28 = t16 * t14;
t27 = t16 * t15;
t12 = t14 ^ 2;
t13 = t15 ^ 2;
t25 = t12 - t13;
t23 = qJD(2) * qJD(3);
t22 = qJD(3) * qJD(1);
t21 = 0.2e1 * t23;
t20 = qJD(3) * t26;
t19 = -0.2e1 * pkin(2) * t23;
t18 = (-t15 * pkin(3) - pkin(2)) * qJD(2);
t7 = -t14 * qJD(4) + t15 * t20;
t5 = t15 * qJD(4) + t14 * t20;
t8 = qJD(4) + t18;
t2 = qJD(2) * t7 - t14 * t22;
t1 = qJD(2) * t5 + t15 * t22;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t27, 0, t1 * t14 + t2 * t15 + (-t14 * t3 + t30) * qJD(3); 0, 0, 0, 0, t14 * t15 * t21, -t25 * t21, t27, -t28, 0, -pkin(5) * t27 + t14 * t19, pkin(5) * t28 + t15 * t19, t1 * t15 - t2 * t14 + (-t14 * t6 - t15 * t3) * qJD(3) + (-t14 * t7 + t15 * t5 + (t10 * t14 - t15 * t9) * qJD(3)) * qJD(2), -t1 * t10 + t2 * t9 + t3 * t7 + t6 * t5 + (t8 + t18) * t14 * t24; 0, 0, 0, 0, -t14 * t29, t25 * t17, 0, 0, 0, t17 * pkin(2) * t14, pkin(2) * t29, (-t24 + t31) * t15 * qJD(2), t31 * t6 + (-qJD(2) * t14 * t8 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t12 - t13) * t17, (-t30 + (t3 + t24) * t14) * qJD(2);];
tauc_reg = t11;
