% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:18
% DurationCPUTime: 0.20s
% Computational Cost: add. (168->49), mult. (354->80), div. (0->0), fcn. (130->2), ass. (0->41)
t15 = sin(qJ(3));
t17 = -pkin(1) - pkin(5);
t11 = t17 * qJD(1) + qJD(2);
t16 = cos(qJ(3));
t38 = t16 * t11;
t3 = (qJD(4) + t38) * qJD(3);
t33 = qJD(3) * pkin(3);
t24 = -qJD(4) + t33;
t4 = -t24 - t38;
t29 = qJD(3) * qJ(4);
t40 = t15 * t11;
t6 = t29 + t40;
t42 = t16 * t6;
t45 = ((-t4 + t38) * t15 - t42) * qJD(3) - t3 * t15;
t27 = 2 * qJD(1);
t43 = 0.2e1 * qJD(3);
t39 = t15 * t16;
t18 = qJD(3) ^ 2;
t37 = t18 * t15;
t36 = t18 * t16;
t14 = t16 ^ 2;
t35 = t15 ^ 2 - t14;
t19 = qJD(1) ^ 2;
t34 = t18 + t19;
t23 = pkin(3) * t16 + qJ(4) * t15;
t2 = t23 * qJD(3) - t16 * qJD(4) + qJD(2);
t1 = qJD(1) * t2;
t8 = t15 * pkin(3) - t16 * qJ(4) + qJ(2);
t5 = qJD(1) * t8;
t32 = t19 * qJ(2);
t31 = t5 * qJD(1);
t30 = qJ(2) * qJD(3);
t28 = qJD(1) * qJD(3);
t26 = qJD(2) * t27;
t22 = t5 * t43;
t21 = -t17 * t18 + 0.2e1 * t1;
t12 = t19 * t39;
t10 = t34 * t16;
t9 = t34 * t15;
t7 = t23 * qJD(1);
t13 = [0, 0, 0, 0, t26, qJ(2) * t26, -0.2e1 * t28 * t39, 0.2e1 * t35 * t28, -t37, -t36, 0, -t17 * t37 + (qJD(2) * t15 + t16 * t30) * t27, -t17 * t36 + (qJD(2) * t16 - t15 * t30) * t27, t21 * t15 + t16 * t22, t45, t15 * t22 - t21 * t16, t1 * t8 - t17 * t45 + t5 * t2; 0, 0, 0, 0, -t19, -t32, 0, 0, 0, 0, 0, -t9, -t10, -t9, 0, t10, -t45 - t31; 0, 0, 0, 0, 0, 0, t12, -t35 * t19, 0, 0, 0, -t16 * t32, t15 * t32, (-t15 * t7 - t16 * t5) * qJD(1), ((t6 - t29) * t16 + (t24 + t4) * t15) * qJD(1), qJD(4) * t43 + (-t15 * t5 + t16 * t7) * qJD(1), t3 * qJ(4) + t6 * qJD(4) - t5 * t7 + (-t42 + (-t4 - t33) * t15) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t14 * t19 - t18, t16 * t31 + (-t6 + t40) * qJD(3);];
tauc_reg = t13;
