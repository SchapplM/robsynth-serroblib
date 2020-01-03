% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP6
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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:15
% DurationCPUTime: 0.22s
% Computational Cost: add. (153->47), mult. (320->84), div. (0->0), fcn. (124->2), ass. (0->43)
t20 = -pkin(1) - pkin(5);
t12 = t20 * qJD(1) + qJD(2);
t48 = -qJ(4) * qJD(1) + t12;
t18 = sin(qJ(3));
t29 = qJ(4) * qJD(3);
t19 = cos(qJ(3));
t32 = t19 * qJD(4);
t35 = qJD(3) * t18;
t1 = -t12 * t35 + (t18 * t29 - t32) * qJD(1);
t33 = t18 * qJD(4);
t34 = qJD(3) * t19;
t2 = t12 * t34 + (-t19 * t29 - t33) * qJD(1);
t39 = qJD(3) * pkin(3);
t7 = t48 * t19;
t3 = t7 + t39;
t6 = t48 * t18;
t47 = -t1 * t19 - t2 * t18 + (t18 * t3 - t19 * t6) * qJD(3);
t15 = qJD(1) * qJD(2);
t46 = 0.2e1 * t15;
t45 = t3 - t7;
t21 = qJD(3) ^ 2;
t44 = t21 * t18;
t43 = t21 * t19;
t28 = qJD(1) * qJD(3);
t26 = t19 * t28;
t42 = pkin(3) * t26 + t15;
t16 = t18 ^ 2;
t17 = t19 ^ 2;
t41 = t16 - t17;
t22 = qJD(1) ^ 2;
t40 = -t21 - t22;
t38 = t22 * qJ(2);
t25 = t18 * pkin(3) + qJ(2);
t9 = t25 * qJD(1) + qJD(4);
t37 = t9 * qJD(1);
t36 = qJ(4) - t20;
t31 = qJ(2) * qJD(3);
t27 = 0.2e1 * qJD(1);
t11 = t36 * t19;
t10 = t36 * t18;
t5 = -qJD(3) * t11 - t33;
t4 = t36 * t35 - t32;
t8 = [0, 0, 0, 0, t46, qJ(2) * t46, -0.2e1 * t18 * t26, 0.2e1 * t41 * t28, -t44, -t43, 0, -t20 * t44 + (qJD(2) * t18 + t19 * t31) * t27, -t20 * t43 + (qJD(2) * t19 - t18 * t31) * t27, (-t18 * t5 - t19 * t4 + (t10 * t19 - t11 * t18) * qJD(3)) * qJD(1) + t47, -t2 * t10 + t6 * t5 - t1 * t11 + t3 * t4 + t42 * t25 + t9 * (pkin(3) * t34 + qJD(2)); 0, 0, 0, 0, -t22, -t38, 0, 0, 0, 0, 0, t40 * t18, t40 * t19, 0, -t37 - t47; 0, 0, 0, 0, 0, 0, t19 * t22 * t18, -t41 * t22, 0, 0, 0, -t19 * t38, t18 * t38, (t39 - t45) * t18 * qJD(1), t45 * t6 + (-t19 * t37 + t1) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t16 - t17) * t22, (t18 * t6 + t19 * t3) * qJD(1) + t42;];
tauc_reg = t8;
