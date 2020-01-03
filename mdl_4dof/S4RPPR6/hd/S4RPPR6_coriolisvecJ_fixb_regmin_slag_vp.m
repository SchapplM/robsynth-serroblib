% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (140->52), mult. (425->92), div. (0->0), fcn. (270->4), ass. (0->44)
t27 = sin(pkin(6));
t30 = cos(qJ(4));
t28 = cos(pkin(6));
t29 = sin(qJ(4));
t49 = t28 * t29;
t14 = t27 * t30 - t49;
t51 = t27 * t28;
t48 = -pkin(5) + qJ(2);
t25 = t27 ^ 2;
t26 = t28 ^ 2;
t47 = t25 + t26;
t46 = qJ(2) * t26;
t13 = t27 * t29 + t28 * t30;
t45 = qJD(1) * t13;
t44 = qJD(1) * t27;
t43 = qJD(3) * t27;
t42 = qJD(1) * qJD(2);
t41 = qJD(1) * t49;
t40 = t30 * t44;
t39 = 0.2e1 * qJD(1) * qJD(3);
t32 = qJD(1) ^ 2;
t16 = t47 * t32;
t38 = t27 * qJ(3) + pkin(1);
t37 = qJ(2) * t42;
t36 = 0.2e1 * t45;
t35 = -t28 * pkin(2) - t38;
t34 = t14 * qJD(2);
t33 = t13 * qJD(2);
t6 = t13 * qJD(4);
t11 = (pkin(2) + pkin(3)) * t28 + t38;
t31 = qJD(4) ^ 2;
t21 = t25 * t37;
t20 = qJ(2) * t44 + qJD(3);
t19 = qJD(4) * t41;
t18 = t48 * t28;
t17 = t48 * t27;
t10 = t40 - t41;
t7 = t14 * qJD(4);
t5 = 0.2e1 * t47 * t42;
t4 = t35 * qJD(1) + qJD(2);
t3 = qJD(4) * t40 - t19;
t2 = qJD(1) * t6;
t1 = t11 * qJD(1) - qJD(2);
t8 = [0, 0, 0, 0, 0, t5, 0.2e1 * t26 * t37 + 0.2e1 * t21, t39 * t51, t5, t25 * t39, t21 + (t20 * qJD(2) - t4 * qJD(3)) * t27 + (0.2e1 * qJD(2) * t46 - t35 * t43) * qJD(1), -t10 * t6 - t2 * t14, -t10 * t7 + t2 * t13 - t14 * t3 + t45 * t6, -t6 * qJD(4), -t7 * qJD(4), 0, t1 * t7 + t11 * t3 + ((-t17 * t29 - t18 * t30) * qJD(4) + t34) * qJD(4) + t36 * t43, -t1 * t6 - t11 * t2 + ((-t17 * t30 + t18 * t29) * qJD(4) - t33) * qJD(4) + (qJD(1) * t14 + t10) * t43; 0, 0, 0, 0, 0, -t16, -qJ(2) * t16, 0, -t16, 0, -t32 * t46 + (-qJD(3) - t20) * t44, 0, 0, 0, 0, 0, t19 + (-t10 - t40) * qJD(4), t36 * qJD(4); 0, 0, 0, 0, 0, 0, 0, -t32 * t51, 0, -t25 * t32, (qJD(2) + t4) * t44, 0, 0, 0, 0, 0, -t31 * t29 - t44 * t45, -t10 * t44 - t31 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t45, t10 ^ 2 - t45 ^ 2, 0, t19 + (t10 - t40) * qJD(4), 0, qJD(1) * t34 - t1 * t10, -qJD(1) * t33 + t1 * t45;];
tauc_reg = t8;
