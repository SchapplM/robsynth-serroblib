% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:21
% EndTime: 2019-03-08 18:25:21
% DurationCPUTime: 0.20s
% Computational Cost: add. (213->73), mult. (274->92), div. (0->0), fcn. (148->10), ass. (0->46)
t29 = cos(qJ(3));
t53 = t29 * pkin(2);
t24 = pkin(7) + qJ(2);
t22 = qJ(3) + t24;
t13 = sin(t22);
t14 = cos(t22);
t52 = g(1) * t14 + g(2) * t13;
t27 = sin(qJ(3));
t28 = cos(qJ(4));
t51 = t27 * t28;
t23 = qJDD(2) + qJDD(3);
t18 = qJDD(4) + t23;
t50 = t28 * t18;
t49 = pkin(2) * qJD(2);
t48 = qJD(3) * t29;
t26 = sin(qJ(4));
t47 = qJD(4) * t26;
t46 = qJD(4) * t28;
t41 = -qJD(3) - qJD(4);
t21 = qJD(2) - t41;
t45 = -qJD(4) + t21;
t44 = qJDD(2) * t27;
t43 = -t18 - qJDD(2);
t15 = qJ(4) + t22;
t11 = sin(t15);
t12 = cos(t15);
t40 = t27 * t49;
t42 = g(1) * t12 + g(2) * t11 + t40 * t47;
t39 = t26 * t48;
t38 = t26 * t44;
t25 = qJD(2) + qJD(3);
t3 = pkin(3) * t25 + t29 * t49;
t37 = t45 * t3;
t17 = qJDD(2) * t53;
t2 = pkin(3) * t23 - qJD(3) * t40 + t17;
t36 = g(1) * t11 - g(2) * t12 + t2 * t28;
t35 = qJD(2) * (-qJD(3) + t25);
t34 = qJD(3) * (-qJD(2) - t25);
t33 = t41 * t21;
t32 = g(1) * t13 - g(2) * t14 + t17;
t31 = (-qJD(2) - t21) * t48;
t16 = pkin(3) + t53;
t30 = qJD(4) * (-t16 * t21 - t3);
t20 = cos(t24);
t19 = sin(t24);
t1 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, qJDD(2), g(1) * t19 - g(2) * t20, g(1) * t20 + g(2) * t19, t23 (t23 * t29 + t27 * t34) * pkin(2) + t32 ((-qJDD(2) - t23) * t27 + t29 * t34) * pkin(2) + t52, t18, t16 * t50 + t26 * t30 + (t26 * t31 + (t43 * t26 + (-qJD(2) * qJD(4) + t33) * t28) * t27) * pkin(2) + t36 (-t16 * t18 - t2) * t26 + t28 * t30 + (t28 * t31 + (-t26 * t33 + t28 * t43) * t27) * pkin(2) + t42; 0, 0, 0, 0, t23, pkin(2) * t27 * t35 + t32 (t29 * t35 - t44) * pkin(2) + t52, t18, -t3 * t47 + (-t21 * t47 + t50) * pkin(3) + (-t38 + (-t27 * t46 - t39 - (-t26 * t29 - t51) * t21) * qJD(2)) * pkin(2) + t36, -t3 * t46 - t26 * t2 + (-t26 * t18 - t21 * t46) * pkin(3) + (-t28 * t44 + (-t28 * t48 + (-t26 * t27 + t28 * t29) * t21) * qJD(2)) * pkin(2) + t42; 0, 0, 0, 0, 0, 0, 0, t18, t26 * t37 + (-t38 + (t45 * t51 - t39) * qJD(2)) * pkin(2) + t36 (-t21 * t40 - t2) * t26 + (t37 + (-qJD(2) * t48 - t44) * pkin(2)) * t28 + t42;];
tau_reg  = t1;
