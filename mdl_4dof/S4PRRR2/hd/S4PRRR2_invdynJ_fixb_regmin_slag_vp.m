% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:25
% EndTime: 2019-07-18 13:27:26
% DurationCPUTime: 0.19s
% Computational Cost: add. (189->72), mult. (274->92), div. (0->0), fcn. (148->10), ass. (0->45)
t27 = cos(qJ(3));
t52 = t27 * pkin(1);
t24 = sin(qJ(3));
t26 = cos(qJ(4));
t51 = t24 * t26;
t20 = qJDD(2) + qJDD(3);
t15 = qJDD(4) + t20;
t50 = t26 * t15;
t49 = pkin(1) * qJD(2);
t48 = qJD(3) * t27;
t23 = sin(qJ(4));
t47 = qJD(4) * t23;
t46 = qJD(4) * t26;
t22 = qJ(2) + qJ(3);
t41 = -qJD(3) - qJD(4);
t16 = qJD(2) - t41;
t45 = -qJD(4) + t16;
t44 = qJDD(2) * t24;
t43 = -t15 - qJDD(2);
t19 = qJ(4) + t22;
t11 = sin(t19);
t12 = cos(t19);
t14 = qJDD(2) * t52;
t39 = t24 * t49;
t2 = t20 * pkin(2) - qJD(3) * t39 + t14;
t42 = g(1) * t12 + g(3) * t11 + t26 * t2;
t17 = sin(t22);
t18 = cos(t22);
t40 = g(1) * t18 + g(3) * t17 + t14;
t38 = t23 * t48;
t37 = -g(1) * t17 + g(3) * t18;
t36 = t23 * t44;
t21 = qJD(2) + qJD(3);
t3 = t21 * pkin(2) + t27 * t49;
t35 = t45 * t3;
t34 = -g(1) * t11 + g(3) * t12 + t39 * t47;
t33 = qJD(2) * (-qJD(3) + t21);
t32 = qJD(3) * (-qJD(2) - t21);
t31 = t41 * t16;
t30 = (-qJD(2) - t16) * t48;
t13 = pkin(2) + t52;
t29 = qJD(4) * (-t13 * t16 - t3);
t28 = cos(qJ(2));
t25 = sin(qJ(2));
t1 = [qJDD(1) + g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, qJDD(2), g(1) * t28 + g(3) * t25, -g(1) * t25 + g(3) * t28, t20, (t20 * t27 + t24 * t32) * pkin(1) + t40, ((-qJDD(2) - t20) * t24 + t27 * t32) * pkin(1) + t37, t15, t13 * t50 + t23 * t29 + (t23 * t30 + (t43 * t23 + (-qJD(2) * qJD(4) + t31) * t26) * t24) * pkin(1) + t42, (-t13 * t15 - t2) * t23 + t26 * t29 + (t26 * t30 + (-t23 * t31 + t43 * t26) * t24) * pkin(1) + t34; 0, 0, 0, 0, t20, t24 * pkin(1) * t33 + t40, (t27 * t33 - t44) * pkin(1) + t37, t15, -t3 * t47 + (-t16 * t47 + t50) * pkin(2) + (-t36 + (-t24 * t46 - t38 - (-t23 * t27 - t51) * t16) * qJD(2)) * pkin(1) + t42, -t3 * t46 - t23 * t2 + (-t23 * t15 - t16 * t46) * pkin(2) + (-t26 * t44 + (-t26 * t48 + (-t23 * t24 + t26 * t27) * t16) * qJD(2)) * pkin(1) + t34; 0, 0, 0, 0, 0, 0, 0, t15, t23 * t35 + (-t36 + (t45 * t51 - t38) * qJD(2)) * pkin(1) + t42, (-t16 * t39 - t2) * t23 + (t35 + (-qJD(2) * t48 - t44) * pkin(1)) * t26 + t34;];
tau_reg  = t1;
