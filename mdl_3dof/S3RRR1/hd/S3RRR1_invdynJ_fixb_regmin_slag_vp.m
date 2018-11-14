% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% tau_reg [3x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3RRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:16:00
% EndTime: 2018-11-14 10:16:00
% DurationCPUTime: 0.20s
% Computational Cost: add. (188->71), mult. (274->92), div. (0->0), fcn. (148->10), ass. (0->45)
t27 = cos(qJ(2));
t52 = t27 * pkin(1);
t22 = qJ(1) + qJ(2);
t17 = sin(t22);
t18 = cos(t22);
t51 = g(1) * t18 + g(2) * t17;
t24 = sin(qJ(2));
t26 = cos(qJ(3));
t50 = t24 * t26;
t20 = qJDD(1) + qJDD(2);
t15 = qJDD(3) + t20;
t49 = t26 * t15;
t48 = pkin(1) * qJD(1);
t47 = qJD(2) * t27;
t23 = sin(qJ(3));
t46 = qJD(3) * t23;
t45 = qJD(3) * t26;
t40 = -qJD(2) - qJD(3);
t16 = qJD(1) - t40;
t44 = -qJD(3) + t16;
t43 = qJDD(1) * t24;
t42 = -t15 - qJDD(1);
t19 = qJ(3) + t22;
t11 = sin(t19);
t12 = cos(t19);
t39 = t24 * t48;
t41 = g(1) * t12 + g(2) * t11 + t39 * t46;
t38 = t23 * t47;
t37 = t23 * t43;
t21 = qJD(1) + qJD(2);
t3 = t21 * pkin(2) + t27 * t48;
t36 = t44 * t3;
t14 = qJDD(1) * t52;
t2 = t20 * pkin(2) - qJD(2) * t39 + t14;
t35 = g(1) * t11 - g(2) * t12 + t26 * t2;
t34 = qJD(1) * (-qJD(2) + t21);
t33 = qJD(2) * (-qJD(1) - t21);
t32 = t40 * t16;
t31 = g(1) * t17 - g(2) * t18 + t14;
t30 = (-qJD(1) - t16) * t47;
t13 = pkin(2) + t52;
t29 = qJD(3) * (-t13 * t16 - t3);
t28 = cos(qJ(1));
t25 = sin(qJ(1));
t1 = [qJDD(1), g(1) * t25 - g(2) * t28, g(1) * t28 + g(2) * t25, t20 (t20 * t27 + t24 * t33) * pkin(1) + t31 ((-qJDD(1) - t20) * t24 + t27 * t33) * pkin(1) + t51, t15, t13 * t49 + t23 * t29 + (t23 * t30 + (t42 * t23 + (-qJD(1) * qJD(3) + t32) * t26) * t24) * pkin(1) + t35 (-t13 * t15 - t2) * t23 + t26 * t29 + (t26 * t30 + (-t23 * t32 + t42 * t26) * t24) * pkin(1) + t41; 0, 0, 0, t20, t24 * pkin(1) * t34 + t31 (t27 * t34 - t43) * pkin(1) + t51, t15, -t3 * t46 + (-t16 * t46 + t49) * pkin(2) + (-t37 + (-t24 * t45 - t38 - (-t23 * t27 - t50) * t16) * qJD(1)) * pkin(1) + t35, -t3 * t45 - t23 * t2 + (-t23 * t15 - t16 * t45) * pkin(2) + (-t26 * t43 + (-t26 * t47 + (-t23 * t24 + t26 * t27) * t16) * qJD(1)) * pkin(1) + t41; 0, 0, 0, 0, 0, 0, t15, t23 * t36 + (-t37 + (t44 * t50 - t38) * qJD(1)) * pkin(1) + t35 (-t16 * t39 - t2) * t23 + (t36 + (-qJD(1) * t47 - t43) * pkin(1)) * t26 + t41;];
tau_reg  = t1;
