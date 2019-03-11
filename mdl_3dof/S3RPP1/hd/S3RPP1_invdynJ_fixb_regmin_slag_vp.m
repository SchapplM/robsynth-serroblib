% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S3RPP1
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% tau_reg [3x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S3RPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:02
% EndTime: 2019-03-08 18:05:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (66->40), mult. (82->36), div. (0->0), fcn. (28->2), ass. (0->22)
t12 = -pkin(1) - qJ(3);
t26 = t12 * qJDD(1);
t10 = (qJD(1) * qJD(2));
t5 = 2 * t10;
t13 = sin(qJ(1));
t14 = cos(qJ(1));
t25 = t14 * pkin(1) + t13 * qJ(2);
t24 = g(1) * t13 - g(2) * t14;
t23 = qJDD(1) * pkin(1);
t11 = qJ(2) * qJDD(1);
t22 = (qJD(3) * qJD(1));
t21 = -qJDD(2) + t24;
t20 = qJDD(3) + t10 + t11;
t19 = qJDD(2) - t23;
t18 = g(1) * t14 + g(2) * t13;
t17 = qJDD(2) + t26;
t16 = -t18 + t5 + 0.2e1 * t11;
t15 = qJD(1) ^ 2;
t4 = t14 * qJ(2);
t2 = qJ(2) * qJD(1) + qJD(3);
t1 = qJD(1) * t12 + qJD(2);
t3 = [qJDD(1), t24, t18, -t21 - 0.2e1 * t23, t16, -t19 * pkin(1) - g(1) * (-t13 * pkin(1) + t4) - g(2) * t25 + (t5 + t11) * qJ(2), qJDD(3) + t16, t21 + (2 * t22) - 0.2e1 * t26 (t17 - t22) * t12 - t1 * qJD(3) + t20 * qJ(2) + t2 * qJD(2) - g(1) * (t12 * t13 + t4) - g(2) * (t14 * qJ(3) + t25); 0, 0, 0, qJDD(1), -t15, -t15 * qJ(2) + t19 - t24, -t15, -qJDD(1) (-qJD(3) - t2) * qJD(1) + t17 - t24; 0, 0, 0, 0, 0, 0, qJDD(1), -t15, t1 * qJD(1) - t18 + t20;];
tau_reg  = t3;
