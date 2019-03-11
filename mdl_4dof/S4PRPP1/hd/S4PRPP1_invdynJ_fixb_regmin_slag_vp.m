% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPP1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:06
% EndTime: 2019-03-08 18:18:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (97->42), mult. (82->36), div. (0->0), fcn. (28->2), ass. (0->24)
t16 = -pkin(2) - qJ(4);
t28 = t16 * qJDD(2);
t13 = (qJD(2) * qJD(3));
t10 = 2 * t13;
t12 = pkin(5) + qJ(2);
t8 = sin(t12);
t9 = cos(t12);
t27 = t9 * pkin(2) + t8 * qJ(3);
t26 = g(1) * t8 - g(2) * t9;
t25 = qJDD(2) * pkin(2);
t14 = qJ(3) * qJDD(2);
t24 = (qJD(4) * qJD(2));
t23 = -qJDD(3) + t26;
t22 = qJDD(4) + t13 + t14;
t21 = qJDD(3) - t25;
t20 = g(1) * t9 + g(2) * t8;
t19 = qJDD(3) + t28;
t18 = t10 + 0.2e1 * t14 - t20;
t17 = qJD(2) ^ 2;
t15 = qJDD(1) - g(3);
t7 = qJ(3) * qJD(2) + qJD(4);
t3 = t9 * qJ(3);
t1 = qJD(2) * t16 + qJD(3);
t2 = [t15, 0, 0, 0, 0, 0, t15, 0, 0, t15; 0, qJDD(2), t26, t20, -t23 - 0.2e1 * t25, t18, -t21 * pkin(2) - g(1) * (-t8 * pkin(2) + t3) - g(2) * t27 + (t10 + t14) * qJ(3), qJDD(4) + t18, t23 + (2 * t24) - 0.2e1 * t28 (t19 - t24) * t16 - t1 * qJD(4) + t22 * qJ(3) + t7 * qJD(3) - g(1) * (t16 * t8 + t3) - g(2) * (t9 * qJ(4) + t27); 0, 0, 0, 0, qJDD(2), -t17, -t17 * qJ(3) + t21 - t26, -t17, -qJDD(2) (-qJD(4) - t7) * qJD(2) + t19 - t26; 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t17, t1 * qJD(2) - t20 + t22;];
tau_reg  = t2;
