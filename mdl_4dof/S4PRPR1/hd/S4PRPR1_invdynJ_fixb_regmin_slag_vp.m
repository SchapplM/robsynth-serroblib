% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:08
% EndTime: 2019-03-08 18:21:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (141->46), mult. (153->47), div. (0->0), fcn. (86->4), ass. (0->31)
t12 = pkin(6) + qJ(2);
t10 = cos(t12);
t9 = sin(t12);
t44 = g(1) * t9 - g(2) * t10;
t43 = -qJDD(3) + t44;
t37 = pkin(2) * qJDD(2);
t42 = t37 + t43;
t17 = -pkin(2) - pkin(3);
t3 = t17 * qJDD(2) + qJDD(3);
t30 = qJD(2) - qJD(4);
t35 = qJD(4) + t30;
t41 = t35 * qJ(3) * qJD(2) - t3;
t11 = qJDD(2) - qJDD(4);
t31 = qJD(2) * qJD(3);
t7 = t17 * qJD(2) + qJD(3);
t40 = (qJD(4) * t17 + qJD(3)) * t30 + t31 + qJD(4) * t7 + (t11 + qJDD(2)) * qJ(3);
t39 = (qJD(2) + t30) * qJ(3) * qJD(4) - t17 * t11 - t3;
t32 = qJ(3) * qJDD(2);
t27 = t30 ^ 2;
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t1 = -t10 * t16 - t9 * t15;
t2 = t10 * t15 - t9 * t16;
t25 = -g(1) * t2 + g(2) * t1;
t24 = -g(1) * t1 - g(2) * t2;
t23 = g(1) * t10 + g(2) * t9;
t21 = -t23 + 0.2e1 * t31;
t20 = -t35 * t7 - t31 - t32;
t18 = qJD(2) ^ 2;
t14 = qJDD(1) - g(3);
t4 = [t14, 0, 0, 0, 0, 0, t14, 0, 0, 0; 0, qJDD(2), t44, t23, 0.2e1 * t37 + t43, t21 + 0.2e1 * t32, t42 * pkin(2) + (t21 + t32) * qJ(3), t11, t40 * t15 + t39 * t16 + t25, -t39 * t15 + t40 * t16 - t24; 0, 0, 0, 0, -qJDD(2), -t18, -t18 * qJ(3) - t42, 0, -t16 * t11 - t15 * t27, t15 * t11 - t16 * t27; 0, 0, 0, 0, 0, 0, 0, -t11, t20 * t15 - t41 * t16 - t25, t41 * t15 + t20 * t16 + t24;];
tau_reg  = t4;
