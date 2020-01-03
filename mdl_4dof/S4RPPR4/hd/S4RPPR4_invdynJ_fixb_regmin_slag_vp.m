% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:55
% EndTime: 2019-12-31 16:38:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (144->49), mult. (237->71), div. (0->0), fcn. (124->8), ass. (0->41)
t33 = qJD(3) * qJD(1);
t17 = sin(pkin(6));
t8 = t17 * pkin(1) + qJ(3);
t38 = t8 * qJDD(1);
t3 = t33 + t38;
t13 = qJ(1) + pkin(6);
t12 = cos(t13);
t9 = g(2) * t12;
t46 = qJDD(3) + t9;
t18 = cos(pkin(6));
t10 = -t18 * pkin(1) - pkin(2);
t45 = t10 * qJDD(1);
t6 = qJD(1) * t8;
t7 = -pkin(5) + t10;
t44 = 0.2e1 * qJD(4) * t6 + qJDD(4) * t7;
t11 = sin(t13);
t43 = g(1) * t11;
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t42 = t19 * t21;
t15 = t21 ^ 2;
t41 = t19 ^ 2 - t15;
t23 = qJD(4) ^ 2;
t24 = qJD(1) ^ 2;
t40 = -t23 - t24;
t16 = qJDD(2) - g(3);
t37 = qJDD(4) * t19;
t36 = qJDD(4) * t21;
t35 = t21 * qJDD(1);
t34 = qJD(1) * qJD(4);
t31 = -g(1) * t12 - g(2) * t11;
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t30 = g(1) * t20 - g(2) * t22;
t29 = -t6 * qJD(1) - t43;
t28 = qJDD(3) + t45;
t27 = -t7 * qJDD(1) - t29 - t46;
t26 = -t23 * t7 + 0.2e1 * t3 + t31;
t5 = -t23 * t19 + t36;
t4 = -t23 * t21 - t37;
t1 = [qJDD(1), t30, g(1) * t22 + g(2) * t20, (t30 + (t17 ^ 2 + t18 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), -t43 + 0.2e1 * t45 + t46, t31 + 0.2e1 * t33 + 0.2e1 * t38, t3 * t8 + t6 * qJD(3) + t28 * t10 - g(1) * (-t20 * pkin(1) - t11 * pkin(2) + t12 * qJ(3)) - g(2) * (t22 * pkin(1) + t12 * pkin(2) + t11 * qJ(3)), t15 * qJDD(1) - 0.2e1 * t34 * t42, -0.2e1 * t19 * t35 + 0.2e1 * t41 * t34, t5, t4, 0, t26 * t19 + t44 * t21, -t44 * t19 + t26 * t21; 0, 0, 0, t16, 0, 0, t16, 0, 0, 0, 0, 0, t4, -t5; 0, 0, 0, 0, qJDD(1), -t24, t28 + t29 + t9, 0, 0, 0, 0, 0, t40 * t19 + t36, t40 * t21 - t37; 0, 0, 0, 0, 0, 0, 0, t24 * t42, -t41 * t24, t35, -t19 * qJDD(1), qJDD(4), -t16 * t19 - t27 * t21, -t16 * t21 + t27 * t19;];
tau_reg = t1;
