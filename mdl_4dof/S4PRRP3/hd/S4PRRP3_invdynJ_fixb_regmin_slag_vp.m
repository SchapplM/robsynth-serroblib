% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:57
% EndTime: 2019-12-31 16:26:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (215->71), mult. (445->105), div. (0->0), fcn. (234->4), ass. (0->44)
t44 = qJ(4) + pkin(5);
t47 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t44 * qJDD(2);
t23 = sin(qJ(3));
t11 = t44 * t23;
t24 = cos(qJ(3));
t4 = t24 * qJD(1) - qJD(2) * t11;
t42 = qJD(3) * pkin(3);
t3 = t4 + t42;
t46 = t3 - t4;
t45 = g(3) * t24;
t20 = t23 ^ 2;
t21 = t24 ^ 2;
t43 = t20 - t21;
t41 = qJDD(1) - g(3);
t40 = t23 * qJDD(2);
t39 = t24 * qJDD(2);
t38 = qJD(2) * qJD(3);
t14 = t24 * pkin(3) + pkin(2);
t12 = t44 * t24;
t37 = t23 * t38;
t36 = t44 * qJD(3);
t19 = pkin(6) + qJ(2);
t15 = sin(t19);
t16 = cos(t19);
t35 = g(1) * t16 + g(2) * t15;
t34 = g(1) * t15 - g(2) * t16;
t6 = t23 * qJD(1) + qJD(2) * t12;
t33 = t23 * t3 - t24 * t6;
t32 = qJD(2) * t36;
t30 = -0.2e1 * pkin(2) * t38 - pkin(5) * qJDD(3);
t29 = pkin(3) * t37 - qJDD(2) * t14 + qJDD(4);
t25 = qJD(3) ^ 2;
t28 = 0.2e1 * qJDD(2) * pkin(2) - pkin(5) * t25 + t34;
t26 = qJD(2) ^ 2;
t27 = t26 * pkin(2) - qJDD(2) * pkin(5) + t35;
t17 = t24 * qJDD(1);
t10 = qJDD(3) * t24 - t25 * t23;
t9 = qJDD(3) * t23 + t25 * t24;
t8 = -qJD(2) * t14 + qJD(4);
t7 = -t23 * qJD(4) - t24 * t36;
t5 = t24 * qJD(4) - t23 * t36;
t2 = (qJDD(1) - t32) * t23 + t47 * t24;
t1 = qJDD(3) * pkin(3) - t47 * t23 - t24 * t32 + t17;
t13 = [t41, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, -qJD(3) * t33 + t1 * t24 + t2 * t23 - g(3); 0, qJDD(2), t34, t35, t20 * qJDD(2) + 0.2e1 * t24 * t37, 0.2e1 * t23 * t39 - 0.2e1 * t43 * t38, t9, t10, 0, t23 * t30 + t24 * t28, -t23 * t28 + t24 * t30, (-qJD(3) * t3 + qJDD(2) * t12 + t2 + (qJD(3) * t11 + t5) * qJD(2)) * t24 + (-qJD(3) * t6 + qJDD(2) * t11 - t1 + (-qJD(3) * t12 - t7) * qJD(2)) * t23 - t35, t2 * t12 + t6 * t5 - t1 * t11 + t3 * t7 - t29 * t14 + t8 * t23 * t42 - g(1) * (-t15 * t14 + t16 * t44) - g(2) * (t16 * t14 + t15 * t44); 0, 0, 0, 0, -t23 * t26 * t24, t43 * t26, t40, t39, qJDD(3), t23 * t27 + t17 - t45, -t41 * t23 + t27 * t24, -pkin(3) * t40 + (-t42 + t46) * t24 * qJD(2), t46 * t6 + (-t45 + t1 + (-qJD(2) * t8 + t35) * t23) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t20 - t21) * t26, qJD(2) * t33 + t29 - t34;];
tau_reg = t13;
