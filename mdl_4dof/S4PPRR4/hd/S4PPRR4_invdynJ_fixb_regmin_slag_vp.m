% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:41
% EndTime: 2019-12-31 16:18:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (154->40), mult. (354->69), div. (0->0), fcn. (285->10), ass. (0->37)
t27 = qJD(3) ^ 2;
t15 = pkin(7) + qJ(3);
t13 = sin(t15);
t14 = cos(t15);
t19 = sin(pkin(6));
t21 = cos(pkin(6));
t38 = g(1) * t21 + g(2) * t19;
t18 = sin(pkin(7));
t20 = cos(pkin(7));
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t9 = t18 * t23 - t20 * t25;
t51 = -g(3) * t14 - t9 * qJDD(1) + t38 * t13;
t26 = qJD(4) ^ 2;
t50 = (2 * qJDD(3) * pkin(3)) - pkin(5) * t26 + t51;
t49 = t9 * qJD(1);
t10 = t18 * t25 + t20 * t23;
t47 = g(3) * t13 - t10 * qJDD(1) + t14 * t38;
t22 = sin(qJ(4));
t24 = cos(qJ(4));
t45 = t22 * t24;
t16 = t22 ^ 2;
t44 = -t24 ^ 2 + t16;
t43 = qJD(3) * pkin(3);
t42 = qJD(3) * t9;
t40 = t24 * qJDD(3);
t39 = qJD(3) * qJD(4);
t37 = -qJDD(3) * t9 - t27 * t10;
t36 = g(1) * t19 - g(2) * t21 - qJDD(2);
t33 = t10 * t26 - t37;
t32 = 0.2e1 * t42 * qJD(4) - qJDD(4) * t10;
t3 = t49 - t43;
t30 = -pkin(5) * qJDD(4) + (t3 - t49 - t43) * qJD(4);
t28 = -qJDD(3) * pkin(5) + qJD(1) * t42 - t3 * qJD(3) + t47;
t12 = qJDD(4) * t24 - t22 * t26;
t11 = qJDD(4) * t22 + t24 * t26;
t1 = [qJDD(1) - g(3), -g(3) + (t18 ^ 2 + t20 ^ 2) * qJDD(1), 0, t37, qJD(3) * t42 - qJDD(3) * t10, 0, 0, 0, 0, 0, t22 * t32 - t24 * t33, t22 * t33 + t24 * t32; 0, -t36, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11; 0, 0, qJDD(3), t51, t47, qJDD(3) * t16 + 0.2e1 * t39 * t45, 0.2e1 * t22 * t40 - 0.2e1 * t44 * t39, t11, t12, 0, t30 * t22 + t50 * t24, -t50 * t22 + t30 * t24; 0, 0, 0, 0, 0, -t27 * t45, t44 * t27, t22 * qJDD(3), t40, qJDD(4), t22 * t28 - t24 * t36, t22 * t36 + t24 * t28;];
tau_reg = t1;
