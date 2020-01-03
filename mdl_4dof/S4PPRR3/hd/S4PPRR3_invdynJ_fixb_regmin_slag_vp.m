% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:28
% EndTime: 2019-12-31 16:17:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (91->34), mult. (216->65), div. (0->0), fcn. (157->6), ass. (0->28)
t17 = sin(qJ(3));
t15 = cos(pkin(6));
t19 = cos(qJ(3));
t38 = sin(pkin(6));
t3 = -t15 * t19 - t38 * t17;
t4 = t15 * t17 - t38 * t19;
t43 = -g(1) * t3 - g(2) * t4 - t17 * qJDD(2);
t42 = g(1) * t4 - g(2) * t3 + t19 * qJDD(2);
t16 = sin(qJ(4));
t12 = t16 ^ 2;
t18 = cos(qJ(4));
t41 = -t18 ^ 2 + t12;
t20 = qJD(4) ^ 2;
t21 = qJD(3) ^ 2;
t40 = t20 + t21;
t39 = qJD(3) * pkin(3);
t14 = qJDD(1) - g(3);
t35 = qJDD(4) * t16;
t33 = t18 * qJDD(3);
t31 = t19 * qJDD(3);
t29 = qJD(3) * qJD(4);
t28 = t16 * t29;
t24 = -qJDD(3) * pkin(5) + qJD(3) * t39 + t43;
t23 = -pkin(5) * qJDD(4) - 0.2e1 * t39 * qJD(4);
t22 = 0.2e1 * qJDD(3) * pkin(3) - pkin(5) * t20 + t42;
t6 = qJDD(4) * t18 - t20 * t16;
t5 = t20 * t18 + t35;
t1 = [t14, t14, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5; 0, -g(1) * t38 + g(2) * t15 + qJDD(2), 0, -t21 * t17 + t31, -qJDD(3) * t17 - t21 * t19, 0, 0, 0, 0, 0, (-0.2e1 * t28 + t33) * t19 + (-t40 * t18 - t35) * t17, (-qJDD(4) * t17 - 0.2e1 * t19 * t29) * t18 + (t40 * t17 - t31) * t16; 0, 0, qJDD(3), t42, t43, t12 * qJDD(3) + 0.2e1 * t18 * t28, 0.2e1 * t16 * t33 - 0.2e1 * t41 * t29, t5, t6, 0, t23 * t16 + t22 * t18, -t22 * t16 + t23 * t18; 0, 0, 0, 0, 0, -t16 * t21 * t18, t41 * t21, t16 * qJDD(3), t33, qJDD(4), -t14 * t18 + t24 * t16, t14 * t16 + t24 * t18;];
tau_reg = t1;
