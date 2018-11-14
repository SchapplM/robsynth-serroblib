% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S3RPR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% tau_reg [3x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3RPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:29
% EndTime: 2018-11-14 10:14:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (109->44), mult. (153->47), div. (0->0), fcn. (86->4), ass. (0->29)
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t42 = g(1) * t12 - g(2) * t14;
t41 = -qJDD(2) + t42;
t35 = pkin(1) * qJDD(1);
t40 = t35 + t41;
t15 = -pkin(1) - pkin(2);
t3 = t15 * qJDD(1) + qJDD(2);
t28 = qJD(1) - qJD(3);
t33 = qJD(3) + t28;
t39 = t33 * qJ(2) * qJD(1) - t3;
t29 = qJD(1) * qJD(2);
t4 = t15 * qJD(1) + qJD(2);
t9 = qJDD(1) - qJDD(3);
t38 = (qJD(3) * t15 + qJD(2)) * t28 + t29 + qJD(3) * t4 + (t9 + qJDD(1)) * qJ(2);
t37 = (qJD(1) + t28) * qJ(2) * qJD(3) - t15 * t9 - t3;
t30 = qJ(2) * qJDD(1);
t25 = t28 ^ 2;
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t1 = -t12 * t11 - t14 * t13;
t2 = t14 * t11 - t12 * t13;
t23 = -g(1) * t2 + g(2) * t1;
t22 = -g(1) * t1 - g(2) * t2;
t21 = g(1) * t14 + g(2) * t12;
t19 = -t21 + 0.2e1 * t29;
t17 = -t33 * t4 - t29 - t30;
t16 = qJD(1) ^ 2;
t5 = [qJDD(1), t42, t21, 0.2e1 * t35 + t41, t19 + 0.2e1 * t30, t40 * pkin(1) + (t19 + t30) * qJ(2), t9, t11 * t38 + t13 * t37 + t23, -t11 * t37 + t13 * t38 - t22; 0, 0, 0, -qJDD(1), -t16, -t16 * qJ(2) - t40, 0, -t11 * t25 - t13 * t9, t11 * t9 - t13 * t25; 0, 0, 0, 0, 0, 0, -t9, t17 * t11 - t13 * t39 - t23, t11 * t39 + t17 * t13 + t22;];
tau_reg  = t5;
