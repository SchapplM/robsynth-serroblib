% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S3PRR1
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
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% tau_reg [3x7]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3PRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:58
% EndTime: 2018-11-14 10:11:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (89->40), mult. (144->54), div. (0->0), fcn. (106->6), ass. (0->31)
t12 = qJD(2) + qJD(3);
t35 = t12 ^ 2;
t11 = qJDD(2) + qJDD(3);
t34 = pkin(2) * t11;
t14 = sin(qJ(3));
t17 = cos(qJ(2));
t33 = t14 * t17;
t15 = sin(qJ(2));
t32 = qJD(1) * t15;
t31 = qJD(3) * t15;
t30 = t17 * qJD(1);
t29 = -qJD(3) + t12;
t28 = qJDD(1) - g(2);
t27 = t15 * qJDD(1);
t13 = qJ(2) + qJ(3);
t10 = cos(t13);
t9 = sin(t13);
t26 = t14 * qJD(1) * t31 + g(1) * t10 + g(2) * t9;
t25 = qJD(1) * qJD(2);
t16 = cos(qJ(3));
t8 = t17 * qJDD(1);
t2 = qJDD(2) * pkin(2) - t15 * t25 + t8;
t24 = g(1) * t9 - g(2) * t10 + t16 * t2;
t23 = -t12 * t32 - t2;
t22 = t16 * t15 + t33;
t21 = -t14 * t15 + t16 * t17;
t4 = qJD(2) * pkin(2) + t30;
t20 = -t27 + (-pkin(2) * t12 - t4) * qJD(3);
t19 = -t17 * t25 + t29 * t4 - t27;
t18 = qJD(2) ^ 2;
t1 = [t28, 0, t17 * qJDD(2) - t18 * t15, -qJDD(2) * t15 - t18 * t17, 0, t21 * t11 - t35 * t22, -t22 * t11 - t35 * t21; 0, qJDD(2), g(1) * t15 - g(2) * t17 + t8, g(1) * t17 - t28 * t15, t11, t16 * t34 + t20 * t14 + (-qJD(2) * t33 + t22 * t12 - t16 * t31) * qJD(1) + t24 (t23 - t34) * t14 + ((-qJD(2) + t12) * t30 + t20) * t16 + t26; 0, 0, 0, 0, t11, t29 * t16 * t32 + t19 * t14 + t24, t23 * t14 + t19 * t16 + t26;];
tau_reg  = t1;
