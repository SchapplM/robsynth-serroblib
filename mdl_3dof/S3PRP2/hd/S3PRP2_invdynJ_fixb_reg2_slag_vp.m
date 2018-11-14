% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S3PRP2
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
%   pkin=[a2,a3,d2]';
% 
% Output:
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3PRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:07:18
% EndTime: 2018-11-14 10:07:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->38), mult. (91->42), div. (0->0), fcn. (50->2), ass. (0->19)
t11 = sin(qJ(2));
t12 = cos(qJ(2));
t22 = -g(1) * t12 + g(2) * t11;
t20 = qJDD(2) * pkin(2);
t19 = t11 * qJD(1);
t18 = t12 * qJD(1);
t17 = qJDD(2) * qJ(3);
t9 = t12 * qJDD(1);
t16 = t9 + t22;
t6 = qJD(2) * qJ(3) + t19;
t15 = (-qJD(2) * pkin(2) + qJD(3) - t18) * t11 + t6 * t12;
t8 = t11 * qJDD(1);
t14 = g(1) * t11 + g(2) * t12 - t8;
t2 = qJD(2) * t19 + qJDD(3) - t20 - t9;
t13 = qJD(2) ^ 2;
t4 = t12 * qJDD(2) - t13 * t11;
t3 = qJDD(2) * t11 + t13 * t12;
t1 = t17 + t8 + (qJD(3) + t18) * qJD(2);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(1), 0, 0, 0, 0, 0, 0, t4, -t3, 0, -g(1) + (t11 ^ 2 + t12 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t4, 0, t3, t15 * qJD(2) + t1 * t11 - t2 * t12 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t16, t14, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -qJDD(3) + t16 + 0.2e1 * t20, 0, 0.2e1 * qJD(2) * qJD(3) - t14 + 0.2e1 * t17, t1 * qJ(3) + t6 * qJD(3) - t2 * pkin(2) - g(1) * (t12 * pkin(2) + t11 * qJ(3)) - g(2) * (-t11 * pkin(2) + t12 * qJ(3)) - t15 * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t13, -t6 * qJD(2) + t2 - t22;];
tau_reg  = t5;
