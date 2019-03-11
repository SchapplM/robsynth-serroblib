% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S3PRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:06
% EndTime: 2019-03-08 18:04:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (161->54), mult. (316->72), div. (0->0), fcn. (230->6), ass. (0->37)
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t20 = cos(t23);
t41 = g(1) * t20 + g(2) * t19;
t40 = g(1) * t19 - g(2) * t20;
t22 = qJD(2) + qJD(3);
t27 = cos(qJ(2));
t18 = t27 * qJDD(1);
t25 = sin(qJ(2));
t35 = qJD(1) * qJD(2);
t10 = qJDD(2) * pkin(2) - t25 * t35 + t18;
t24 = sin(qJ(3));
t37 = qJD(1) * t25;
t34 = qJD(3) * t37;
t13 = t24 * t34;
t26 = cos(qJ(3));
t14 = qJD(2) * pkin(2) + t27 * qJD(1);
t31 = -t25 * qJDD(1) - t27 * t35;
t30 = qJD(3) * t14 - t31;
t1 = t24 * t10 + t30 * t26 - t13;
t21 = qJDD(2) + qJDD(3);
t39 = pkin(2) * t21;
t36 = qJDD(1) - g(2);
t32 = g(1) * t25 - g(2) * t27;
t12 = t24 * t27 + t26 * t25;
t11 = -t24 * t25 + t26 * t27;
t29 = (-pkin(2) * t22 - t14) * qJD(3) + t31;
t7 = t26 * t10;
t2 = -t30 * t24 - t26 * t34 + t7;
t28 = qJD(2) ^ 2;
t9 = t11 * qJD(1);
t8 = t12 * qJD(1);
t6 = t24 * t14 + t26 * t37;
t5 = t26 * t14 - t24 * t37;
t4 = t22 * t12;
t3 = t22 * t11;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t27 * qJDD(2) - t28 * t25, -qJDD(2) * t25 - t28 * t27, 0, -g(2) + (t25 ^ 2 + t27 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t11 * t21 - t4 * t22, -t12 * t21 - t3 * t22, 0, t1 * t12 + t2 * t11 + t6 * t3 - t5 * t4 - g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t18 + t32, g(1) * t27 - t36 * t25, 0, 0, 0, 0, 0, 0, 0, t21, t8 * t22 + t7 + (-t34 + t39) * t26 + t29 * t24 + t40, t9 * t22 + t13 + (-t10 - t39) * t24 + t29 * t26 + t41, 0, t5 * t8 - t6 * t9 + (t1 * t24 + t2 * t26 + (-t24 * t5 + t26 * t6) * qJD(3) + t32) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t6 * t22 + t2 + t40, t5 * t22 - t1 + t41, 0, 0;];
tau_reg  = t15;
