% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S3PPR1
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
%   pkin=[a2,a3,d3]';
% 
% Output:
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S3PPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:01
% EndTime: 2019-03-08 18:02:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (50->22), mult. (88->29), div. (0->0), fcn. (64->2), ass. (0->13)
t8 = qJDD(1) - g(2);
t12 = qJDD(2) - g(1);
t10 = cos(qJ(3));
t9 = sin(qJ(3));
t5 = -t9 * qJD(1) + t10 * qJD(2);
t6 = t10 * qJD(1) + t9 * qJD(2);
t11 = qJD(3) ^ 2;
t7 = t10 * qJDD(2);
t4 = -t9 * qJDD(3) - t10 * t11;
t3 = -t10 * qJDD(3) + t9 * t11;
t2 = -t6 * qJD(3) - t9 * qJDD(1) + t7;
t1 = t5 * qJD(3) + t10 * qJDD(1) + t9 * qJDD(2);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, t4, t3, 0, t1 * t10 - t2 * t9 - g(2) + (-t10 * t5 - t6 * t9) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, -t3, t4, 0, t1 * t9 + t2 * t10 - g(1) + (t10 * t6 - t5 * t9) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -g(1) * t10 - t8 * t9 + t7, -t8 * t10 - t12 * t9, 0, 0;];
tau_reg  = t13;
