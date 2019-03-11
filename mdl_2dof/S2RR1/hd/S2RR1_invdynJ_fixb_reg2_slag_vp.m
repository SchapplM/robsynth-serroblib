% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% tau_reg [2x(2*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S2RR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:03
% EndTime: 2019-03-08 18:00:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (30->20), mult. (101->40), div. (0->0), fcn. (58->4), ass. (0->23)
t24 = pkin(1) * qJDD(1);
t3 = sin(qJ(2));
t1 = t3 ^ 2;
t5 = cos(qJ(2));
t2 = t5 ^ 2;
t23 = (t1 + t2) * t24;
t22 = t3 * t5;
t21 = t1 - t2;
t20 = qJDD(2) * t3;
t19 = qJDD(2) * t5;
t18 = t5 * qJDD(1);
t8 = qJD(1) ^ 2;
t17 = t8 * t22;
t16 = qJD(1) * qJD(2);
t14 = t16 * t22;
t4 = sin(qJ(1));
t6 = cos(qJ(1));
t13 = g(1) * t6 - g(3) * t4;
t12 = -g(1) * t4 - g(3) * t6;
t7 = qJD(2) ^ 2;
t11 = pkin(1) * t7 + t13;
t10 = t12 + t24;
t9 = [0, 0, 0, 0, 0, qJDD(1), t13, t12, 0, 0, t1 * qJDD(1) + 0.2e1 * t14, -0.2e1 * t21 * t16 + 0.2e1 * t3 * t18, -t7 * t5 - t20, t2 * qJDD(1) - 0.2e1 * t14, t7 * t3 - t19, 0, pkin(1) * t20 + t11 * t5, pkin(1) * t19 - t11 * t3, t12 + 0.2e1 * t23 (t12 + t23) * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t21 * t8, -t3 * qJDD(1), t17, -t18, qJDD(2), g(2) * t5 + t10 * t3, -g(2) * t3 + t10 * t5, 0, 0;];
tau_reg  = t9;
