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
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (30->20), mult. (101->40), div. (0->0), fcn. (58->4), ass. (0->23)
t28 = pkin(1) * qJDD(1);
t6 = sin(qJ(1));
t8 = cos(qJ(1));
t27 = g(1) * t6 + g(3) * t8;
t5 = sin(qJ(2));
t3 = t5 ^ 2;
t7 = cos(qJ(2));
t4 = t7 ^ 2;
t26 = (t3 + t4) * t28;
t25 = t5 * t7;
t23 = t3 - t4;
t22 = qJDD(2) * t5;
t21 = qJDD(2) * t7;
t20 = t7 * qJDD(1);
t19 = qJD(1) * qJD(2);
t10 = qJD(1) ^ 2;
t18 = t10 * t25;
t16 = t19 * t25;
t15 = -g(1) * t8 + g(3) * t6;
t9 = qJD(2) ^ 2;
t13 = pkin(1) * t9 + t15;
t12 = t27 + t28;
t1 = [0, 0, 0, 0, 0, qJDD(1), t15, t27, 0, 0, t3 * qJDD(1) + 0.2e1 * t16, -0.2e1 * t23 * t19 + 0.2e1 * t5 * t20, -t9 * t7 - t22, t4 * qJDD(1) - 0.2e1 * t16, t9 * t5 - t21, 0, pkin(1) * t22 + t13 * t7, pkin(1) * t21 - t13 * t5, t27 + 0.2e1 * t26, (t27 + t26) * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t23 * t10, -t5 * qJDD(1), t18, -t20, qJDD(2), g(2) * t7 + t12 * t5, -g(2) * t5 + t12 * t7, 0, 0;];
tau_reg = t1;
