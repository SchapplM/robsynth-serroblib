% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [2x(2*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S2RR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:10
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (54->31), mult. (130->45), div. (0->0), fcn. (84->4), ass. (0->26)
t12 = sin(qJ(2));
t10 = t12 ^ 2;
t14 = cos(qJ(2));
t11 = t14 ^ 2;
t29 = t10 - t11;
t17 = qJD(1) ^ 2;
t7 = t14 * t17 * t12;
t28 = t12 * (qJDD(2) + t7);
t27 = t14 * (qJDD(2) - t7);
t26 = t10 * t17;
t25 = t11 * t17;
t24 = qJDD(1) * pkin(1);
t23 = t12 * qJDD(1);
t22 = t14 * qJDD(1);
t21 = qJD(1) * qJD(2);
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t20 = -g(1) * t15 + t13 * g(3);
t19 = g(1) * t13 + g(3) * t15;
t3 = -t19 - t24;
t1 = -t14 * g(2) + t12 * t3;
t2 = g(2) * t12 + t14 * t3;
t18 = t12 * t1 + t14 * t2;
t16 = qJD(2) ^ 2;
t4 = pkin(1) * t17 - t20;
t5 = [0, 0, 0, 0, 0, qJDD(1), t20, t19, 0, 0, (0.2e1 * t14 * t21 + t23) * t12, 0.2e1 * t12 * t22 - 0.2e1 * t29 * t21, -t28 - t14 * (t16 - t26), (-0.2e1 * t12 * t21 + t22) * t14, -t12 * (-t16 + t25) - t27, 0, -t14 * t4 - pkin(1) * (t14 * (-t16 - t25) - t28), t12 * t4 - pkin(1) * (-t27 - t12 * (-t16 - t26)), -(-t10 - t11) * t24 - t18, -pkin(1) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t29 * t17, -t23, t7, -t22, qJDD(2), -t1, -t2, 0, 0;];
tauJ_reg = t5;
