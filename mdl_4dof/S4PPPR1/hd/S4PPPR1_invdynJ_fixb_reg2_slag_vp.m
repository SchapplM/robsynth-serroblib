% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPPR1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:09:05
% EndTime: 2019-03-08 18:09:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->28), mult. (106->37), div. (0->0), fcn. (86->4), ass. (0->20)
t13 = sin(pkin(5));
t14 = cos(pkin(5));
t21 = -g(1) * t13 + g(2) * t14;
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t20 = -t15 * qJDD(2) + t16 * qJDD(3);
t19 = -g(1) * t14 - g(2) * t13;
t9 = t16 * qJD(2) + t15 * qJD(3);
t8 = -t15 * qJD(2) + t16 * qJD(3);
t18 = t16 * qJDD(2) + t15 * qJDD(3);
t17 = qJD(4) ^ 2;
t12 = qJDD(1) - g(3);
t7 = -t15 * qJDD(4) - t16 * t17;
t6 = -t16 * qJDD(4) + t15 * t17;
t5 = qJDD(2) + t21;
t4 = t13 * t16 + t14 * t15;
t3 = -t13 * t15 + t14 * t16;
t2 = -t9 * qJD(4) + t20;
t1 = t8 * qJD(4) + t18;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, t7, t6, 0, t1 * t16 - t2 * t15 + (-t15 * t9 - t16 * t8) * qJD(4) + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t19, 0, 0, 0, 0, 0, 0, -t6, t7, 0, t1 * t15 + t2 * t16 + (-t15 * t8 + t16 * t9) * qJD(4) + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), -g(1) * t3 - g(2) * t4 + t20, g(1) * t4 - g(2) * t3 - t18, 0, 0;];
tau_reg  = t10;
