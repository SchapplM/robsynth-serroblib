% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
% 
% Output:
% tau_reg [4x6]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:26
% EndTime: 2019-03-08 18:11:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->21), mult. (90->25), div. (0->0), fcn. (84->6), ass. (0->16)
t18 = qJD(4) ^ 2;
t12 = sin(pkin(5));
t13 = cos(pkin(5));
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t3 = -t15 * t12 - t14 * t13;
t17 = -t14 * t12 + t15 * t13;
t16 = t3 * qJDD(4) - t17 * t18;
t11 = qJDD(1) - g(2);
t10 = pkin(5) + qJ(4);
t9 = cos(t10);
t8 = sin(t10);
t5 = t13 * qJDD(1) + t12 * qJDD(2);
t4 = -t12 * qJDD(1) + t13 * qJDD(2);
t1 = t17 * qJDD(4) + t3 * t18;
t2 = [t11, t11, -t4 * t12 + t5 * t13 - g(2), 0, t16, -t1; 0, qJDD(2) - g(1), t5 * t12 + t4 * t13 - g(1), 0, t1, t16; 0, 0, qJDD(3) + g(3), 0, 0, 0; 0, 0, 0, qJDD(4), -g(1) * t9 + g(2) * t8 - t14 * t5 + t15 * t4, g(1) * t8 + g(2) * t9 - t14 * t4 - t15 * t5;];
tau_reg  = t2;
