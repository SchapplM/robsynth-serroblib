% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [2x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S2RR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:27
% EndTime: 2018-11-16 16:44:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (23->16), mult. (70->35), div. (0->0), fcn. (45->4), ass. (0->18)
t3 = sin(qJ(2));
t5 = cos(qJ(2));
t18 = t3 * t5;
t1 = t3 ^ 2;
t17 = -t5 ^ 2 + t1;
t16 = qJDD(2) * t3;
t15 = qJDD(2) * t5;
t14 = t5 * qJDD(1);
t13 = qJD(1) * qJD(2);
t4 = sin(qJ(1));
t6 = cos(qJ(1));
t12 = g(1) * t6 - g(3) * t4;
t11 = -g(1) * t4 - g(3) * t6;
t7 = qJD(2) ^ 2;
t10 = pkin(1) * t7 + t12;
t9 = pkin(1) * qJDD(1) + t11;
t8 = qJD(1) ^ 2;
t2 = [qJDD(1), t12, t11, t1 * qJDD(1) + 0.2e1 * t13 * t18, -0.2e1 * t17 * t13 + 0.2e1 * t3 * t14, -t7 * t5 - t16, t7 * t3 - t15, 0, pkin(1) * t16 + t10 * t5, pkin(1) * t15 - t10 * t3; 0, 0, 0, -t8 * t18, t17 * t8, -t3 * qJDD(1), -t14, qJDD(2), g(2) * t5 + t9 * t3, -g(2) * t3 + t9 * t5;];
tau_reg  = t2;
