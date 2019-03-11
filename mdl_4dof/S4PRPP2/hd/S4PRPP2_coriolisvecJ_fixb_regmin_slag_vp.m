% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% tauc_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:59
% EndTime: 2019-03-08 18:19:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (80->37), mult. (213->49), div. (0->0), fcn. (154->4), ass. (0->27)
t18 = sin(pkin(5));
t19 = cos(pkin(5));
t20 = sin(qJ(2));
t21 = cos(qJ(2));
t13 = t18 * t21 + t19 * t20;
t27 = t13 * qJD(1);
t12 = t18 * t20 - t19 * t21;
t8 = t13 * qJD(2);
t6 = qJD(1) * t8;
t26 = t6 * t12;
t24 = t21 * qJD(1);
t16 = qJD(2) * pkin(2) + t24;
t25 = qJD(1) * t20;
t17 = t19 * t25;
t5 = t18 * t16 + t17;
t23 = t18 * t25;
t4 = t19 * t16 - t23;
t22 = qJD(2) ^ 2;
t15 = t19 * qJD(2) * t24;
t11 = t12 * qJD(1);
t10 = t12 * qJD(2);
t9 = t18 * t24 + t17;
t7 = -qJD(2) * t23 + t15;
t3 = t15 + (qJD(4) - t23) * qJD(2);
t2 = qJD(2) * qJ(4) + t5;
t1 = -qJD(2) * pkin(3) + qJD(4) - t4;
t14 = [0, 0, -t22 * t20, -t22 * t21, -t5 * t10 + t7 * t13 - t4 * t8 + t26, -t8 * qJD(2), -t10 * qJD(2), t1 * t8 - t2 * t10 + t3 * t13 + t26; 0, 0, 0, 0, t5 * t11 + t4 * t9 + (t18 * t7 - t6 * t19) * pkin(2) (t9 - t27) * qJD(2), t15 + (0.2e1 * qJD(4) + t11 - t23) * qJD(2), t3 * (t18 * pkin(2) + qJ(4)) + t6 * (-t19 * pkin(2) - pkin(3)) - t1 * t9 + (qJD(4) + t11) * t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t22 (-t2 + t27) * qJD(2);];
tauc_reg  = t14;
