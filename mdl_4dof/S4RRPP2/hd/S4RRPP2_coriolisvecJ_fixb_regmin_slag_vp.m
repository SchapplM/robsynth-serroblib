% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauc_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:05
% EndTime: 2019-03-08 18:34:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->33), mult. (192->39), div. (0->0), fcn. (58->2), ass. (0->30)
t15 = qJD(1) + qJD(2);
t33 = t15 ^ 2;
t32 = -pkin(2) - pkin(3);
t16 = sin(qJ(2));
t28 = pkin(1) * qJD(1);
t24 = t16 * t28;
t11 = t15 * qJ(3) + t24;
t17 = cos(qJ(2));
t27 = pkin(1) * qJD(2);
t25 = t17 * t27;
t12 = qJD(3) + t25;
t21 = qJD(1) * t27;
t13 = t17 * t21;
t14 = t15 * qJD(3);
t9 = t13 + t14;
t31 = t9 * (t16 * pkin(1) + qJ(3)) + t11 * t12;
t30 = t9 * qJ(3) + t11 * qJD(3);
t29 = t11 * t17;
t26 = t16 * t27;
t23 = t17 * t28;
t22 = -t17 * pkin(1) - pkin(2);
t20 = t12 * t15 + t9;
t19 = qJD(3) - t23;
t18 = t15 * t23 - t13;
t6 = (-qJD(1) - t15) * t26;
t5 = (-qJD(2) + t15) * t24;
t10 = -t15 * pkin(2) + t19;
t4 = t32 * t15 + t19;
t1 = -t11 * t15 + t16 * t21;
t2 = [0, 0, 0, 0, t6, -t15 * t25 - t13, t6, t20 (t22 * qJD(1) + t10) * t26 + t31, t6, t20 (t4 + (-pkin(3) + t22) * qJD(1)) * t26 + t31; 0, 0, 0, 0, t5, t18, t5, 0.2e1 * t14 - t18 (-t29 + (-pkin(2) * qJD(2) - t10) * t16) * t28 + t30, t5, t19 * t15 + t9 (-t29 + (qJD(2) * t32 - t4) * t16) * t28 + t30; 0, 0, 0, 0, 0, 0, 0, -t33, t1, 0, -t33, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t2;
