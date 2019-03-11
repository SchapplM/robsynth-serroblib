% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:49
% EndTime: 2019-03-08 18:21:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (117->43), mult. (296->71), div. (0->0), fcn. (228->6), ass. (0->28)
t17 = qJD(2) + qJD(4);
t29 = qJD(4) - t17;
t18 = sin(pkin(6));
t28 = pkin(2) * t18;
t21 = sin(qJ(2));
t27 = qJD(1) * t21;
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t19 = cos(pkin(6));
t23 = cos(qJ(2));
t13 = t18 * t23 + t19 * t21;
t8 = t13 * qJD(2);
t6 = qJD(1) * t8;
t12 = -t18 * t21 + t19 * t23;
t10 = t12 * qJD(2);
t7 = qJD(1) * t10;
t26 = -t20 * t7 - t22 * t6;
t15 = qJD(2) * pkin(2) + t23 * qJD(1);
t4 = t19 * t15 - t18 * t27;
t2 = qJD(2) * pkin(3) + t4;
t5 = t18 * t15 + t19 * t27;
t25 = -t20 * t2 - t22 * t5;
t24 = qJD(2) ^ 2;
t16 = t19 * pkin(2) + pkin(3);
t11 = t12 * qJD(1);
t9 = t13 * qJD(1);
t1 = qJD(4) * t20 * t5;
t3 = [0, 0, -t24 * t21, -t24 * t23, t5 * t10 - t6 * t12 + t7 * t13 - t4 * t8, 0 (-t20 * t10 - t22 * t8 + (-t12 * t20 - t13 * t22) * qJD(4)) * t17 -(t22 * t10 - t20 * t8 + (t12 * t22 - t13 * t20) * qJD(4)) * t17; 0, 0, 0, 0, -t5 * t11 + t4 * t9 + (t18 * t7 - t19 * t6) * pkin(2), 0 -(-t20 * t11 - t22 * t9) * t17 + ((-t16 * t20 - t22 * t28) * t17 + t25) * qJD(4) + t26, t1 - t22 * t7 + t20 * t6 + (t22 * t11 - t20 * t9) * t17 + (-(t16 * t22 - t20 * t28) * t17 - t22 * t2) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t29 * t25 + t26, t1 + (-t5 * t17 + t6) * t20 + (-t29 * t2 - t7) * t22;];
tauc_reg  = t3;
