% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% tauc_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S3PRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:00
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (100->30), mult. (250->48), div. (0->0), fcn. (176->4), ass. (0->28)
t12 = qJD(2) + qJD(3);
t15 = cos(qJ(3));
t16 = cos(qJ(2));
t26 = t16 * qJD(1);
t24 = qJD(2) * t26;
t11 = qJD(2) * pkin(2) + t26;
t27 = qJD(3) * t11;
t13 = sin(qJ(3));
t14 = sin(qJ(2));
t28 = qJD(1) * t14;
t25 = t13 * t28;
t30 = t12 * t25;
t1 = (t24 + t27) * t15 - t30;
t29 = t15 * t14;
t23 = t13 * t16 + t29;
t22 = -t13 * t14 + t15 * t16;
t21 = (-pkin(2) * t12 - t11) * qJD(3);
t20 = t23 * qJD(2);
t18 = (-qJD(3) * t29 - t20) * qJD(1);
t2 = -t13 * t27 + t18;
t17 = qJD(2) ^ 2;
t8 = t22 * qJD(1);
t7 = t23 * qJD(1);
t6 = t13 * t11 + t15 * t28;
t5 = t15 * t11 - t25;
t4 = -t23 * qJD(3) - t20;
t3 = t12 * t22;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t14, -t17 * t16, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t12, -t3 * t12, 0, t1 * t23 + t2 * t22 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t12 + t13 * t21 + t18, t8 * t12 + (t21 - t24) * t15 + t30, 0, t5 * t7 - t6 * t8 + (t1 * t13 + t2 * t15 + (-t13 * t5 + t15 * t6) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t12 + t2, t5 * t12 - t1, 0, 0;];
tauc_reg  = t9;
