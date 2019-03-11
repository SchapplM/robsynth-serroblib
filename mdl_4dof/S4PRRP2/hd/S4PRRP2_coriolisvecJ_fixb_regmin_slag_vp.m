% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% tauc_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (114->34), mult. (268->51), div. (0->0), fcn. (184->4), ass. (0->30)
t14 = qJD(2) + qJD(3);
t17 = cos(qJ(3));
t18 = cos(qJ(2));
t28 = t18 * qJD(1);
t26 = qJD(2) * t28;
t13 = qJD(2) * pkin(2) + t28;
t29 = qJD(3) * t13;
t15 = sin(qJ(3));
t16 = sin(qJ(2));
t30 = qJD(1) * t16;
t27 = t15 * t30;
t31 = t14 * t27;
t1 = (t26 + t29) * t17 - t31;
t32 = t17 * t16;
t25 = t15 * t18 + t32;
t22 = t25 * qJD(2);
t20 = (-qJD(3) * t32 - t22) * qJD(1);
t2 = -t15 * t29 + t20;
t33 = t2 * pkin(3);
t6 = t17 * t13 - t27;
t24 = -t15 * t16 + t17 * t18;
t23 = (-pkin(2) * t14 - t13) * qJD(3);
t19 = qJD(2) ^ 2;
t9 = t24 * qJD(1);
t8 = t25 * qJD(1);
t7 = t15 * t13 + t17 * t30;
t5 = t14 * pkin(3) + t6;
t4 = -t25 * qJD(3) - t22;
t3 = t14 * t24;
t10 = [0, 0, -t19 * t16, -t19 * t18, 0, t4 * t14, -t3 * t14, t1 * t25 + t2 * t24 + t7 * t3 + t5 * t4; 0, 0, 0, 0, 0, t8 * t14 + t15 * t23 + t20, t9 * t14 + (t23 - t26) * t17 + t31, t33 + t5 * t8 - t7 * t9 + (t1 * t15 + t2 * t17 + (-t15 * t5 + t17 * t7) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, t7 * t14 + t2, t6 * t14 - t1, t33 + (t5 - t6) * t7; 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t10;
