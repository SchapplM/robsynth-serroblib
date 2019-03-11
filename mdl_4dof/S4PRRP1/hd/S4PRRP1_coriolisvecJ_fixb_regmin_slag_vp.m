% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauc_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:59
% EndTime: 2019-03-08 18:23:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->23), mult. (107->32), div. (0->0), fcn. (33->2), ass. (0->20)
t11 = cos(qJ(3));
t18 = pkin(2) * qJD(3);
t13 = qJD(2) * t18;
t7 = t11 * t13;
t9 = qJD(2) + qJD(3);
t8 = t9 * qJD(4);
t3 = t7 + t8;
t19 = pkin(2) * qJD(2);
t10 = sin(qJ(3));
t17 = t10 * t18;
t16 = t11 * t18;
t15 = t10 * t19;
t14 = t11 * t19;
t12 = t9 * t14 - t7;
t6 = qJD(4) + t16;
t5 = t9 * qJ(4) + t15;
t4 = -t9 * pkin(3) + qJD(4) - t14;
t2 = (-qJD(2) - t9) * t17;
t1 = (-qJD(3) + t9) * t15;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t2, -t9 * t16 - t7, t2, t6 * t9 + t3, t3 * (t10 * pkin(2) + qJ(4)) + t5 * t6 + (t4 + (-t11 * pkin(2) - pkin(3)) * qJD(2)) * t17; 0, 0, 0, 0, 0, t1, t12, t1, -t12 + 0.2e1 * t8, t3 * qJ(4) + t5 * qJD(4) + (-t11 * t5 + (-pkin(3) * qJD(3) - t4) * t10) * t19; 0, 0, 0, 0, 0, 0, 0, 0, -t9 ^ 2, t10 * t13 - t5 * t9;];
tauc_reg  = t20;
