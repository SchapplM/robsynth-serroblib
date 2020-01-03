% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (53->31), mult. (160->49), div. (0->0), fcn. (77->4), ass. (0->28)
t7 = sin(qJ(2));
t24 = t7 * qJD(1);
t22 = qJD(2) * qJ(3);
t3 = t22 + t24;
t18 = -t3 + t24;
t31 = qJD(4) * (-t18 + t22);
t15 = t18 * qJD(2);
t9 = cos(qJ(2));
t30 = t3 * t9;
t6 = sin(qJ(4));
t8 = cos(qJ(4));
t29 = t6 ^ 2 - t8 ^ 2;
t12 = qJD(2) ^ 2;
t28 = t12 * t7;
t27 = t12 * t9;
t11 = qJD(4) ^ 2;
t26 = t11 + t12;
t25 = qJD(2) * pkin(2);
t23 = t9 * qJD(1);
t21 = qJD(2) * qJD(4);
t20 = 0.2e1 * t21;
t19 = t26 * t9;
t17 = -0.2e1 * t6 * t21;
t16 = qJD(3) - t23;
t1 = (qJD(3) + t23) * qJD(2);
t13 = t16 * qJD(2) - (-pkin(2) - pkin(5)) * t11 + t1;
t2 = t16 - t25;
t4 = [0, 0, -t28, -t27, t28, t27, t1 * t7 + (t30 + (t2 - t23) * t7) * qJD(2), 0, 0, 0, 0, 0, t8 * t7 * t20 + t6 * t19, t7 * t17 + t8 * t19; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t1 * qJ(3) + t3 * qJD(3) + (-t30 + (-t2 - t25) * t7) * qJD(1), t8 * t17, t29 * t20, -t11 * t6, -t11 * t8, 0, t13 * t6 + t8 * t31, t13 * t8 - t6 * t31; 0, 0, 0, 0, 0, -t12, t15, 0, 0, 0, 0, 0, -t26 * t6, -t26 * t8; 0, 0, 0, 0, 0, 0, 0, t8 * t12 * t6, -t29 * t12, 0, 0, 0, t8 * t15, -t6 * t15;];
tauc_reg = t4;
