% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP4
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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (143->51), mult. (378->84), div. (0->0), fcn. (167->2), ass. (0->36)
t19 = cos(qJ(3));
t18 = sin(qJ(3));
t35 = qJD(2) * t18;
t30 = pkin(5) * t35;
t9 = t19 * qJD(1) - t30;
t41 = qJD(4) - t9;
t32 = qJD(2) * qJD(3);
t40 = -0.2e1 * t32;
t4 = -qJD(3) * pkin(3) + t41;
t10 = t19 * qJD(2) * pkin(5) + t18 * qJD(1);
t6 = qJD(3) * qJ(4) + t10;
t11 = -t19 * pkin(3) - t18 * qJ(4) - pkin(2);
t5 = qJD(2) * t11;
t39 = 0.2e1 * qJD(3) * t5;
t21 = qJD(2) ^ 2;
t38 = t19 * t21;
t20 = qJD(3) ^ 2;
t37 = t20 * t18;
t15 = t20 * t19;
t29 = t19 * t32;
t33 = qJD(1) * qJD(3);
t7 = pkin(5) * t29 + t18 * t33;
t16 = t18 ^ 2;
t36 = -t19 ^ 2 + t16;
t24 = pkin(3) * t18 - qJ(4) * t19;
t3 = qJD(3) * t24 - t18 * qJD(4);
t1 = qJD(2) * t3;
t31 = t18 * t38;
t27 = pkin(2) * t40;
t25 = t10 * qJD(3) - t7;
t23 = -pkin(5) * t20 - 0.2e1 * t1;
t14 = t19 * t33;
t2 = t14 + (qJD(4) - t30) * qJD(3);
t22 = t7 * t18 + t2 * t19 + (-t18 * t6 + t19 * t4) * qJD(3);
t8 = t24 * qJD(2);
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t15, -t37, 0, t15, t2 * t18 - t7 * t19 + (t18 * t4 + t19 * t6) * qJD(3); 0, 0, 0, 0, 0.2e1 * t18 * t29, t36 * t40, t15, -t37, 0, -pkin(5) * t15 + t18 * t27, pkin(5) * t37 + t19 * t27, t18 * t39 + t19 * t23, t22, t18 * t23 - t19 * t39, pkin(5) * t22 + t1 * t11 + t5 * t3; 0, 0, 0, 0, -t31, t36 * t21, 0, 0, 0, t21 * pkin(2) * t18 + t25, pkin(2) * t38 - t14 + (t9 + t30) * qJD(3), (-t18 * t5 + t19 * t8) * qJD(2) + t25, 0, t14 + (0.2e1 * qJD(4) - t9) * qJD(3) + (t19 * t5 + (-pkin(5) * qJD(3) + t8) * t18) * qJD(2), -t7 * pkin(3) + t2 * qJ(4) - t4 * t10 + t41 * t6 - t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t16 * t21 - t20, -t6 * qJD(3) + t35 * t5 + t7;];
tauc_reg = t12;
