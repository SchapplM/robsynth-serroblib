% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP3
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
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:05
% EndTime: 2021-01-14 22:27:06
% DurationCPUTime: 0.17s
% Computational Cost: add. (155->62), mult. (421->97), div. (0->0), fcn. (192->2), ass. (0->45)
t35 = qJD(3) * pkin(3);
t17 = cos(qJ(3));
t13 = t17 * qJD(1);
t16 = sin(qJ(3));
t37 = -qJ(4) - pkin(5);
t9 = t37 * t16;
t4 = qJD(2) * t9 + t13;
t3 = t4 + t35;
t44 = t3 - t4;
t14 = t16 ^ 2;
t43 = pkin(3) * t14;
t42 = t17 * pkin(3);
t10 = t37 * t17;
t32 = t16 * qJD(1);
t6 = -qJD(2) * t10 + t32;
t41 = t17 * t6;
t19 = qJD(2) ^ 2;
t40 = t17 * t19;
t18 = qJD(3) ^ 2;
t39 = t18 * t16;
t38 = t18 * t17;
t15 = t17 ^ 2;
t36 = t14 - t15;
t12 = -pkin(2) - t42;
t33 = qJD(2) * t12;
t8 = qJD(4) + t33;
t34 = -qJD(4) - t8;
t31 = t16 * qJD(3);
t30 = t17 * qJD(4);
t29 = qJD(2) * qJD(3);
t28 = qJD(3) * qJD(1);
t27 = t16 * t40;
t26 = 0.2e1 * t29;
t25 = qJ(4) * t31;
t24 = t16 * t29;
t23 = qJD(3) * t37;
t22 = -0.2e1 * pkin(2) * t29;
t21 = t17 * t26;
t20 = t17 * t23;
t7 = -t16 * qJD(4) + t20;
t11 = pkin(5) * t24;
t5 = t16 * t23 + t30;
t2 = qJD(2) * t7 - t16 * t28;
t1 = t17 * t28 - t11 + (-t25 + t30) * qJD(2);
t45 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, -t39, -t38, 0, t1 * t16 + t2 * t17 + (-t16 * t3 + t41) * qJD(3); 0, 0, 0, 0, t16 * t21, -t36 * t26, t38, -t39, 0, -pkin(5) * t38 + t16 * t22, pkin(5) * t39 + t17 * t22, (t7 + (t8 + (t12 - 0.2e1 * t42) * qJD(2)) * t16) * qJD(3), (t17 * t8 - t5 + (t12 * t17 + 0.2e1 * t43) * qJD(2)) * qJD(3), t1 * t17 - t2 * t16 + (-t16 * t6 - t17 * t3) * qJD(3) + (-t16 * t7 + t17 * t5 + (t10 * t16 - t17 * t9) * qJD(3)) * qJD(2), -t1 * t10 + t2 * t9 + t3 * t7 + t6 * t5 + (t8 + t33) * pkin(3) * t31; 0, 0, 0, 0, -t27, t36 * t19, 0, 0, 0, t19 * pkin(2) * t16, -qJD(2) * pkin(5) * t31 + pkin(2) * t40 + t11, pkin(3) * t27 + (t6 - t32) * qJD(3) + (t34 * t16 + t20) * qJD(2), -t19 * t43 + t11 + (t4 - t13) * qJD(3) + (t34 * t17 + t25) * qJD(2), (-t35 + t44) * t17 * qJD(2), t44 * t6 + (-t8 * t16 * qJD(2) + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t24, t21, (-t14 - t15) * t19, (-t41 + (t3 + t35) * t16) * qJD(2);];
tauc_reg = t45;
