% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:59
% DurationCPUTime: 0.21s
% Computational Cost: add. (175->58), mult. (511->94), div. (0->0), fcn. (235->2), ass. (0->45)
t22 = sin(qJ(3));
t41 = qJD(2) * t22;
t36 = pkin(5) * t41;
t23 = cos(qJ(3));
t39 = t23 * qJD(1);
t11 = -t36 + t39;
t48 = qJD(4) - t11;
t4 = -qJD(3) * pkin(3) + t48;
t40 = t22 * qJD(1);
t12 = t23 * qJD(2) * pkin(5) + t40;
t6 = qJD(3) * qJ(4) + t12;
t13 = -t23 * pkin(3) - t22 * qJ(4) - pkin(2);
t5 = qJD(2) * t13;
t47 = 0.2e1 * qJD(3) * t5;
t8 = t12 * qJD(3);
t46 = t8 * t22;
t45 = t8 * t23;
t25 = qJD(2) ^ 2;
t44 = t23 * t25;
t24 = qJD(3) ^ 2;
t43 = t24 * t22;
t19 = t24 * t23;
t20 = t22 ^ 2;
t42 = t23 ^ 2 - t20;
t29 = pkin(3) * t22 - qJ(4) * t23;
t3 = t29 * qJD(3) - t22 * qJD(4);
t1 = qJD(2) * t3;
t37 = qJD(2) * qJD(3);
t35 = t22 * t37;
t10 = t29 * qJD(2);
t34 = -pkin(5) * qJD(3) + t10;
t32 = -0.2e1 * pkin(2) * t37;
t30 = t23 * t35;
t28 = -pkin(5) * t24 - 0.2e1 * t1;
t18 = qJD(3) * t39;
t2 = t18 + (qJD(4) - t36) * qJD(3);
t27 = t2 * t23 + t46 + (-t22 * t6 + t23 * t4) * qJD(3);
t7 = -pkin(5) * t35 + t18;
t26 = t46 + t7 * t23 + (-t11 * t23 - t12 * t22) * qJD(3);
t17 = t22 * t44;
t16 = -0.2e1 * t30;
t15 = 0.2e1 * t30;
t14 = t42 * t25;
t9 = t42 * t37;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t19, 0, t7 * t22 - t45 + (-t11 * t22 + t12 * t23) * qJD(3), 0, 0, 0, 0, 0, 0, -t43, 0, t19, t2 * t22 - t45 + (t22 * t4 + t23 * t6) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0.2e1 * t9, t19, t16, -t43, 0, -pkin(5) * t19 + t22 * t32, pkin(5) * t43 + t23 * t32, t26, t26 * pkin(5), t15, t19, -0.2e1 * t9, 0, t43, t16, t22 * t47 + t28 * t23, t27, t28 * t22 - t23 * t47, t27 * pkin(5) + t1 * t13 + t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t14, 0, t17, 0, 0, t25 * pkin(2) * t22, pkin(2) * t44 - t18 + (t11 + t36) * qJD(3), 0, 0, -t17, 0, t14, 0, 0, t17, (t12 - t40) * qJD(3) + (-t22 * t5 + t34 * t23) * qJD(2), 0, t18 + (0.2e1 * qJD(4) - t11) * qJD(3) + (t34 * t22 + t23 * t5) * qJD(2), -t8 * pkin(3) + t2 * qJ(4) - t5 * t10 - t4 * t12 + t48 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t20 * t25 - t24, t5 * t41 + (t12 - t6) * qJD(3);];
tauc_reg = t21;
