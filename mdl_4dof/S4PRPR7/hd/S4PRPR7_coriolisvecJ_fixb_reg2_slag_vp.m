% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (88->41), mult. (253->62), div. (0->0), fcn. (120->4), ass. (0->38)
t7 = sin(qJ(4));
t5 = t7 ^ 2;
t9 = cos(qJ(4));
t6 = t9 ^ 2;
t33 = t5 + t6;
t8 = sin(qJ(2));
t38 = t33 * t8;
t11 = -pkin(2) - pkin(5);
t37 = qJD(2) * t11;
t27 = t8 * qJD(1);
t25 = qJD(2) * qJ(3);
t4 = t25 + t27;
t20 = -t4 + t27;
t36 = qJD(4) * (-t20 + t25);
t16 = t20 * qJD(2);
t10 = cos(qJ(2));
t26 = t10 * qJD(1);
t2 = (qJD(3) + t26) * qJD(2);
t35 = t2 * t8;
t34 = t5 - t6;
t32 = t10 * t4;
t13 = qJD(2) ^ 2;
t31 = t13 * t8;
t30 = t13 * t10;
t12 = qJD(4) ^ 2;
t29 = t12 + t13;
t28 = qJD(2) * pkin(2);
t24 = qJD(2) * qJD(4);
t23 = t9 * t13 * t7;
t22 = t9 * t24;
t21 = t10 * t29;
t19 = t7 * t22;
t18 = qJD(3) - t26;
t17 = t2 * qJ(3) + t4 * qJD(3);
t14 = t18 * qJD(2) - t11 * t12 + t2;
t3 = t18 - t28;
t1 = t18 + t37;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t30, t35 + (-t20 * t10 + t3 * t8) * qJD(2), 0, 0, 0, 0, 0, 0, t7 * t21 + 0.2e1 * t8 * t22, -0.2e1 * t7 * t8 * t24 + t9 * t21, -t33 * t31, t35 + (t32 + (t1 - t26) * t38) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), (-t32 + (-t3 - t28) * t8) * qJD(1) + t17, -0.2e1 * t19, 0.2e1 * t34 * t24, -t12 * t7, 0.2e1 * t19, -t12 * t9, 0, t14 * t7 + t9 * t36, t14 * t9 - t7 * t36, 0, (-t32 + (-t1 + t37) * t38) * qJD(1) + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t16, 0, 0, 0, 0, 0, 0, -t29 * t7, -t29 * t9, 0, (t33 * t27 - t4) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t34 * t13, 0, -t23, 0, 0, t9 * t16, -t7 * t16, 0, 0;];
tauc_reg = t15;
