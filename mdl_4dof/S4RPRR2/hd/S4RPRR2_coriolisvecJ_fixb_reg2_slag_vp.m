% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:14
% DurationCPUTime: 0.25s
% Computational Cost: add. (454->58), mult. (1089->86), div. (0->0), fcn. (588->6), ass. (0->54)
t26 = cos(pkin(7)) * pkin(1) + pkin(2);
t23 = t26 * qJD(1);
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t65 = sin(pkin(7)) * pkin(1);
t52 = qJD(1) * t65;
t13 = t37 * t23 - t35 * t52;
t11 = t13 * qJD(3);
t36 = cos(qJ(4));
t14 = t35 * t23 + t37 * t52;
t29 = qJD(1) + qJD(3);
t10 = t29 * pkin(6) + t14;
t34 = sin(qJ(4));
t4 = t36 * qJD(2) - t34 * t10;
t2 = t4 * qJD(4) + t36 * t11;
t5 = t34 * qJD(2) + t36 * t10;
t3 = -t5 * qJD(4) - t34 * t11;
t68 = t2 * t36 - t3 * t34 + (-t34 * t5 - t36 * t4) * qJD(4);
t57 = t35 * t26 + t37 * t65;
t67 = t14 * qJD(3);
t54 = qJD(4) * t36;
t64 = t29 * pkin(3);
t9 = -t13 - t64;
t66 = t34 * t67 + t9 * t54;
t63 = t13 * t29;
t62 = t14 * t29;
t43 = t37 * t26 - t35 * t65;
t15 = t43 * qJD(3);
t61 = t15 * t29;
t16 = t57 * qJD(3);
t60 = t16 * t29;
t38 = qJD(4) ^ 2;
t58 = t38 * t34;
t30 = t34 ^ 2;
t31 = t36 ^ 2;
t56 = t30 - t31;
t55 = t30 + t31;
t28 = t29 ^ 2;
t53 = t34 * t28 * t36;
t51 = -t29 * t9 - t11;
t49 = t34 * t29 * t54;
t48 = t34 * t4 - t36 * t5;
t47 = pkin(6) * t38 - t62;
t19 = pkin(6) + t57;
t46 = t19 * t38 + t60;
t45 = qJD(4) * (t13 - t64);
t18 = -pkin(3) - t43;
t44 = qJD(4) * (t18 * t29 - t15);
t27 = t38 * t36;
t21 = -0.2e1 * t49;
t20 = 0.2e1 * t49;
t17 = -0.2e1 * t56 * t29 * qJD(4);
t6 = t9 * qJD(4) * t34;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 - t60, -t11 - t61, 0, t11 * t57 - t13 * t16 + t14 * t15 - t43 * t67, t20, t17, t27, t21, -t58, 0, t6 + t34 * t44 + (-t67 - t46) * t36, t46 * t34 + t36 * t44 + t66, t55 * t61 + t68, -t48 * t15 + t9 * t16 + t18 * t67 + t19 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t27, 0, -t48 * qJD(4) + t2 * t34 + t3 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 + t62, -t11 + t63, 0, 0, t20, t17, t27, t21, -t58, 0, t6 + t34 * t45 + (-t67 - t47) * t36, t47 * t34 + t36 * t45 + t66, -t55 * t63 + t68, -pkin(3) * t67 + pkin(6) * t68 + t48 * t13 - t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t56 * t28, 0, t53, 0, 0, t51 * t34, t51 * t36, 0, 0;];
tauc_reg = t1;
