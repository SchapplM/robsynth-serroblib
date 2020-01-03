% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:43
% DurationCPUTime: 0.34s
% Computational Cost: add. (435->86), mult. (1197->132), div. (0->0), fcn. (853->6), ass. (0->64)
t47 = cos(qJ(2));
t58 = t47 * qJD(1);
t51 = qJD(3) - t58;
t43 = cos(pkin(7));
t46 = cos(qJ(4));
t68 = t46 * t43;
t42 = sin(pkin(7));
t44 = sin(qJ(4));
t69 = t44 * t42;
t27 = -t68 + t69;
t29 = (qJD(3) + t58) * qJD(2);
t75 = t27 * t29;
t45 = sin(qJ(2));
t18 = t27 * t45;
t59 = t45 * qJD(1);
t34 = qJD(2) * qJ(3) + t59;
t64 = t42 ^ 2 + t43 ^ 2;
t74 = t64 * t34;
t28 = t46 * t42 + t44 * t43;
t62 = qJD(2) * t28;
t73 = t62 ^ 2;
t65 = pkin(5) + qJ(3);
t30 = t65 * t42;
t31 = t65 * t43;
t10 = -t44 * t30 + t46 * t31;
t72 = t10 * qJD(4) + t51 * t28;
t49 = t27 * t47;
t9 = -t46 * t30 - t44 * t31;
t71 = -qJD(1) * t49 + t27 * qJD(3) - t9 * qJD(4);
t56 = qJD(2) * t68;
t57 = qJD(2) * t69;
t20 = -t56 + t57;
t70 = t62 * t20;
t48 = qJD(2) ^ 2;
t67 = t48 * t45;
t66 = t48 * t47;
t63 = qJD(2) * pkin(2);
t38 = -t43 * pkin(3) - pkin(2);
t61 = qJD(2) * t38;
t60 = qJD(2) * t45;
t55 = t64 * t47;
t54 = t64 * t29;
t53 = pkin(5) * qJD(2) + t34;
t52 = t64 * qJD(3);
t13 = t53 * t42;
t14 = t53 * t43;
t5 = -t46 * t13 - t44 * t14;
t6 = -t44 * t13 + t46 * t14;
t50 = t28 * t29;
t17 = t28 * t45;
t25 = t28 * qJD(4);
t39 = qJD(2) * t59;
t33 = qJD(4) * t56;
t32 = t51 - t63;
t26 = t51 + t61;
t24 = t27 * qJD(4);
t19 = t20 ^ 2;
t16 = qJD(2) * t25;
t15 = qJD(4) * t57 - t33;
t8 = qJD(4) * t18 - t47 * t62;
t7 = -qJD(2) * t49 - qJD(4) * t17;
t2 = -t6 * qJD(4) - t50;
t1 = t5 * qJD(4) - t75;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t67, t42 * t67, t64 * t66, t45 * t54 + (t32 * t45 + (-t59 + t74) * t47) * qJD(2), 0, 0, 0, 0, 0, 0, t8 * qJD(4) - t47 * t16 + t20 * t60, -t7 * qJD(4) + t47 * t15 + t60 * t62, -t17 * t15 + t18 * t16 - t7 * t20 - t62 * t8, -t1 * t18 - t2 * t17 + t5 * t8 + t6 * t7 + (t26 - t58) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 + (-qJD(1) * t55 + t52) * qJD(2), t34 * t52 + qJ(3) * t54 + ((-t32 - t63) * t45 - t34 * t55) * qJD(1), -t15 * t28 - t24 * t62, t15 * t27 - t28 * t16 + t24 * t20 - t25 * t62, -t24 * qJD(4), t16 * t27 + t20 * t25, -t25 * qJD(4), 0, t38 * t16 + t26 * t25 - t72 * qJD(4) + (qJD(2) * t27 - t20) * t59, t71 * qJD(4) - t38 * t15 - t26 * t24, -t1 * t27 - t10 * t16 + t9 * t15 - t2 * t28 + t71 * t20 + t5 * t24 - t6 * t25 + t62 * t72, t1 * t10 + t2 * t9 - t71 * t6 - t72 * t5 + (-t26 + t61) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t48, -qJD(2) * t74 + t39, 0, 0, 0, 0, 0, 0, 0.2e1 * t62 * qJD(4), t33 + (-t20 - t57) * qJD(4), -t19 - t73, t6 * t20 + t5 * t62 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t19 + t73, t33 + (t20 - t57) * qJD(4), -t70, 0, 0, -t26 * t62 - t50, t26 * t20 + t75, 0, 0;];
tauc_reg = t3;
