% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:57
% DurationCPUTime: 0.41s
% Computational Cost: add. (451->90), mult. (1139->135), div. (0->0), fcn. (788->4), ass. (0->64)
t45 = sin(pkin(8));
t43 = t45 ^ 2;
t46 = cos(pkin(8));
t72 = t46 ^ 2 + t43;
t82 = t72 * qJD(2) * qJD(3);
t48 = cos(qJ(5));
t69 = qJD(2) * t45;
t59 = t48 * t69;
t47 = sin(qJ(5));
t68 = qJD(2) * t46;
t60 = t47 * t68;
t24 = t59 - t60;
t28 = t45 * t48 - t46 * t47;
t57 = t45 * qJ(4) + pkin(2);
t81 = (t46 * pkin(3) + t57) * qJD(2);
t27 = t45 * t47 + t46 * t48;
t70 = qJD(2) * t27;
t80 = t70 ^ 2;
t79 = t24 ^ 2;
t35 = qJD(5) * t60;
t14 = qJD(5) * t59 - t35;
t20 = t27 * qJD(5);
t78 = -t28 * t14 + t20 * t70;
t77 = t24 * t70;
t64 = qJ(3) * qJD(2);
t30 = t45 * qJD(1) + t46 * t64;
t76 = t30 * t46;
t73 = -pkin(6) + qJ(3);
t25 = (pkin(3) + pkin(4)) * t46 + t57;
t71 = qJD(2) * t25;
t67 = qJD(3) * t45;
t66 = qJD(4) * t45;
t65 = t20 * qJD(5);
t62 = qJD(2) * qJD(4);
t61 = qJ(3) * t82 + qJD(3) * t76;
t58 = t45 * t62;
t29 = t46 * qJD(1) - t45 * t64;
t55 = 0.2e1 * t70;
t26 = qJD(4) - t29;
t12 = -pkin(6) * t69 + t26;
t16 = -pkin(6) * t68 + t30;
t5 = t48 * t12 - t47 * t16;
t6 = t47 * t12 + t48 * t16;
t13 = qJD(2) * t20;
t21 = t28 * qJD(5);
t54 = -t27 * t13 + t24 * t21;
t32 = t73 * t45;
t33 = t73 * t46;
t9 = t48 * t32 - t47 * t33;
t10 = t47 * t32 + t48 * t33;
t52 = t24 * qJD(3);
t51 = t27 * qJD(3);
t50 = qJD(2) ^ 2;
t49 = qJD(5) ^ 2;
t31 = t72 * t50;
t19 = 0.2e1 * t82;
t17 = qJD(3) - t81;
t15 = t21 * qJD(5);
t11 = -qJD(3) + t71;
t4 = t28 * qJD(3) - t10 * qJD(5);
t3 = t9 * qJD(5) + t51;
t2 = -t6 * qJD(5) + t52;
t1 = qJD(2) * t51 + t5 * qJD(5);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t65, t54 + t78, t1 * t28 - t2 * t27 - t6 * t20 - t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t29 * t67 + t61, 0, 0, 0, 0, 0, 0, 0.2e1 * t46 * t58, t19, 0.2e1 * t43 * t62, t26 * t67 + (-t17 + t81) * t66 + t61, -t13 * t28 - t24 * t20, -t54 + t78, -t65, t14 * t27 + t21 * t70, -t15, 0, t4 * qJD(5) + t11 * t21 + t25 * t14 + t55 * t66, -t3 * qJD(5) - t11 * t20 - t25 * t13 + (qJD(2) * t28 + t24) * t66, -t1 * t27 - t10 * t14 + t9 * t13 - t2 * t28 + t5 * t20 - t6 * t21 - t4 * t24 - t3 * t70, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + (t11 + t71) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, (t29 * t45 - t76) * qJD(2), 0, 0, 0, 0, 0, 0, 0, -t31, 0, (-t76 + (-qJD(4) - t26) * t45) * qJD(2), 0, 0, 0, 0, 0, 0, t35 + (-t24 - t59) * qJD(5), t55 * qJD(5), t79 + t80, -t5 * t24 - t6 * t70 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t50 * t46, 0, -t43 * t50, (qJD(3) + t17) * t69, 0, 0, 0, 0, 0, 0, -t49 * t47 - t69 * t70, -t24 * t69 - t49 * t48, t48 * t13 - t47 * t14 + (t24 * t47 - t48 * t70) * qJD(5), -t11 * t69 + t1 * t47 + t2 * t48 + (-t47 * t5 + t48 * t6) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t79 - t80, 0, -t77, t35 + (t24 - t59) * qJD(5), 0, -t11 * t24 + t52, -(qJD(3) - t11) * t70, 0, 0;];
tauc_reg = t7;
