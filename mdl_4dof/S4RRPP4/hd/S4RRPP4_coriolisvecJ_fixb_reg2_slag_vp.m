% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.38s
% Computational Cost: add. (286->105), mult. (787->147), div. (0->0), fcn. (335->2), ass. (0->80)
t82 = pkin(2) + pkin(3);
t61 = t82 * qJD(2);
t43 = sin(qJ(2));
t71 = qJD(1) * t43;
t35 = pkin(5) * t71;
t65 = qJ(4) * qJD(1);
t17 = t43 * t65 - t35;
t66 = qJD(3) - t17;
t7 = -t61 + t66;
t44 = cos(qJ(2));
t58 = t43 * qJ(3) + pkin(1);
t15 = t82 * t44 + t58;
t72 = qJD(1) * t15;
t4 = qJD(4) + t72;
t73 = qJD(4) + t4;
t83 = t73 * t43;
t76 = pkin(5) * qJD(1);
t36 = t44 * t76;
t19 = -t44 * t65 + t36;
t40 = qJD(2) * qJ(3);
t12 = t19 + t40;
t63 = qJD(1) * qJD(2);
t41 = t43 ^ 2;
t42 = t44 ^ 2;
t77 = -t41 + t42;
t13 = 0.2e1 * t77 * t63;
t47 = qJD(1) ^ 2;
t81 = t44 * t47;
t46 = qJD(2) ^ 2;
t80 = t46 * t43;
t37 = t46 * t44;
t79 = pkin(5) - qJ(4);
t59 = t44 * t63;
t69 = t43 * qJD(3);
t78 = qJ(3) * t59 + qJD(1) * t69;
t75 = qJ(3) * t44;
t74 = qJD(2) * pkin(2);
t22 = -t44 * pkin(2) - t58;
t14 = qJD(1) * t22;
t70 = qJD(2) * t43;
t68 = t43 * qJD(4);
t67 = t44 * qJD(4);
t64 = qJ(4) * qJD(2);
t62 = t44 * t64;
t26 = t79 * t44;
t60 = t43 * t63;
t1 = -t82 * t60 + t78;
t49 = -t82 * t43 + t75;
t3 = t49 * qJD(2) + t69;
t57 = qJD(1) * t3 + t1;
t56 = t4 + t72;
t55 = 0.2e1 * t14;
t54 = qJD(3) - t74;
t53 = -0.2e1 * t60;
t52 = t43 * t59;
t51 = pkin(2) * t43 - t75;
t11 = t51 * qJD(2) - t69;
t6 = pkin(2) * t60 - t78;
t50 = -pkin(5) * t46 - qJD(1) * t11 - t6;
t39 = qJD(2) * qJD(3);
t20 = -pkin(5) * t60 + t39;
t21 = t35 + t54;
t24 = t36 + t40;
t48 = t20 * t44 + (t21 * t44 + (-t24 + t36) * t43) * qJD(2);
t38 = 0.2e1 * t39;
t33 = pkin(5) * t59;
t32 = t43 * t81;
t30 = qJ(4) * t60;
t29 = -t41 * t47 - t46;
t28 = -0.2e1 * t52;
t27 = 0.2e1 * t52;
t25 = t79 * t43;
t23 = t77 * t47;
t18 = t51 * qJD(1);
t10 = qJD(2) * t26 - t68;
t9 = -t79 * t70 - t67;
t8 = t49 * qJD(1);
t5 = t33 + (-t62 - t68) * qJD(1);
t2 = t30 + t39 + (-pkin(5) * t70 - t67) * qJD(1);
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t13, t37, t28, -t80, 0, pkin(1) * t53 - pkin(5) * t37, -0.2e1 * pkin(1) * t59 + pkin(5) * t80, 0, 0, t27, t37, -t13, 0, t80, t28, t50 * t44 + t55 * t70, t48, -t55 * t44 * qJD(2) + t50 * t43, t48 * pkin(5) + t14 * t11 + t6 * t22, t27, -t13, -t37, t28, -t80, 0, t57 * t44 + (-t56 * t43 - t10) * qJD(2), t57 * t43 + (t56 * t44 + t9) * qJD(2), -t2 * t44 - t5 * t43 + (t12 * t43 - t44 * t7) * qJD(2) + (-t10 * t43 - t44 * t9 + (-t25 * t44 + t26 * t43) * qJD(2)) * qJD(1), t1 * t15 + t7 * t10 + t12 * t9 + t2 * t26 + t5 * t25 + t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t23, 0, t32, 0, 0, t47 * pkin(1) * t43, pkin(1) * t81, 0, 0, -t32, 0, t23, 0, 0, t32, (-t14 * t43 + t18 * t44) * qJD(1), ((t24 - t40) * t43 + (-t21 + t54) * t44) * qJD(1), t38 + (t14 * t44 + t18 * t43) * qJD(1), t20 * qJ(3) + t24 * qJD(3) - t14 * t18 + (t24 * t43 + (-t21 - t74) * t44) * t76, -t32, t23, 0, t32, 0, 0, t19 * qJD(2) - t33 + ((-t8 + t64) * t44 + t83) * qJD(1), -t17 * qJD(2) + t30 + t38 + (-t73 * t44 + (-pkin(5) * qJD(2) - t8) * t43) * qJD(1), 0, t2 * qJ(3) + t66 * t12 - t7 * t19 - t4 * t8 - t5 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, t29, -t24 * qJD(2) + t14 * t71 + t33, 0, 0, 0, 0, 0, 0, -t32, t29, 0, -t12 * qJD(2) + t33 + (-t62 - t83) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0.2e1 * t59, (-t41 - t42) * t47, (t12 * t44 + (t7 - t61) * t43) * qJD(1) + t78;];
tauc_reg = t16;
