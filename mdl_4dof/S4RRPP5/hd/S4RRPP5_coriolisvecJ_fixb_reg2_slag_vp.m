% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:38
% DurationCPUTime: 0.37s
% Computational Cost: add. (286->103), mult. (779->139), div. (0->0), fcn. (326->2), ass. (0->82)
t51 = cos(qJ(2));
t73 = qJD(1) * t51;
t38 = pkin(5) * t73;
t20 = pkin(3) * t73 + t38;
t84 = -qJD(4) - t20;
t70 = qJD(2) * qJ(3);
t12 = t70 - t84;
t49 = -pkin(2) - qJ(4);
t83 = qJD(2) * t49;
t82 = pkin(3) + pkin(5);
t50 = sin(qJ(2));
t74 = qJD(1) * t50;
t37 = pkin(5) * t74;
t76 = qJD(2) * pkin(2);
t23 = qJD(3) + t37 - t76;
t81 = t23 * t51;
t26 = -t38 - t70;
t80 = t26 * t50;
t48 = t51 ^ 2;
t53 = qJD(1) ^ 2;
t79 = t48 * t53;
t78 = t51 * t53;
t52 = qJD(2) ^ 2;
t43 = t52 * t50;
t44 = t52 * t51;
t69 = qJD(1) * qJD(2);
t66 = t51 * t69;
t33 = pkin(5) * t66;
t77 = pkin(3) * t66 + t33;
t75 = t50 * qJ(3);
t15 = t49 * t51 - pkin(1) - t75;
t4 = qJD(1) * t15;
t60 = -t51 * pkin(2) - t75;
t24 = -pkin(1) + t60;
t14 = qJD(1) * t24;
t72 = t50 * qJD(3);
t18 = -pkin(3) * t74 - t37;
t71 = qJD(3) - t18;
t30 = t82 * t51;
t68 = t82 * qJD(2);
t67 = t50 * t69;
t34 = pkin(2) * t67;
t59 = -qJ(3) * t51 + qJ(4) * t50;
t55 = t59 * qJD(2) - t51 * qJD(4) - t72;
t1 = t55 * qJD(1) + t34;
t40 = t50 * t76;
t2 = t40 + t55;
t65 = -qJD(1) * t2 - t1;
t64 = 0.2e1 * t4;
t63 = -0.2e1 * pkin(1) * t69;
t19 = t50 * t68;
t61 = t50 * t66;
t58 = -0.2e1 * qJD(2) * t14;
t56 = -t51 * t70 - t72;
t11 = t40 + t56;
t5 = t56 * qJD(1) + t34;
t57 = pkin(5) * t52 + qJD(1) * t11 + t5;
t46 = qJD(2) * qJD(3);
t22 = pkin(5) * t67 - t46;
t54 = -t22 * t51 + (t81 + (t26 + t38) * t50) * qJD(2);
t47 = t50 ^ 2;
t45 = 0.2e1 * t46;
t42 = t47 * t53;
t41 = pkin(2) * t74;
t36 = qJD(3) * t73;
t32 = t50 * t78;
t31 = -t42 - t52;
t29 = t82 * t50;
t28 = -0.2e1 * t61;
t27 = 0.2e1 * t61;
t25 = -t42 + t79;
t21 = qJD(2) * t30;
t17 = -qJ(3) * t73 + t41;
t16 = (-t47 + t48) * t69;
t13 = 0.2e1 * t16;
t10 = t59 * qJD(1) + t41;
t9 = -qJD(2) * qJD(4) + t77;
t8 = -qJD(1) * t19 + t46;
t7 = t71 + t83;
t6 = t14 * t74;
t3 = t4 * t73;
t35 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t13, t44, t28, -t43, 0, -pkin(5) * t44 + t50 * t63, pkin(5) * t43 + t51 * t63, 0, 0, 0, -t44, t43, t27, t13, t28, t54, t50 * t58 + t57 * t51, -t57 * t50 + t51 * t58, t54 * pkin(5) + t14 * t11 + t5 * t24, 0, t43, t44, t28, -0.2e1 * t16, t27, t9 * t50 + t8 * t51 + (-t12 * t50 + t51 * t7) * qJD(2) + (-t19 * t51 + t21 * t50 + (t29 * t51 - t30 * t50) * qJD(2)) * qJD(1), t65 * t50 + (-t64 * t51 - t19) * qJD(2), t65 * t51 + (t64 * t50 - t21) * qJD(2), t1 * t15 - t12 * t19 + t4 * t2 + t7 * t21 + t9 * t29 + t8 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t25, 0, t32, 0, 0, t53 * pkin(1) * t50, pkin(1) * t78, 0, 0, 0, 0, 0, -t32, -t25, t32, t36 + (t60 * qJD(2) - t80 - t81) * qJD(1), -t17 * t73 + t6, t45 + (t14 * t51 + t17 * t50) * qJD(1), -t22 * qJ(3) - t26 * qJD(3) - t14 * t17 + (-t80 + (-t23 - t76) * t51) * qJD(1) * pkin(5), 0, 0, 0, t32, t25, -t32, t36 + (-t18 - t7 + t83) * t73, -t18 * qJD(2) + t3 + t45 + (t10 - t68) * t74, (0.2e1 * qJD(4) + t20) * qJD(2) + (t10 * t51 - t4 * t50) * qJD(1) - t77, t8 * qJ(3) - t4 * t10 + t71 * t12 + t9 * t49 + t84 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t31, t26 * qJD(2) + t33 + t6, 0, 0, 0, 0, 0, 0, 0, t31, -t32, t4 * t74 + (-qJD(4) - t12) * qJD(2) + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t52 - t79, t3 + t46 + (-t82 * t74 + t7) * qJD(2);];
tauc_reg = t35;
