% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:14
% DurationCPUTime: 0.53s
% Computational Cost: add. (414->128), mult. (1066->168), div. (0->0), fcn. (489->2), ass. (0->87)
t91 = pkin(3) + pkin(4);
t69 = t91 * qJD(3);
t47 = sin(qJ(3));
t77 = qJD(2) * t47;
t39 = pkin(6) * t77;
t48 = cos(qJ(3));
t25 = t48 * qJD(1) - t39;
t11 = qJ(5) * t77 + t25;
t72 = qJD(4) - t11;
t4 = -t69 + t72;
t93 = qJD(4) - t25;
t66 = t47 * qJ(4) + pkin(2);
t19 = t91 * t48 + t66;
t78 = qJD(2) * t19;
t6 = qJD(5) + t78;
t79 = qJD(5) + t6;
t92 = t79 * t47;
t75 = t47 * qJD(1);
t26 = t48 * qJD(2) * pkin(6) + t75;
t44 = qJD(3) * qJ(4);
t20 = t26 + t44;
t85 = pkin(6) - qJ(5);
t30 = t85 * t48;
t13 = qJD(2) * t30 + t75;
t8 = t44 + t13;
t17 = -qJD(3) * pkin(3) + t93;
t70 = qJD(2) * qJD(3);
t45 = t47 ^ 2;
t46 = t48 ^ 2;
t82 = -t45 + t46;
t16 = 0.2e1 * t82 * t70;
t90 = t48 * t8;
t22 = t26 * qJD(3);
t89 = t22 * t47;
t88 = t22 * t48;
t51 = qJD(2) ^ 2;
t87 = t48 * t51;
t50 = qJD(3) ^ 2;
t86 = t50 * t47;
t41 = t50 * t48;
t67 = t48 * t70;
t74 = t47 * qJD(4);
t84 = qJ(4) * t67 + qJD(2) * t74;
t71 = qJD(1) * qJD(3);
t38 = t48 * t71;
t43 = qJD(3) * qJD(4);
t83 = t38 + 0.2e1 * t43;
t81 = pkin(6) * qJD(3);
t80 = qJ(4) * t48;
t27 = -t48 * pkin(3) - t66;
t18 = qJD(2) * t27;
t76 = qJD(3) * t47;
t73 = t48 * qJD(5);
t68 = t47 * t70;
t2 = -t91 * t68 + t84;
t54 = -t91 * t47 + t80;
t5 = t54 * qJD(3) + t74;
t65 = qJD(2) * t5 + t2;
t56 = pkin(3) * t47 - t80;
t24 = t56 * qJD(2);
t64 = t24 - t81;
t63 = t6 + t78;
t62 = t85 * qJD(3);
t61 = 0.2e1 * t18;
t59 = -0.2e1 * t68;
t58 = t47 * t67;
t57 = qJD(3) * t30;
t21 = -pkin(6) * t68 + t38;
t15 = t56 * qJD(3) - t74;
t7 = pkin(3) * t68 - t84;
t55 = -pkin(6) * t50 - qJD(2) * t15 - t7;
t14 = -t47 * qJD(5) + t57;
t10 = t21 + t43;
t53 = t10 * t48 + t89 + (t17 * t48 - t20 * t47) * qJD(3);
t52 = t21 * t48 + t89 + (-t25 * t48 - t26 * t47) * qJD(3);
t36 = t47 * t87;
t34 = qJ(5) * t68;
t33 = -t45 * t51 - t50;
t32 = -0.2e1 * t58;
t31 = 0.2e1 * t58;
t29 = t85 * t47;
t28 = t82 * t51;
t12 = -t47 * t62 - t73;
t9 = t54 * qJD(2);
t3 = qJD(2) * t14 + t47 * t71;
t1 = t34 + t38 + t43 + (-pkin(6) * t76 - t73) * qJD(2);
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t41, 0, t21 * t47 - t88 + (-t25 * t47 + t26 * t48) * qJD(3), 0, 0, 0, 0, 0, 0, -t86, 0, t41, t10 * t47 - t88 + (t17 * t47 + t20 * t48) * qJD(3), 0, 0, 0, 0, 0, 0, -t86, t41, 0, t1 * t47 - t3 * t48 + (t4 * t47 + t90) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t16, t41, t32, -t86, 0, pkin(2) * t59 - pkin(6) * t41, -0.2e1 * pkin(2) * t67 + pkin(6) * t86, t52, t52 * pkin(6), t31, t41, -t16, 0, t86, t32, t55 * t48 + t61 * t76, t53, -t61 * t48 * qJD(3) + t55 * t47, t53 * pkin(6) + t18 * t15 + t7 * t27, t31, -t16, -t41, t32, -t86, 0, t65 * t48 + (-t63 * t47 - t14) * qJD(3), t65 * t47 + (t63 * t48 + t12) * qJD(3), -t1 * t48 - t3 * t47 + (-t4 * t48 + t47 * t8) * qJD(3) + (-t12 * t48 - t14 * t47 + (-t29 * t48 + t30 * t47) * qJD(3)) * qJD(2), t1 * t30 + t8 * t12 + t4 * t14 + t2 * t19 + t3 * t29 + t6 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t28, 0, t36, 0, 0, t51 * pkin(2) * t47, pkin(2) * t87 - t38 + (t25 + t39) * qJD(3), 0, 0, -t36, 0, t28, 0, 0, t36, (t26 - t75) * qJD(3) + (-t18 * t47 + t64 * t48) * qJD(2), 0, -t25 * qJD(3) + (t18 * t48 + t64 * t47) * qJD(2) + t83, -t22 * pkin(3) + t10 * qJ(4) - t17 * t26 - t18 * t24 + t93 * t20, -t36, t28, 0, t36, 0, 0, (t13 - t75) * qJD(3) + (t92 + (-t9 - t62) * t48) * qJD(2), -t11 * qJD(3) + t34 + (-t79 * t48 + (-t9 - t81) * t47) * qJD(2) + t83, 0, t1 * qJ(4) - t4 * t13 - t3 * t91 - t6 * t9 + t72 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t33, t18 * t77 + (-t20 + t26) * qJD(3), 0, 0, 0, 0, 0, 0, -t36, t33, 0, (-t8 + t75) * qJD(3) + (t57 - t92) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0.2e1 * t67, (-t45 - t46) * t51, (t90 + (t4 - t69) * t47) * qJD(2) + t84;];
tauc_reg = t23;
