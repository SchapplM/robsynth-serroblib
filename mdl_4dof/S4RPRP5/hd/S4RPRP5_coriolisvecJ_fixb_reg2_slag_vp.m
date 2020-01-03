% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:03
% DurationCPUTime: 0.38s
% Computational Cost: add. (615->111), mult. (1769->137), div. (0->0), fcn. (1234->4), ass. (0->75)
t52 = cos(pkin(6));
t87 = cos(qJ(3));
t73 = t87 * t52;
t68 = qJD(1) * t73;
t51 = sin(pkin(6));
t53 = sin(qJ(3));
t84 = t53 * t51;
t74 = qJD(1) * t84;
t29 = -t68 + t74;
t26 = t29 ^ 2;
t37 = t87 * t51 + t53 * t52;
t91 = t37 * qJD(1);
t90 = t91 ^ 2;
t93 = -t26 - t90;
t92 = -t90 + t26;
t83 = pkin(5) + qJ(2);
t41 = t83 * t51;
t38 = qJD(1) * t41;
t42 = t83 * t52;
t39 = qJD(1) * t42;
t19 = -t53 * t38 + t87 * t39;
t60 = t37 * qJD(2);
t58 = qJD(1) * t60;
t6 = t19 * qJD(3) + t58;
t64 = -t87 * t41 - t53 * t42;
t89 = t6 * t64;
t48 = -t52 * pkin(2) - pkin(1);
t40 = t48 * qJD(1) + qJD(2);
t9 = t29 * pkin(3) - qJ(4) * t91 + t40;
t88 = t9 * t91;
t86 = t91 * t29;
t85 = t53 * t39;
t15 = qJD(3) * qJ(4) + t19;
t82 = t15 - t19;
t69 = qJD(3) * t87;
t81 = qJD(2) * t68 - t38 * t69;
t80 = t51 ^ 2 + t52 ^ 2;
t34 = t37 * qJD(3);
t79 = qJD(3) * t34;
t63 = t73 - t84;
t10 = t63 * qJD(2) + t64 * qJD(3);
t78 = t10 * qJD(3);
t21 = -t53 * t41 + t87 * t42;
t11 = t21 * qJD(3) + t60;
t77 = t11 * qJD(3);
t18 = -t87 * t38 - t85;
t76 = qJD(4) - t18;
t75 = qJD(1) * qJD(2);
t72 = t80 * qJD(1) ^ 2;
t71 = t51 * t75;
t70 = t18 + t85;
t44 = qJD(3) * t68;
t22 = qJD(3) * t74 - t44;
t23 = qJD(1) * t34;
t67 = t23 * pkin(3) + t22 * qJ(4);
t66 = -t23 * t63 + t29 * t34;
t65 = 0.2e1 * t80 * t75;
t62 = t53 * t71 - t81;
t59 = -t10 * t29 + t11 * t91 - t21 * t23 + t22 * t64 + t6 * t37;
t33 = qJD(3) * t84 - t52 * t69;
t57 = -t22 * t63 - t37 * t23 + t33 * t29 - t34 * t91;
t56 = 0.2e1 * t91 * qJD(3);
t55 = qJD(2) * t91;
t24 = t33 * qJD(3);
t17 = -pkin(3) * t63 - t37 * qJ(4) + t48;
t16 = pkin(3) * t91 + t29 * qJ(4);
t14 = t44 + (t29 - t74) * qJD(3);
t13 = -t44 + (t29 + t74) * qJD(3);
t12 = -qJD(3) * pkin(3) + t76;
t7 = t34 * pkin(3) + t33 * qJ(4) - t37 * qJD(4);
t5 = (-qJD(3) * t39 - t71) * t53 + t81;
t4 = (qJD(4) - t85) * qJD(3) - t62;
t3 = -qJD(4) * t91 + t67;
t1 = -t22 * t37 - t33 * t91;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, qJ(2) * t65, t1, t57, -t24, t66, -t79, 0, t48 * t23 + t40 * t34 - t77, -t48 * t22 - t40 * t33 - t78, t18 * t33 - t19 * t34 + t5 * t63 + t59, t19 * t10 - t18 * t11 + t5 * t21 - t89, t1, -t24, -t57, 0, t79, t66, t17 * t23 + t7 * t29 - t3 * t63 + t9 * t34 - t77, -t12 * t33 - t15 * t34 + t4 * t63 + t59, t17 * t22 - t3 * t37 + t9 * t33 - t7 * t91 + t78, t15 * t10 + t12 * t11 + t3 * t17 + t4 * t21 + t9 * t7 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -qJ(2) * t72, 0, 0, 0, 0, 0, 0, t56, -t13, t93, t18 * t91 + t19 * t29, 0, 0, 0, 0, 0, 0, t56, t93, t13, t15 * t29 + (-qJD(4) - t12) * t91 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t92, t14, -t86, 0, 0, -t40 * t91 - t55, t70 * qJD(3) + t40 * t29 + t62, 0, 0, t86, t14, t92, 0, 0, -t86, -t16 * t29 - t55 - t88, pkin(3) * t22 - t23 * qJ(4) + t82 * t91 + (t12 - t76) * t29, t16 * t91 - t9 * t29 + (0.2e1 * qJD(4) - t70) * qJD(3) - t62, -t6 * pkin(3) + t4 * qJ(4) - t12 * t19 + t76 * t15 - t9 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t14, -qJD(3) ^ 2 - t90, -t82 * qJD(3) + t58 + t88;];
tauc_reg = t2;
