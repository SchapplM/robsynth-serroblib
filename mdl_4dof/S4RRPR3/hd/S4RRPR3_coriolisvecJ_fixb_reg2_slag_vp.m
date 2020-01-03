% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:35
% DurationCPUTime: 0.32s
% Computational Cost: add. (471->69), mult. (1105->113), div. (0->0), fcn. (615->6), ass. (0->71)
t37 = sin(pkin(7));
t40 = sin(qJ(2));
t67 = pkin(1) * qJD(1);
t63 = t40 * t67;
t60 = t37 * t63;
t66 = pkin(1) * qJD(2);
t62 = qJD(1) * t66;
t38 = cos(pkin(7));
t42 = cos(qJ(2));
t72 = t38 * t42;
t14 = -qJD(2) * t60 + t62 * t72;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t34 = qJD(1) + qJD(2);
t25 = t34 * pkin(2) + t42 * t67;
t11 = t37 * t25 + t38 * t63;
t9 = t34 * pkin(6) + t11;
t4 = t41 * qJD(3) - t39 * t9;
t2 = t4 * qJD(4) + t41 * t14;
t5 = t39 * qJD(3) + t41 * t9;
t3 = -t5 * qJD(4) - t39 * t14;
t80 = t2 * t41 - t3 * t39 + (-t39 * t5 - t4 * t41) * qJD(4);
t73 = t38 * t40;
t53 = t37 * t42 + t73;
t49 = pkin(1) * t53;
t19 = qJD(2) * t49;
t13 = qJD(1) * t19;
t65 = qJD(4) * t41;
t10 = t38 * t25 - t60;
t8 = -t34 * pkin(3) - t10;
t79 = t13 * t39 + t8 * t65;
t18 = qJD(1) * t49;
t78 = t18 * t34;
t77 = t19 * t34;
t74 = t37 * t40;
t48 = pkin(1) * (t72 - t74);
t20 = qJD(1) * t48;
t76 = t20 * t34;
t21 = qJD(2) * t48;
t75 = t21 * t34;
t43 = qJD(4) ^ 2;
t71 = t43 * t39;
t31 = t42 * pkin(1) + pkin(2);
t70 = pkin(1) * t73 + t37 * t31;
t35 = t39 ^ 2;
t36 = t41 ^ 2;
t69 = t35 - t36;
t68 = t35 + t36;
t33 = t34 ^ 2;
t64 = t39 * t33 * t41;
t61 = -t34 * t8 - t14;
t59 = t39 * t34 * t65;
t58 = t39 * t4 - t41 * t5;
t57 = (-qJD(2) + t34) * t67;
t56 = (-qJD(1) - t34) * t66;
t17 = pkin(6) + t70;
t55 = t17 * t43 + t77;
t29 = t37 * pkin(2) + pkin(6);
t54 = t29 * t43 - t78;
t50 = -pkin(1) * t74 + t38 * t31;
t16 = -pkin(3) - t50;
t52 = qJD(4) * (t16 * t34 - t21);
t30 = -t38 * pkin(2) - pkin(3);
t51 = qJD(4) * (t30 * t34 + t20);
t47 = t53 * t62;
t32 = t43 * t41;
t24 = -0.2e1 * t59;
t23 = 0.2e1 * t59;
t15 = -0.2e1 * t69 * t34 * qJD(4);
t6 = t8 * qJD(4) * t39;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t56, t42 * t56, 0, 0, 0, 0, 0, 0, 0, 0, -t47 - t77, -t14 - t75, 0, -t10 * t19 + t11 * t21 - t13 * t50 + t14 * t70, t23, t15, t32, t24, -t71, 0, t6 + t39 * t52 + (-t13 - t55) * t41, t55 * t39 + t41 * t52 + t79, t68 * t75 + t80, t13 * t16 + t17 * t80 + t8 * t19 - t58 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t57, t42 * t57, 0, 0, 0, 0, 0, 0, 0, 0, -t47 + t78, -t14 + t76, 0, t10 * t18 - t11 * t20 + (-t13 * t38 + t14 * t37) * pkin(2), t23, t15, t32, t24, -t71, 0, t6 + t39 * t51 + (-t13 - t54) * t41, t54 * t39 + t41 * t51 + t79, -t68 * t76 + t80, t13 * t30 - t8 * t18 + t58 * t20 + t29 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t32, 0, -t58 * qJD(4) + t2 * t39 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t69 * t33, 0, t64, 0, 0, t61 * t39, t61 * t41, 0, 0;];
tauc_reg = t1;
