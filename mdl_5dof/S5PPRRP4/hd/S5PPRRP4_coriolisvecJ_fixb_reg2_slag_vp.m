% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (369->107), mult. (961->151), div. (0->0), fcn. (531->4), ass. (0->85)
t66 = (qJD(3) * qJD(4));
t90 = -2 * t66;
t41 = sin(qJ(3));
t44 = qJD(4) ^ 2;
t45 = qJD(3) ^ 2;
t89 = (t44 + t45) * t41;
t43 = cos(qJ(3));
t27 = qJD(3) * pkin(6) + t41 * qJD(2);
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t72 = t42 * qJD(1);
t10 = -t40 * t27 - t72;
t35 = t40 * qJD(1);
t83 = t42 * t27;
t11 = -t35 + t83;
t49 = t10 * t40 - t11 * t42;
t88 = t49 * t43;
t67 = qJD(2) * qJD(3);
t31 = t41 * t67;
t60 = t40 * t66;
t17 = pkin(4) * t60 + t31;
t76 = qJD(4) * pkin(4);
t57 = qJ(5) * qJD(3) + t27;
t8 = -t57 * t40 - t72;
t5 = t8 + t76;
t9 = t57 * t42 - t35;
t51 = t40 * t5 - t42 * t9;
t87 = t51 * qJD(3) + t17;
t86 = t5 - t8;
t38 = t40 ^ 2;
t85 = pkin(4) * t38;
t84 = t42 * pkin(4);
t36 = t44 * t40;
t37 = t44 * t42;
t82 = t45 * t43;
t81 = -qJ(5) - pkin(6);
t39 = t42 ^ 2;
t80 = t38 - t39;
t79 = t38 + t39;
t77 = qJD(3) * pkin(3);
t75 = qJD(3) * t40;
t74 = qJD(3) * t42;
t73 = qJD(4) * t40;
t71 = t42 * qJD(5);
t70 = t43 * qJD(2);
t69 = qJ(5) * qJD(4);
t68 = qJD(1) * qJD(4);
t29 = t40 * t45 * t42;
t62 = t43 * t67;
t65 = t27 * t73 + (-t62 + t68) * t42;
t34 = -pkin(3) - t84;
t64 = t40 * t69;
t63 = t42 * t69;
t61 = qJD(4) * t70;
t59 = t42 * t66;
t58 = qJD(4) * t81;
t56 = t43 * t90;
t55 = -qJD(5) - t70;
t54 = t40 * t59;
t53 = t79 * t70;
t32 = t40 * t68;
t52 = -qJD(4) * t83 + t32;
t50 = -t17 + t31;
t28 = -t70 - t77;
t48 = qJD(4) * (t28 - t77);
t1 = (-t64 + t71) * qJD(3) - t65;
t2 = (t55 * t40 - t63) * qJD(3) + t52;
t47 = t1 * t42 - t2 * t40 + (-t40 * t9 - t42 * t5) * qJD(4);
t4 = -t40 * t62 + t52;
t46 = -t65 * t42 - t4 * t40 + (-t10 * t42 - t11 * t40) * qJD(4);
t25 = t42 * t61;
t24 = t40 * t61;
t23 = -0.2e1 * t54;
t22 = 0.2e1 * t54;
t21 = t81 * t42;
t20 = t81 * t40;
t19 = t80 * t45;
t16 = t79 * t82;
t15 = qJD(3) * t34 + qJD(5) - t70;
t14 = t80 * t90;
t13 = -t40 * qJD(5) + t42 * t58;
t12 = t40 * t58 + t71;
t7 = t40 * t56 - t42 * t89;
t6 = t40 * t89 + t42 * t56;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t37, 0, t49 * qJD(4) - t4 * t42 + t40 * t65, 0, 0, 0, 0, 0, 0, t36, t37, 0, t51 * qJD(4) - t1 * t40 - t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t41, -t82, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6, t16, -qJD(3) * t88 + ((t28 - t70) * qJD(3) + t46) * t41, 0, 0, 0, 0, 0, 0, t7, t6, t16, -t87 * t43 + (qJD(3) * t15 + t47) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t14, t37, t23, -t36, 0, -pkin(6) * t37 + t40 * t48 + t24, pkin(6) * t36 + t42 * t48 + t25, -qJD(3) * t53 + t46, (t88 + (-t28 - t77) * t41) * qJD(2) + t46 * pkin(6), t22, t14, t37, t23, -t36, 0, t24 + t50 * t42 + (t13 + (t15 + (t34 - t84) * qJD(3)) * t40) * qJD(4), t25 - t50 * t40 + (t15 * t42 - t12 + (t34 * t42 + t85) * qJD(3)) * qJD(4), (t12 * t42 - t13 * t40 + (-t20 * t42 + t21 * t40) * qJD(4) - t53) * qJD(3) + t47, t15 * pkin(4) * t73 - t1 * t21 + t9 * t12 + t5 * t13 + t17 * t34 + t2 * t20 + (-t15 * t41 + t51 * t43) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t19, 0, t29, 0, 0, t32 + (t11 - t83) * qJD(4) + (-t28 - t70) * t75, t10 * qJD(4) - t28 * t74 + t65, 0, 0, -t29, t19, 0, t29, 0, 0, pkin(4) * t29 + t32 + (t9 - t83) * qJD(4) + (-t63 + (-t15 + t55) * t40) * qJD(3), -t45 * t85 + t8 * qJD(4) + (t64 + (-qJD(5) - t15) * t42) * qJD(3) + t65, (-t76 + t86) * t74, t86 * t9 + (-t15 * t75 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t60, 0.2e1 * t59, -t79 * t45, t87;];
tauc_reg = t3;
