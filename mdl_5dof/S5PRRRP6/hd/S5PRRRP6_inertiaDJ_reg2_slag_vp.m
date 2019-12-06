% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:18
% EndTime: 2019-12-05 16:52:21
% DurationCPUTime: 0.74s
% Computational Cost: add. (702->117), mult. (1937->191), div. (0->0), fcn. (1673->6), ass. (0->67)
t84 = cos(qJ(4));
t62 = t84 * qJD(4);
t87 = t84 * qJD(3) + t62;
t86 = qJD(3) + qJD(4);
t50 = 2 * qJD(5);
t85 = -pkin(7) - pkin(6);
t45 = sin(qJ(4));
t46 = sin(qJ(3));
t83 = t45 * t46;
t48 = cos(qJ(3));
t82 = t45 * t48;
t43 = t46 ^ 2;
t44 = t48 ^ 2;
t81 = t43 + t44;
t80 = qJD(4) * t45;
t79 = t46 * qJD(3);
t47 = sin(qJ(2));
t78 = t47 * qJD(2);
t77 = t48 * qJD(3);
t49 = cos(qJ(2));
t76 = t49 * qJD(2);
t28 = t84 * t46 + t82;
t18 = t86 * t28;
t55 = t84 * t48 - t83;
t75 = -0.2e1 * t55 * t18;
t74 = -0.2e1 * pkin(2) * qJD(3);
t73 = pkin(3) * t80;
t72 = pkin(3) * t79;
t71 = t46 * t77;
t70 = t47 * t79;
t69 = t47 * t76;
t68 = t46 * t76;
t67 = t49 * t79;
t41 = -t48 * pkin(3) - pkin(2);
t33 = t85 * t48;
t61 = t85 * t84;
t57 = qJD(3) * t61;
t58 = t46 * t61;
t65 = t85 * qJD(3);
t10 = -qJD(4) * t58 - t33 * t80 - t46 * t57 - t65 * t82;
t11 = -t33 * t62 - t48 * t57 + (qJD(4) * t85 + t65) * t83;
t19 = -t45 * t33 - t58;
t20 = -t84 * t33 + t85 * t83;
t66 = -t20 * t10 + t19 * t11;
t64 = t81 * t49;
t60 = t84 * t76;
t17 = -t87 * t48 + t86 * t83;
t59 = -t17 * t55 - t28 * t18;
t56 = t49 * t17 + t28 * t78;
t21 = t28 * t47;
t22 = t55 * t47;
t8 = t18 * t47 + t45 * t68 - t48 * t60;
t9 = -t45 * t70 + (t45 * t76 + t87 * t47) * t48 + (-t47 * t80 + t60) * t46;
t54 = -t22 * t10 + t21 * t11 + t9 * t19 - t8 * t20;
t53 = -t21 * t17 - t22 * t18 + t9 * t28 - t55 * t8;
t52 = -0.2e1 * t10 * t55 + 0.2e1 * t11 * t28 - 0.2e1 * t19 * t17 - 0.2e1 * t20 * t18;
t51 = 0.2e1 * t21 * t9 - 0.2e1 * t22 * t8 - 0.2e1 * t69;
t42 = pkin(3) * t62;
t40 = -t84 * pkin(3) - pkin(4);
t38 = t45 * pkin(3) + qJ(5);
t37 = -0.2e1 * t73;
t34 = t42 + qJD(5);
t16 = -pkin(4) * t55 - t28 * qJ(5) + t41;
t15 = -0.2e1 * t28 * t17;
t13 = -t49 * t18 - t55 * t78;
t3 = t18 * pkin(4) + t17 * qJ(5) - t28 * qJD(5) + t72;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t81) * t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t76, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t78 - t67, t46 * t78 - t49 * t77, qJD(2) * t64, (-pkin(2) * t47 + pkin(6) * t64) * qJD(2), 0, 0, 0, 0, 0, 0, t13, t56, t53, -pkin(3) * t67 + t41 * t78 + t54, 0, 0, 0, 0, 0, 0, t13, t53, -t56, t16 * t78 - t49 * t3 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t71, 0.2e1 * (-t43 + t44) * qJD(3), 0, -0.2e1 * t71, 0, 0, t46 * t74, t48 * t74, 0, 0, t15, 0.2e1 * t59, 0, t75, 0, 0, 0.2e1 * t41 * t18 - 0.2e1 * t55 * t72, -0.2e1 * t41 * t17 + 0.2e1 * t28 * t72, t52, 0.2e1 * t41 * t72 + 0.2e1 * t66, t15, 0, -0.2e1 * t59, 0, 0, t75, 0.2e1 * t16 * t18 - 0.2e1 * t3 * t55, t52, 0.2e1 * t16 * t17 - 0.2e1 * t3 * t28, 0.2e1 * t16 * t3 + 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t77 - t68, -t48 * t76 + t70, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, 0, (-t84 * t9 - t45 * t8 + (t21 * t45 + t84 * t22) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, -t9, 0, -t8, t21 * t73 + t22 * t34 - t8 * t38 + t9 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t79, 0, -pkin(6) * t77, pkin(6) * t79, 0, 0, 0, 0, -t17, 0, -t18, 0, -t11, t10, (t84 * t17 - t18 * t45 + (t28 * t45 + t55 * t84) * qJD(4)) * pkin(3), (-t84 * t11 - t10 * t45 + (t19 * t45 + t84 * t20) * qJD(4)) * pkin(3), 0, -t17, 0, 0, t18, 0, -t11, -t40 * t17 - t38 * t18 + t28 * t73 + t34 * t55, -t10, -t10 * t38 + t11 * t40 + t19 * t73 + t20 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -0.2e1 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0.2e1 * t34, 0.2e1 * t38 * t34 + 0.2e1 * t40 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, -t8, -t9 * pkin(4) - t8 * qJ(5) + t22 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t18, 0, -t11, t10, 0, 0, 0, -t17, 0, 0, t18, 0, -t11, pkin(4) * t17 - t18 * qJ(5) + qJD(5) * t55, -t10, -t11 * pkin(4) - t10 * qJ(5) + t20 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t42, 0, 0, 0, 0, 0, 0, 0, 0, -t73, 0, t50 + t42, -pkin(4) * t73 + t34 * qJ(5) + t38 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, qJ(5) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
