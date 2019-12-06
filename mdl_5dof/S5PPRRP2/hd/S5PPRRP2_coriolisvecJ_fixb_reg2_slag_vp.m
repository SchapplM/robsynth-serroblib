% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:12
% DurationCPUTime: 0.53s
% Computational Cost: add. (559->103), mult. (1526->137), div. (0->0), fcn. (1089->6), ass. (0->77)
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t87 = -t42 * t39 + t44 * t40;
t21 = t87 * qJD(1);
t43 = cos(qJ(4));
t63 = t43 * qJD(2);
t26 = t44 * t39 + t42 * t40;
t22 = t26 * qJD(1);
t17 = qJD(3) * pkin(6) + t22;
t41 = sin(qJ(4));
t81 = t41 * t17;
t11 = t63 - t81;
t86 = qJD(5) - t11;
t8 = -qJD(4) * pkin(4) + t86;
t12 = t41 * qJD(2) + t43 * t17;
t9 = qJD(4) * qJ(5) + t12;
t45 = qJD(4) ^ 2;
t85 = pkin(6) * t45;
t23 = t87 * qJD(3);
t18 = qJD(1) * t23;
t80 = t41 * t18;
t5 = t12 * qJD(4) + t80;
t84 = t5 * t41;
t83 = t5 * t43;
t24 = t26 * qJD(3);
t19 = qJD(1) * t24;
t82 = t19 * t87;
t79 = t41 * t43;
t76 = t45 * t41;
t36 = t45 * t43;
t75 = qJD(4) * t63 + t43 * t18;
t65 = t22 * qJD(3);
t67 = qJD(4) * t41;
t74 = t21 * t67 + t43 * t65;
t51 = pkin(4) * t41 - qJ(5) * t43;
t20 = t51 * qJD(4) - t41 * qJD(5);
t73 = t20 - t22;
t37 = t41 ^ 2;
t38 = t43 ^ 2;
t72 = -t37 + t38;
t71 = t37 + t38;
t70 = qJD(3) * pkin(3);
t16 = -t21 - t70;
t69 = qJD(3) * t16;
t29 = -t43 * pkin(4) - t41 * qJ(5) - pkin(3);
t68 = qJD(3) * t29;
t66 = qJD(4) * t43;
t64 = t23 * qJD(3);
t61 = qJD(3) * qJD(4);
t7 = (t22 + t20) * qJD(3);
t60 = -t7 - t85;
t59 = -t19 - t85;
t58 = t11 + t81;
t57 = t16 - t70;
t10 = -t21 + t68;
t56 = t10 + t68;
t54 = t61 * t79;
t53 = t71 * t21 * qJD(3);
t52 = t41 * t8 + t43 * t9;
t50 = t11 * t41 - t12 * t43;
t3 = (qJD(5) - t81) * qJD(4) + t75;
t48 = t3 * t43 + t84 + (-t41 * t9 + t43 * t8) * qJD(4);
t4 = -t17 * t67 + t75;
t47 = t4 * t43 + t84 + (-t11 * t43 - t12 * t41) * qJD(4);
t46 = qJD(3) ^ 2;
t34 = t46 * t79;
t32 = -0.2e1 * t54;
t31 = 0.2e1 * t54;
t30 = t72 * t46;
t28 = t51 * qJD(3);
t27 = t72 * t61;
t6 = t71 * t64;
t2 = -t23 * t67 - t26 * t36 + (-t24 * t43 - t67 * t87) * qJD(3);
t1 = -t23 * t66 + t26 * t76 + (t24 * t41 - t66 * t87) * qJD(3);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * qJD(3), -t64, 0, t18 * t26 - t21 * t24 + t22 * t23 - t82, 0, 0, 0, 0, 0, 0, t2, t1, t6, t16 * t24 - t50 * t23 + t47 * t26 - t82, 0, 0, 0, 0, 0, 0, t2, t6, -t1, t10 * t24 + t52 * t23 + t48 * t26 - t7 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t36, 0, -t50 * qJD(4) + t4 * t41 - t83, 0, 0, 0, 0, 0, 0, -t76, 0, t36, t52 * qJD(4) + t3 * t41 - t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0.2e1 * t27, t36, t32, -t76, 0, t59 * t43 + t57 * t67 + t74, (-t59 - t65) * t41 + (t21 + t57) * t66, -t53 + t47, -t19 * pkin(3) + t47 * pkin(6) - t16 * t22 + t50 * t21, t31, t36, -0.2e1 * t27, 0, t76, t32, t56 * t67 + (-qJD(3) * t20 + t60) * t43 + t74, -t53 + t48, (-t21 - t56) * t66 + (-t73 * qJD(3) + t60) * t41, t48 * pkin(6) + t73 * t10 - t52 * t21 + t7 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t30, 0, t34, 0, 0, (-t18 - t69) * t41, t58 * qJD(4) - t43 * t69 - t75, 0, 0, -t34, 0, t30, 0, 0, t34, -t80 + (-t10 * t41 + t28 * t43) * qJD(3), 0, (t10 * t43 + t28 * t41) * qJD(3) + (0.2e1 * qJD(5) - t58) * qJD(4) + t75, -t5 * pkin(4) + t3 * qJ(5) - t10 * t28 - t8 * t12 + t86 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t37 * t46 - t45, (qJD(3) * t10 + t18) * t41 + (t12 - t9) * qJD(4);];
tauc_reg = t13;
