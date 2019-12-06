% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:17
% DurationCPUTime: 0.64s
% Computational Cost: add. (567->118), mult. (1564->165), div. (0->0), fcn. (1106->6), ass. (0->84)
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t99 = -t52 * t49 + t54 * t50;
t26 = t99 * qJD(1);
t32 = t54 * t49 + t52 * t50;
t29 = t32 * qJD(3);
t28 = t99 * qJD(3);
t22 = qJD(1) * t28;
t100 = -qJD(2) * qJD(4) - t22;
t53 = cos(qJ(4));
t27 = t32 * qJD(1);
t21 = qJD(3) * pkin(6) + t27;
t64 = qJ(5) * qJD(3) + t21;
t59 = t64 * t53;
t51 = sin(qJ(4));
t77 = t51 * qJD(2);
t10 = t59 + t77;
t83 = qJD(4) * pkin(4);
t45 = t53 * qJD(2);
t9 = -t51 * t64 + t45;
t8 = t9 + t83;
t98 = t8 - t9;
t47 = t51 ^ 2;
t97 = pkin(4) * t47;
t96 = t53 * pkin(4);
t23 = qJD(1) * t29;
t95 = t23 * t99;
t56 = qJD(3) ^ 2;
t92 = t53 * t56;
t55 = qJD(4) ^ 2;
t89 = t55 * t51;
t46 = t55 * t53;
t88 = -qJ(5) - pkin(6);
t79 = t27 * qJD(3);
t81 = qJD(4) * t51;
t87 = t26 * t81 + t53 * t79;
t48 = t53 ^ 2;
t86 = t47 - t48;
t85 = t47 + t48;
t84 = qJD(3) * pkin(3);
t82 = qJD(3) * t53;
t80 = qJD(4) * t53;
t78 = t28 * qJD(3);
t76 = t53 * qJD(5);
t44 = -pkin(3) - t96;
t13 = qJD(3) * t44 + qJD(5) - t26;
t75 = -qJD(5) - t13;
t72 = qJD(3) * qJD(4);
t71 = t100 * t53 + t21 * t81;
t67 = t51 * t72;
t14 = pkin(4) * t67 + t23;
t70 = qJ(5) * t81;
t69 = -pkin(6) * t55 - t23;
t68 = t85 * t26;
t66 = t53 * t72;
t65 = qJD(4) * t88;
t63 = t51 * t66;
t62 = t10 * t53 - t51 * t8;
t11 = -t51 * t21 + t45;
t12 = t53 * t21 + t77;
t61 = t11 * t51 - t12 * t53;
t20 = -t26 - t84;
t60 = qJD(4) * (t20 - t84);
t3 = (-t70 + t76) * qJD(3) - t71;
t4 = (-qJD(3) * qJD(5) - t22) * t51 - t10 * qJD(4);
t58 = t3 * t53 - t4 * t51 + (-t10 * t51 - t53 * t8) * qJD(4);
t6 = -t12 * qJD(4) - t51 * t22;
t57 = -t71 * t53 - t6 * t51 + (-t11 * t53 - t12 * t51) * qJD(4);
t41 = t51 * t92;
t37 = -0.2e1 * t63;
t36 = 0.2e1 * t63;
t35 = t88 * t53;
t34 = t88 * t51;
t33 = t86 * t56;
t30 = -0.2e1 * t86 * t72;
t25 = -t51 * qJD(5) + t53 * t65;
t24 = t51 * t65 + t76;
t18 = t26 * t80;
t7 = t85 * t78;
t2 = -t28 * t81 - t32 * t46 + (-t29 * t53 - t81 * t99) * qJD(3);
t1 = -t28 * t80 + t32 * t89 + (t29 * t51 - t80 * t99) * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29 * qJD(3), -t78, 0, t22 * t32 - t26 * t29 + t27 * t28 - t95, 0, 0, 0, 0, 0, 0, t2, t1, t7, t20 * t29 - t28 * t61 + t32 * t57 - t95, 0, 0, 0, 0, 0, 0, t2, t1, t7, t13 * t29 - t14 * t99 + t28 * t62 + t32 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t46, 0, -qJD(4) * t61 - t51 * t71 + t6 * t53, 0, 0, 0, 0, 0, 0, -t89, -t46, 0, qJD(4) * t62 + t3 * t51 + t4 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 - t23, 0, 0, 0, t36, t30, t46, t37, -t89, 0, t51 * t60 + t53 * t69 + t87, t18 + t53 * t60 + (-t69 - t79) * t51, -qJD(3) * t68 + t57, -t23 * pkin(3) + pkin(6) * t57 - t20 * t27 + t26 * t61, t36, t30, t46, t37, -t89, 0, -t14 * t53 + (t25 + (t13 + (t44 - t96) * qJD(3)) * t51) * qJD(4) + t87, t18 + (t14 - t79) * t51 + (t13 * t53 - t24 + (t44 * t53 + t97) * qJD(3)) * qJD(4), (t24 * t53 - t25 * t51 - t68 + (-t34 * t53 + t35 * t51) * qJD(4)) * qJD(3) + t58, t14 * t44 - t3 * t35 + t4 * t34 + (t26 * t51 + t25) * t8 + (pkin(4) * t81 - t27) * t13 + (-t26 * t53 + t24) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t33, 0, t41, 0, 0, (-qJD(3) * t20 - t22) * t51, t11 * qJD(4) - t20 * t82 + t71, 0, 0, -t41, t33, 0, t41, 0, 0, (t10 - t59) * qJD(4) + (pkin(4) * t92 + qJD(3) * t75 + t100) * t51, -t56 * t97 + t9 * qJD(4) + (t53 * t75 + t70) * qJD(3) + t71, (-t83 + t98) * t82, t98 * t10 + (-t13 * t51 * qJD(3) + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, 0.2e1 * t66, -t85 * t56, -qJD(3) * t62 + t14;];
tauc_reg = t5;
