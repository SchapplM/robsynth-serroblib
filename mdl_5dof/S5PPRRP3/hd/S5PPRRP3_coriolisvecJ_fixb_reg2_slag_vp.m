% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRP3
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:19
% DurationCPUTime: 0.64s
% Computational Cost: add. (580->122), mult. (1572->170), div. (0->0), fcn. (1084->6), ass. (0->94)
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t44 = sin(pkin(8));
t87 = qJD(1) * t44;
t28 = t47 * qJD(2) + t49 * t87;
t21 = qJD(3) * pkin(6) + t28;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t45 = cos(pkin(8));
t86 = qJD(1) * t45;
t59 = t46 * t21 + t48 * t86;
t107 = qJD(5) + t59;
t7 = -qJD(4) * pkin(4) + t107;
t82 = qJD(3) * t48;
t105 = t49 * qJD(2) - t47 * t87;
t23 = t105 * qJD(3);
t50 = qJD(4) ^ 2;
t51 = qJD(3) ^ 2;
t106 = (t50 + t51) * t47;
t38 = t46 * t86;
t97 = t48 * t21;
t12 = -t38 + t97;
t8 = qJD(4) * qJ(5) + t12;
t104 = pkin(6) * t50;
t99 = t44 * t49;
t24 = t45 * t48 + t46 * t99;
t37 = qJD(4) * t38;
t79 = qJD(4) * t48;
t98 = t46 * t23;
t4 = t21 * t79 - t37 + t98;
t103 = t4 * t24;
t102 = t4 * t46;
t22 = t28 * qJD(3);
t101 = t22 * t47;
t100 = t44 * t47;
t96 = t50 * t46;
t95 = t51 * t47;
t94 = t51 * t49;
t80 = qJD(4) * t46;
t93 = t105 * t80 + t22 * t48;
t61 = pkin(4) * t46 - qJ(5) * t48;
t19 = qJD(4) * t61 - t46 * qJD(5);
t92 = t19 - t28;
t42 = t46 ^ 2;
t43 = t48 ^ 2;
t91 = -t42 + t43;
t90 = t42 + t43;
t88 = qJD(3) * pkin(3);
t32 = -t48 * pkin(4) - t46 * qJ(5) - pkin(3);
t83 = qJD(3) * t32;
t13 = -t105 + t83;
t85 = qJD(3) * t13;
t20 = -t105 - t88;
t84 = qJD(3) * t20;
t81 = qJD(3) * t49;
t76 = qJD(3) * qJD(4);
t75 = t48 * t99;
t74 = t44 * t94;
t73 = qJD(3) * t100;
t5 = (t19 + t28) * qJD(3);
t72 = -t5 - t104;
t71 = -t22 - t104;
t70 = t46 * t76;
t69 = t20 - t88;
t68 = t13 + t83;
t66 = t48 * t73;
t65 = -0.2e1 * t49 * t76;
t64 = t48 * t70;
t63 = t90 * t23;
t62 = t46 * t7 + t48 * t8;
t60 = -t12 * t48 - t46 * t59;
t57 = t37 + (t12 - t97) * qJD(4);
t9 = -qJD(4) * t75 + (qJD(4) * t45 + t73) * t46;
t55 = t9 * qJD(4) + t70 * t100 - t48 * t74;
t16 = t48 * t23;
t2 = t16 + (qJD(5) - t59) * qJD(4);
t54 = t2 * t48 + t102 + (-t46 * t8 + t48 * t7) * qJD(4);
t3 = -qJD(4) * t59 + t16;
t53 = t3 * t48 + t102 + (-t12 * t46 + t48 * t59) * qJD(4);
t10 = -qJD(4) * t24 - t66;
t25 = -t45 * t46 + t75;
t52 = t10 * t82 + (-t46 * t9 + (t24 * t48 - t25 * t46) * qJD(4)) * qJD(3);
t41 = t50 * t48;
t40 = t46 * t51 * t48;
t36 = -0.2e1 * t64;
t35 = 0.2e1 * t64;
t34 = t91 * t51;
t30 = t61 * qJD(3);
t29 = t91 * t76;
t26 = t90 * t94;
t15 = -t48 * t106 + t46 * t65;
t14 = t46 * t106 + t48 * t65;
t1 = -t46 * t74 + (t10 - t66) * qJD(4);
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t44 * t95, 0, (t101 + t23 * t49 + (-t105 * t49 - t28 * t47) * qJD(3)) * t44, 0, 0, 0, 0, 0, 0, t55, -t1, t52, t12 * t10 - t59 * t9 + t103 + t3 * t25 + (t20 * t81 + t101) * t44, 0, 0, 0, 0, 0, 0, t55, t52, t1, t8 * t10 + t2 * t25 + t103 - t7 * t9 + (t13 * t81 + t47 * t5) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94, 0, -t22 * t49 + t23 * t47 + (-t105 * t47 + t28 * t49) * qJD(3), 0, 0, 0, 0, 0, 0, t15, t14, t26, (-qJD(3) * t60 - t22) * t49 + (t53 + t84) * t47, 0, 0, 0, 0, 0, 0, t15, t26, -t14, (qJD(3) * t62 - t5) * t49 + (t54 + t85) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0.2e1 * t29, t41, t36, -t96, 0, t48 * t71 + t69 * t80 + t93, (-t71 - t22) * t46 + (t105 + t69) * t79, -t63 + t53, -t22 * pkin(3) + pkin(6) * t53 + t105 * t60 - t20 * t28, t35, t41, -0.2e1 * t29, 0, t96, t36, t68 * t80 + (-qJD(3) * t19 + t72) * t48 + t93, -t63 + t54, (-t105 - t68) * t79 + (-t92 * qJD(3) + t72) * t46, t54 * pkin(6) - t105 * t62 + t92 * t13 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t34, 0, t40, 0, 0, (-t23 - t84) * t46 + t57, -t20 * t82 - t16, 0, 0, -t40, 0, t34, 0, 0, t40, -t98 + (-t13 * t46 + t30 * t48) * qJD(3) + t57, 0, t16 + (t13 * t48 + t30 * t46) * qJD(3) + 0.2e1 * qJD(4) * qJD(5), -t4 * pkin(4) + t2 * qJ(5) + t107 * t8 - t7 * t12 - t13 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, -t42 * t51 - t50, -t37 + (t23 + t85) * t46 + (-t8 + t97) * qJD(4);];
tauc_reg = t6;
