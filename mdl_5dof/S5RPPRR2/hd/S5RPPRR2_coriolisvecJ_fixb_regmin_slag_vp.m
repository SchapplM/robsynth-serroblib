% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:40:01
% DurationCPUTime: 0.82s
% Computational Cost: add. (811->127), mult. (1944->184), div. (0->0), fcn. (1454->6), ass. (0->89)
t68 = cos(pkin(8));
t73 = cos(qJ(4));
t101 = t73 * t68;
t67 = sin(pkin(8));
t71 = sin(qJ(4));
t48 = -t71 * t67 + t101;
t72 = cos(qJ(5));
t47 = t73 * t67 + t71 * t68;
t99 = qJD(1) * t47;
t102 = t72 * t99;
t92 = qJD(1) * t101;
t98 = qJD(1) * t67;
t93 = t71 * t98;
t44 = t92 - t93;
t70 = sin(qJ(5));
t17 = t70 * t44 + t102;
t64 = qJD(4) + qJD(5);
t106 = t17 * t64;
t45 = t47 * qJD(4);
t36 = qJD(1) * t45;
t53 = qJD(4) * t93;
t37 = qJD(4) * t92 - t53;
t97 = qJD(5) * t70;
t2 = -qJD(5) * t102 - t72 * t36 - t70 * t37 - t44 * t97;
t123 = t2 + t106;
t82 = t72 * t44 - t70 * t99;
t122 = t82 * t17;
t121 = t48 * qJD(3);
t107 = t82 * t64;
t3 = qJD(5) * t82 - t70 * t36 + t72 * t37;
t120 = -t3 + t107;
t119 = -t17 ^ 2 + t82 ^ 2;
t69 = -pkin(1) - qJ(3);
t115 = t69 * qJD(1);
t54 = qJD(2) + t115;
t86 = -pkin(6) * qJD(1) + t54;
t38 = t86 * t67;
t39 = t86 * t68;
t83 = -t73 * t38 - t71 * t39;
t13 = -pkin(7) * t99 - t83;
t66 = qJD(1) * qJ(2);
t60 = qJD(3) + t66;
t51 = pkin(3) * t98 + t60;
t23 = pkin(4) * t99 + t51;
t80 = t121 * qJD(1);
t5 = t36 * pkin(7) + qJD(4) * t83 - t80;
t118 = t23 * t17 + t13 * t97 + (-t13 * t64 - t5) * t70;
t116 = -t71 * t38 + t73 * t39;
t114 = qJD(5) - t64;
t100 = t67 ^ 2 + t68 ^ 2;
t113 = t100 * qJD(3);
t78 = t47 * qJD(3);
t4 = -t37 * pkin(7) - qJD(1) * t78 + t116 * qJD(4);
t112 = -t23 * t82 - t70 * t4 + t72 * t5;
t65 = qJD(1) * qJD(2);
t90 = 0.2e1 * t65;
t111 = pkin(4) * t44;
t21 = t72 * t47 + t70 * t48;
t94 = t73 * qJD(4);
t95 = t71 * qJD(4);
t46 = -t67 * t95 + t68 * t94;
t6 = -qJD(5) * t21 - t72 * t45 - t70 * t46;
t110 = t6 * t64;
t22 = -t70 * t47 + t72 * t48;
t7 = qJD(5) * t22 - t70 * t45 + t72 * t46;
t109 = t7 * t64;
t108 = -pkin(6) + t69;
t103 = t72 * t13;
t58 = t67 * pkin(3) + qJ(2);
t96 = t46 * qJD(4);
t12 = -t44 * pkin(7) + t116;
t11 = qJD(4) * pkin(4) + t12;
t88 = -pkin(4) * t64 - t11;
t85 = qJD(1) * t100;
t49 = t108 * t67;
t50 = t108 * t68;
t81 = -t73 * t49 - t71 * t50;
t76 = -t49 * t95 + t50 * t94 - t78;
t75 = qJD(4) * t81 - t121;
t74 = qJD(1) ^ 2;
t40 = t45 * qJD(4);
t31 = t46 * pkin(4) + qJD(2);
t29 = t47 * pkin(4) + t58;
t24 = t37 * pkin(4) + t65;
t15 = -t47 * pkin(7) - t81;
t14 = -t48 * pkin(7) - t71 * t49 + t73 * t50;
t10 = t45 * pkin(7) + t75;
t9 = -t46 * pkin(7) + t76;
t1 = [0, 0, 0, 0, t90, qJ(2) * t90, t67 * t90, t68 * t90, 0.2e1 * qJD(3) * t85, (t60 + t66) * qJD(2) + (-t54 - t115) * t113, -t36 * t48 - t44 * t45, t36 * t47 - t48 * t37 - t44 * t46 + t45 * t99, -t40, -t96, 0, 0.2e1 * t99 * qJD(2) + t75 * qJD(4) + t58 * t37 + t51 * t46, -t58 * t36 - t51 * t45 - t76 * qJD(4) + (qJD(1) * t48 + t44) * qJD(2), t2 * t22 + t6 * t82, -t6 * t17 - t2 * t21 - t22 * t3 - t7 * t82, t110, -t109, 0, t31 * t17 + t29 * t3 + t24 * t21 + t23 * t7 + (t72 * t10 - t70 * t9 + (-t14 * t70 - t15 * t72) * qJD(5)) * t64, t31 * t82 + t29 * t2 + t24 * t22 + t23 * t6 - (t70 * t10 + t72 * t9 + (t14 * t72 - t15 * t70) * qJD(5)) * t64; 0, 0, 0, 0, -t74, -t74 * qJ(2), -t74 * t67, -t74 * t68, 0, (-t60 - t113) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t99 - t40, -qJD(1) * t44 - t96, 0, 0, 0, 0, 0, -qJD(1) * t17 + t110, -qJD(1) * t82 - t109; 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t74, t54 * t85 + t65, 0, 0, 0, 0, 0, -t53 + (t44 + t92) * qJD(4), -0.2e1 * t99 * qJD(4), 0, 0, 0, 0, 0, t3 + t107, t2 - t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t99, t44 ^ 2 - t99 ^ 2, 0, t53 + (t44 - t92) * qJD(4), 0, -t51 * t44 - t80, (qJD(3) + t51) * t99, t122, t119, t123, t120, 0, -t17 * t111 - (-t70 * t12 - t103) * t64 + (t88 * t70 - t103) * qJD(5) + t112, -t82 * t111 + (t88 * qJD(5) + t12 * t64 - t4) * t72 + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t119, t123, t120, 0, t114 * (-t70 * t11 - t103) + t112, (-t114 * t11 - t4) * t72 + t118;];
tauc_reg = t1;
