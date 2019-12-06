% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:22
% EndTime: 2019-12-05 15:41:25
% DurationCPUTime: 0.68s
% Computational Cost: add. (570->159), mult. (1029->192), div. (0->0), fcn. (602->6), ass. (0->104)
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t69 = t47 * pkin(4) - t49 * qJ(5);
t24 = qJ(3) + t69;
t101 = qJD(2) * t24;
t48 = sin(qJ(2));
t97 = t48 * qJD(1);
t9 = t97 + t101;
t103 = t9 * qJD(2);
t118 = pkin(2) + pkin(6);
t50 = cos(qJ(2));
t96 = t50 * qJD(1);
t76 = qJD(3) - t96;
t22 = -t118 * qJD(2) + t76;
t113 = t47 * t22;
t93 = qJD(4) * qJ(5);
t10 = t93 + t113;
t109 = t49 * t22;
t85 = qJD(1) * qJD(2);
t37 = t48 * t85;
t87 = t50 * qJDD(1);
t67 = qJDD(3) + t37 - t87;
t11 = -t118 * qJDD(2) + t67;
t6 = t47 * t11;
t81 = qJDD(4) * qJ(5);
t2 = t81 + t6 + (qJD(5) + t109) * qJD(4);
t99 = qJDD(4) * pkin(4);
t122 = qJD(4) * t113 - t99;
t7 = t49 * t11;
t3 = qJDD(5) - t7 + t122;
t77 = qJD(4) * pkin(4) - qJD(5);
t8 = -t77 - t109;
t57 = t2 * t47 - t3 * t49 + (t10 * t49 + t8 * t47) * qJD(4);
t124 = -t57 + t103;
t70 = pkin(4) * t49 + qJ(5) * t47;
t13 = qJD(4) * t70 - t49 * qJD(5) + qJD(3);
t88 = t48 * qJDD(1);
t92 = qJDD(2) * t24;
t1 = t88 + t92 + (t13 + t96) * qJD(2);
t52 = qJD(4) ^ 2;
t45 = sin(pkin(7));
t46 = cos(pkin(7));
t75 = g(1) * t46 + g(2) * t45;
t55 = -t50 * (t75 + t85) - g(3) * t48 + t118 * t52;
t123 = qJD(2) * t13 + t1 + t55 + t92;
t53 = qJD(2) ^ 2;
t26 = -t50 * qJDD(2) + t53 * t48;
t94 = qJD(2) * qJ(3);
t29 = t94 + t97;
t121 = qJD(4) * (t29 + t94 - t97);
t95 = qJDD(1) - g(3);
t119 = -t48 * t95 + t50 * t75;
t42 = g(3) * t50;
t115 = t45 * t48;
t114 = t46 * t48;
t112 = t47 * t48;
t111 = t47 * t49;
t110 = t48 * t49;
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t107 = t43 - t44;
t106 = t43 + t44;
t105 = t52 + t53;
t104 = qJ(3) * t50;
t100 = qJDD(2) * pkin(2);
t98 = t29 * qJD(2);
t91 = qJDD(2) * t48;
t90 = qJDD(4) * t47;
t89 = qJDD(4) * t118;
t38 = t49 * qJDD(2);
t84 = qJD(2) * qJD(3);
t83 = qJD(2) * qJD(4);
t82 = qJDD(2) * qJ(3);
t80 = -g(1) * t114 - g(2) * t115 + t42;
t79 = g(3) * (t50 * pkin(2) + t48 * qJ(3));
t78 = t106 * qJDD(2);
t71 = t10 * t47 - t8 * t49;
t68 = (-qJD(2) * pkin(2) + t76) * t48 + t29 * t50;
t17 = t46 * t112 + t45 * t49;
t19 = -t45 * t112 + t46 * t49;
t66 = g(1) * t17 - g(2) * t19 - t6;
t65 = -t80 + t87;
t16 = -t46 * t110 + t45 * t47;
t18 = t45 * t110 + t46 * t47;
t64 = g(1) * t16 - g(2) * t18 + t49 * t42 + t7;
t62 = t105 * t50 + t91;
t61 = -qJDD(4) * t50 + 0.2e1 * t48 * t83;
t59 = -qJDD(5) + t64;
t15 = t67 - t100;
t58 = (t9 - t97 + t101) * qJD(4);
t12 = t82 + t88 + (qJD(3) + t96) * qJD(2);
t54 = t12 + t55 + t82 + t84;
t39 = qJDD(4) * t49;
t33 = t53 * t111;
t32 = t49 * t89;
t31 = t46 * t104;
t30 = t45 * t104;
t25 = t53 * t50 + t91;
t23 = t70 * qJD(2);
t21 = -t105 * t47 + t39;
t20 = t105 * t49 + t90;
t5 = -t47 * t61 + t49 * t62;
t4 = t47 * t62 + t49 * t61;
t14 = [t95, 0, -t26, -t25, t26, t25, qJD(2) * t68 + t12 * t48 - t15 * t50 - g(3), 0, 0, 0, 0, 0, t4, t5, t4, -t26 * t106, -t5, -g(3) + (qJD(2) * t71 + t1) * t48 + t124 * t50; 0, qJDD(2), t65, t119, qJDD(3) - t65 - 0.2e1 * t100, -t119 + 0.2e1 * t82 + 0.2e1 * t84, t12 * qJ(3) + t29 * qJD(3) - t15 * pkin(2) - g(1) * (-pkin(2) * t114 + t31) - g(2) * (-pkin(2) * t115 + t30) - t79 - t68 * qJD(1), t44 * qJDD(2) - 0.2e1 * t83 * t111, 0.2e1 * t107 * t83 - 0.2e1 * t47 * t38, -t52 * t47 + t39, -t52 * t49 - t90, 0, t49 * t121 + t47 * t54 - t32, (t89 - t121) * t47 + t54 * t49, t123 * t47 + t49 * t58 - t32, t106 * t37 + t118 * t78 - t57 - t80, (t58 - t89) * t47 - t123 * t49, t1 * t24 + t9 * t13 - g(1) * t31 - g(2) * t30 - t79 + (-g(3) * pkin(6) - t9 * qJD(1) - t75 * t69) * t50 - t57 * t118 + (-g(3) * t69 - t71 * qJD(1) + t75 * t118) * t48; 0, 0, 0, 0, qJDD(2), -t53, t15 + t80 - t98, 0, 0, 0, 0, 0, t21, -t20, t21, -t78, t20, t80 - t124; 0, 0, 0, 0, 0, 0, 0, t33, -t107 * t53, t38, -t47 * qJDD(2), qJDD(4), -t49 * t98 + t64, (t98 - t42) * t47 + t66, 0.2e1 * t99 + (-t23 * t47 - t49 * t9) * qJD(2) + t59, -t70 * qJDD(2) + ((t10 - t93) * t49 + (t77 + t8) * t47) * qJD(2), t47 * t42 + 0.2e1 * t81 + 0.2e1 * qJD(4) * qJD(5) + (t23 * t49 - t47 * t9) * qJD(2) - t66, t2 * qJ(5) - t3 * pkin(4) - t9 * t23 - t8 * t113 - g(1) * (-t16 * pkin(4) + t17 * qJ(5)) - g(2) * (t18 * pkin(4) - t19 * qJ(5)) + t70 * t42 + (qJD(5) - t109) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t33, t38, -t44 * t53 - t52, -t10 * qJD(4) + t103 * t49 + t122 - t59;];
tau_reg = t14;
