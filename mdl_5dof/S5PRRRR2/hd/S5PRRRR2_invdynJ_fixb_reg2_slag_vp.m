% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:24
% EndTime: 2019-07-18 13:30:25
% DurationCPUTime: 0.89s
% Computational Cost: add. (1408->172), mult. (2409->217), div. (0->0), fcn. (1334->12), ass. (0->122)
t127 = pkin(2) * qJD(2);
t70 = sin(qJ(3));
t115 = t70 * t127;
t106 = qJD(4) * t115;
t126 = pkin(3) * qJD(4);
t112 = qJD(3) * t127;
t74 = cos(qJ(3));
t144 = t74 * pkin(2);
t53 = qJDD(2) * t144;
t62 = qJDD(2) + qJDD(3);
t22 = t62 * pkin(3) - t70 * t112 + t53;
t73 = cos(qJ(4));
t19 = t73 * t22;
t63 = qJD(2) + qJD(3);
t69 = sin(qJ(4));
t118 = qJD(3) + qJD(4);
t119 = t70 * qJDD(2);
t124 = qJD(2) * t74;
t83 = (t118 * t124 + t119) * pkin(2);
t7 = t73 * t106 + (t63 * t126 + t83) * t69 - t19;
t67 = qJ(2) + qJ(3);
t60 = qJ(4) + t67;
t50 = cos(t60);
t146 = g(2) * t50;
t113 = -t7 - t146;
t49 = sin(t60);
t150 = -g(1) * t50 - g(2) * t49;
t91 = pkin(2) * t124 + t63 * pkin(3);
t17 = t73 * t115 + t69 * t91;
t56 = qJD(2) + t118;
t13 = t56 * pkin(6) + t17;
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t8 = t72 * qJD(1) - t68 * t13;
t125 = t8 * qJD(5);
t134 = t73 * t74;
t137 = t70 * t73;
t44 = pkin(2) * t137;
t87 = t73 * t91;
t104 = -qJD(4) * t87 - qJDD(2) * t44 - t112 * t134 + (t106 - t22) * t69;
t55 = qJDD(4) + t62;
t4 = t55 * pkin(6) - t104;
t2 = t68 * qJDD(1) + t72 * t4 + t125;
t57 = t72 * qJDD(1);
t136 = t72 * t13;
t9 = t68 * qJD(1) + t136;
t3 = -t9 * qJD(5) - t68 * t4 + t57;
t80 = t2 * t72 - t3 * t68 + (-t68 * t9 - t72 * t8) * qJD(5);
t79 = t150 + t80;
t64 = t68 ^ 2;
t65 = t72 ^ 2;
t128 = t64 + t65;
t110 = t128 * t56;
t42 = g(1) * t49;
t149 = t42 - t146;
t123 = qJD(4) * t69;
t94 = t69 * t74 + t137;
t24 = t94 * t127;
t139 = t24 * t56;
t51 = t69 * pkin(3) + pkin(6);
t76 = qJD(5) ^ 2;
t148 = -(t56 * t123 - t73 * t55) * pkin(3) - t51 * t76 + t139;
t58 = sin(t67);
t46 = g(1) * t58;
t59 = cos(t67);
t145 = g(2) * t59;
t16 = t69 * t115 - t87;
t143 = t16 * t24;
t142 = t16 * t56;
t141 = t16 * t69;
t140 = t17 * t56;
t138 = t69 * t70;
t135 = t72 * t55;
t121 = qJD(5) * t16;
t133 = t68 * t121 + t72 * t42;
t52 = pkin(3) + t144;
t27 = t69 * t52 + t44;
t131 = g(1) * t59 + g(2) * t58;
t48 = pkin(3) * t59;
t75 = cos(qJ(2));
t130 = t75 * pkin(2) + t48;
t129 = t64 - t65;
t122 = qJD(4) * t73;
t120 = qJD(5) * t72;
t66 = qJDD(1) - g(3);
t117 = -t113 * t68 + t16 * t120;
t54 = t56 ^ 2;
t116 = t68 * t54 * t72;
t111 = t128 * t55;
t109 = qJD(2) * (-qJD(3) + t63);
t108 = qJD(3) * (-qJD(2) - t63);
t107 = t68 * t56 * t120;
t105 = t46 + t53 - t145;
t71 = sin(qJ(2));
t103 = -t71 * pkin(2) - pkin(3) * t58;
t101 = g(1) * t71 - g(2) * t75;
t100 = -t7 * t73 + t46;
t99 = t8 * t68 - t9 * t72;
t98 = pkin(6) * t76 - t140;
t11 = t52 * t123 + (t94 * qJD(3) + t70 * t122) * pkin(2);
t26 = pkin(2) * t138 - t73 * t52;
t97 = t16 * t11 + t7 * t26;
t96 = -t11 * t56 - t26 * t55;
t93 = t134 - t138;
t92 = -pkin(6) * qJDD(5) - t121;
t89 = t104 - t150;
t23 = pkin(6) + t27;
t86 = t23 * t76 - t96;
t10 = t52 * t122 + (t93 * qJD(3) - t70 * t123) * pkin(2);
t85 = -qJDD(5) * t23 + (t26 * t56 - t10) * qJD(5);
t82 = -qJD(1) * qJD(5) - t142 - t150 - t4;
t25 = t93 * t127;
t81 = -qJDD(5) * t51 + (t25 + (-qJD(4) - t56) * t73 * pkin(3)) * qJD(5);
t78 = t149 - t7;
t40 = t50 * pkin(6);
t39 = t49 * pkin(6);
t31 = qJDD(5) * t72 - t76 * t68;
t30 = qJDD(5) * t68 + t76 * t72;
t21 = t65 * t55 - 0.2e1 * t107;
t20 = t64 * t55 + 0.2e1 * t107;
t12 = -0.2e1 * t129 * t56 * qJD(5) + 0.2e1 * t68 * t135;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, t31, -t30, 0, -t99 * qJD(5) + t2 * t68 + t3 * t72 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t101, g(1) * t75 + g(2) * t71, 0, 0, 0, 0, 0, 0, 0, t62, (t70 * t108 + t62 * t74) * pkin(2) + t105, ((-qJDD(2) - t62) * t70 + t74 * t108) * pkin(2) + t131, 0, (t101 + (t70 ^ 2 + t74 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), 0, 0, 0, 0, 0, t55, t78 + t96, -t10 * t56 - t27 * t55 + t89, 0, -g(1) * t103 - g(2) * t130 + t17 * t10 - t104 * t27 + t97, t20, t12, t30, t21, t31, 0, t85 * t68 + (t113 - t86) * t72 + t133, t85 * t72 + (t86 - t42) * t68 + t117, t10 * t110 + t23 * t111 + t79, -g(1) * (t103 + t40) - g(2) * (t39 + t130) - t99 * t10 + t80 * t23 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t70 * pkin(2) * t109 + t105, (t74 * t109 - t119) * pkin(2) + t131, 0, 0, 0, 0, 0, 0, 0, t55, t139 + t19 + (pkin(3) * t55 - t106) * t73 + ((-t56 - t63) * t126 - t83) * t69 + t149, t25 * t56 + (-t56 * t122 - t55 * t69) * pkin(3) + t89, 0, -t143 - t17 * t25 + (-t145 - t104 * t69 + (t17 * t73 + t141) * qJD(4) + t100) * pkin(3), t20, t12, t30, t21, t31, 0, t81 * t68 + (t113 + t148) * t72 + t133, t81 * t72 + (-t148 - t42) * t68 + t117, t51 * t111 + t79 + (t122 * pkin(3) - t25) * t110, -t143 - g(1) * t40 - g(2) * (t39 + t48) + t99 * t25 + t80 * t51 + ((-t99 * t73 + t141) * qJD(4) + t100) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t78 + t140, t89 - t142, 0, 0, t20, t12, t30, t21, t31, 0, t92 * t68 + (t113 - t98) * t72 + t133, t92 * t72 + (t98 - t42) * t68 + t117, pkin(6) * t111 + t16 * t110 + t79, (-t17 - t99) * t16 + t79 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t129 * t54, t68 * t55, t116, t135, qJDD(5), -g(3) * t72 + t57 + (t9 - t136) * qJD(5) + t82 * t68, t125 + (qJD(5) * t13 - t66) * t68 + t82 * t72, 0, 0;];
tau_reg  = t1;
