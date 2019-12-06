% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:25
% EndTime: 2019-12-05 18:22:29
% DurationCPUTime: 0.80s
% Computational Cost: add. (1065->178), mult. (1719->226), div. (0->0), fcn. (978->12), ass. (0->120)
t81 = -qJ(5) - pkin(7);
t124 = pkin(1) * qJD(2);
t110 = qJD(1) * t124;
t83 = sin(qJ(2));
t118 = qJDD(1) * t83;
t86 = cos(qJ(2));
t148 = pkin(1) * t118 + t86 * t110;
t78 = qJ(1) + qJ(2);
t67 = pkin(8) + t78;
t56 = sin(t67);
t57 = cos(t67);
t147 = g(2) * t57 + g(3) * t56;
t70 = sin(t78);
t71 = cos(t78);
t146 = g(2) * t71 + g(3) * t70;
t99 = -g(2) * t56 + g(3) * t57;
t84 = sin(qJ(1));
t145 = pkin(1) * t84;
t144 = pkin(1) * t86;
t87 = cos(qJ(1));
t143 = pkin(1) * t87;
t142 = pkin(2) * t70;
t141 = pkin(2) * t71;
t80 = cos(pkin(8));
t140 = pkin(2) * t80;
t85 = cos(qJ(4));
t139 = pkin(4) * t85;
t138 = g(1) * t85;
t125 = pkin(1) * qJD(1);
t114 = t86 * t125;
t75 = qJD(1) + qJD(2);
t39 = pkin(2) * t75 + t114;
t115 = t83 * t125;
t52 = t80 * t115;
t79 = sin(pkin(8));
t24 = t79 * t39 + t52;
t107 = -t81 * t75 + t24;
t82 = sin(qJ(4));
t11 = t85 * qJD(3) - t107 * t82;
t123 = qJD(4) * pkin(4);
t8 = t11 + t123;
t134 = -t11 + t8;
t51 = t79 * t115;
t23 = t80 * t39 - t51;
t19 = -pkin(3) * t75 - t23;
t65 = qJDD(1) * t144;
t74 = qJDD(1) + qJDD(2);
t28 = pkin(2) * t74 - t83 * t110 + t65;
t13 = -t148 * t79 + t28 * t80;
t9 = -pkin(3) * t74 - t13;
t133 = t19 * qJD(4) * t85 + t9 * t82;
t132 = t79 * t83;
t131 = t80 * t83;
t130 = t82 * t74;
t129 = t85 * t74;
t64 = pkin(2) + t144;
t128 = pkin(1) * t131 + t79 * t64;
t76 = t82 ^ 2;
t77 = t85 ^ 2;
t127 = -t76 - t77;
t126 = t76 - t77;
t31 = pkin(7) + t128;
t122 = -qJ(5) - t31;
t58 = pkin(2) * t79 + pkin(7);
t121 = -qJ(5) - t58;
t120 = qJD(4) * t82;
t119 = qJDD(3) - g(1);
t117 = t19 * t120 + t147 * t85;
t14 = t148 * t80 + t79 * t28;
t116 = t65 + t146;
t113 = pkin(4) * t120;
t111 = t75 * t120;
t63 = pkin(3) + t139;
t109 = -g(2) * t70 + g(3) * t71;
t108 = -pkin(1) * t132 + t64 * t80;
t106 = qJD(4) * t122;
t105 = qJD(4) * t121;
t104 = qJD(1) * (-qJD(2) + t75);
t103 = qJD(2) * (-qJD(1) - t75);
t10 = pkin(7) * t74 + t14;
t89 = qJ(5) * t74 + qJD(3) * qJD(4) + qJD(5) * t75 + t10;
t97 = t107 * qJD(4);
t3 = (qJDD(3) - t97) * t82 + t89 * t85;
t102 = t3 * t85 - t99;
t30 = -pkin(3) - t108;
t12 = t82 * qJD(3) + t107 * t85;
t98 = t12 * t85 - t8 * t82;
t96 = t56 * t81 - t57 * t63 - t141;
t95 = -t56 * t63 - t57 * t81 - t142;
t33 = (t79 * t86 + t131) * t124;
t88 = qJD(4) ^ 2;
t94 = -t30 * t74 - t31 * t88 - t33 * t75;
t32 = t79 * t114 + t52;
t59 = -pkin(3) - t140;
t93 = t32 * t75 - t58 * t88 - t59 * t74;
t92 = -t19 * t75 - t10 + t99;
t35 = (t80 * t86 - t132) * t124;
t91 = -qJDD(4) * t31 + (t30 * t75 - t35) * qJD(4);
t34 = t80 * t114 - t51;
t90 = -qJDD(4) * t58 + (t59 * t75 + t34) * qJD(4);
t4 = pkin(4) * t111 - t63 * t74 + qJDD(5) - t13;
t73 = t75 ^ 2;
t72 = t85 * qJ(5);
t68 = t85 * qJD(5);
t66 = t85 * qJDD(3);
t41 = qJDD(4) * t85 - t82 * t88;
t40 = qJDD(4) * t82 + t85 * t88;
t38 = t58 * t85 + t72;
t37 = t121 * t82;
t29 = 0.2e1 * t85 * t111 + t74 * t76;
t27 = -t82 * qJD(5) + t85 * t105;
t26 = t82 * t105 + t68;
t22 = t31 * t85 + t72;
t21 = t122 * t82;
t18 = -0.2e1 * t126 * t75 * qJD(4) + 0.2e1 * t82 * t129;
t15 = -t63 * t75 + qJD(5) - t23;
t6 = (-qJD(5) - t35) * t82 + t85 * t106;
t5 = t82 * t106 + t85 * t35 + t68;
t2 = qJDD(4) * pkin(4) - t89 * t82 - t85 * t97 + t66;
t1 = [qJDD(1), g(2) * t87 + g(3) * t84, -g(2) * t84 + g(3) * t87, t74, (t83 * t103 + t74 * t86) * pkin(1) + t116, ((-qJDD(1) - t74) * t83 + t86 * t103) * pkin(1) + t109, t14 * t128 + t24 * t35 + t13 * t108 - t23 * t33 - g(2) * (-t141 - t143) - g(3) * (-t142 - t145), t29, t18, t40, t41, 0, t91 * t82 + (-t9 + t94) * t85 + t117, t91 * t85 + (-t147 - t94) * t82 + t133, (t22 * t74 + t5 * t75 + (-t21 * t75 - t8) * qJD(4)) * t85 + (-t21 * t74 - t6 * t75 - t2 + (-t22 * t75 - t12) * qJD(4)) * t82 + t102, t3 * t22 + t12 * t5 + t2 * t21 + t8 * t6 + t4 * (t30 - t139) + t15 * (t33 + t113) - g(2) * (t96 - t143) - g(3) * (t95 - t145); 0, 0, 0, t74, t83 * pkin(1) * t104 + t116, (t86 * t104 - t118) * pkin(1) + t109, t23 * t32 - t24 * t34 + (t13 * t80 + t14 * t79 + t146) * pkin(2), t29, t18, t40, t41, 0, t90 * t82 + (-t9 + t93) * t85 + t117, t90 * t85 + (-t147 - t93) * t82 + t133, (-qJD(4) * t8 + t38 * t74) * t85 + (-qJD(4) * t12 - t37 * t74 - t2) * t82 + (t26 * t85 - t27 * t82 + t127 * t34 + (-t37 * t85 - t38 * t82) * qJD(4)) * t75 + t102, t3 * t38 + t2 * t37 + t4 * (-t63 - t140) - g(2) * t96 - g(3) * t95 + (t34 * t82 + t27) * t8 + (-t32 + t113) * t15 + (-t34 * t85 + t26) * t12; 0, 0, 0, 0, 0, 0, t119, 0, 0, 0, 0, 0, t41, -t40, 0, t98 * qJD(4) + t2 * t85 + t3 * t82 - g(1); 0, 0, 0, 0, 0, 0, 0, -t82 * t73 * t85, t126 * t73, t130, t129, qJDD(4), t92 * t82 - t138 + t66, -t119 * t82 + t92 * t85, -pkin(4) * t130 + (-t123 + t134) * t85 * t75, t134 * t12 + (-t138 + t2 + (-t15 * t75 + t99) * t82) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t73, -t98 * t75 - t147 + t4;];
tau_reg = t1;
