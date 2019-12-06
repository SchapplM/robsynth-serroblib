% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:55
% EndTime: 2019-12-05 15:12:57
% DurationCPUTime: 0.92s
% Computational Cost: add. (1510->165), mult. (3165->218), div. (0->0), fcn. (2636->14), ass. (0->113)
t75 = sin(pkin(9));
t77 = cos(pkin(9));
t81 = sin(qJ(3));
t84 = cos(qJ(3));
t148 = -t81 * t75 + t84 * t77;
t37 = t148 * qJD(1);
t71 = pkin(9) + qJ(3);
t67 = qJ(4) + t71;
t60 = sin(t67);
t78 = cos(pkin(8));
t134 = t60 * t78;
t76 = sin(pkin(8));
t135 = t60 * t76;
t150 = g(1) * t134 + g(2) * t135;
t61 = cos(t67);
t142 = g(3) * t61;
t70 = qJDD(3) + qJDD(4);
t141 = t70 * pkin(4);
t106 = t148 * qJDD(1);
t42 = t84 * t75 + t81 * t77;
t40 = t42 * qJD(3);
t88 = -qJD(1) * t40 + t106;
t24 = qJDD(3) * pkin(3) + t88;
t146 = t42 * qJDD(1);
t39 = t148 * qJD(3);
t25 = qJD(1) * t39 + t146;
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t111 = -t83 * t24 + t80 * t25;
t38 = t42 * qJD(1);
t128 = t83 * t38;
t33 = qJD(3) * pkin(3) + t37;
t18 = t80 * t33 + t128;
t8 = -t18 * qJD(4) - t111;
t6 = -t141 - t8;
t149 = t6 + t142;
t72 = qJD(3) + qJD(4);
t16 = t72 * pkin(7) + t18;
t79 = sin(qJ(5));
t82 = cos(qJ(5));
t11 = t82 * qJD(2) - t79 * t16;
t12 = t79 * qJD(2) + t82 * t16;
t121 = t11 * qJD(5);
t123 = qJD(4) * t83;
t124 = qJD(4) * t80;
t110 = -t33 * t123 + t38 * t124 - t80 * t24 - t83 * t25;
t5 = t70 * pkin(7) - t110;
t2 = t79 * qJDD(2) + t82 * t5 + t121;
t120 = t12 * qJD(5);
t66 = t82 * qJDD(2);
t3 = -t79 * t5 - t120 + t66;
t91 = t2 * t82 - t3 * t79 + (-t11 * t82 - t12 * t79) * qJD(5);
t131 = t80 * t38;
t22 = t83 * t37 - t131;
t147 = pkin(3) * t123 - t22;
t113 = -g(1) * t76 + g(2) * t78;
t145 = -t142 + t150;
t21 = t80 * t37 + t128;
t62 = t80 * pkin(3) + pkin(7);
t139 = t83 * pkin(3);
t63 = -pkin(4) - t139;
t85 = qJD(5) ^ 2;
t144 = (pkin(3) * t124 - t21) * t72 + t62 * t85 + t63 * t70;
t59 = g(3) * t60;
t140 = t72 * pkin(4);
t103 = t148 * t83 - t80 * t42;
t9 = t103 * qJD(4) + t83 * t39 - t80 * t40;
t138 = t9 * t72;
t17 = t83 * t33 - t131;
t137 = t17 * t72;
t136 = t18 * t72;
t133 = t61 * t76;
t132 = t61 * t78;
t129 = t82 * t70;
t73 = t79 ^ 2;
t74 = t82 ^ 2;
t126 = t73 - t74;
t125 = t73 + t74;
t122 = qJD(5) * t82;
t15 = -t17 - t140;
t119 = t15 * t122 + t149 * t79;
t69 = t72 ^ 2;
t118 = t79 * t69 * t82;
t117 = t15 * qJD(5) * t79 + t150 * t82;
t116 = -g(1) * t132 - g(2) * t133 - t59;
t112 = t125 * t70;
t109 = t79 * t72 * t122;
t107 = g(1) * t78 + g(2) * t76;
t27 = t148 * t80 + t83 * t42;
t10 = t27 * qJD(4) + t80 * t39 + t83 * t40;
t105 = -t10 * t72 + t103 * t70;
t104 = t11 * t79 - t12 * t82;
t100 = t107 * t60;
t99 = t110 - t116;
t98 = pkin(7) * t85 - t136 - t141;
t97 = -qJD(5) * t16 + t113;
t96 = t27 * t85 - t105;
t95 = -pkin(7) * qJDD(5) + (t17 - t140) * qJD(5);
t94 = -qJDD(5) * t27 + (-t103 * t72 - t9) * qJD(5);
t64 = sin(t71);
t65 = cos(t71);
t93 = -g(3) * t65 + t107 * t64;
t90 = -qJDD(5) * t62 + (t63 * t72 - t147) * qJD(5);
t89 = t116 + t91;
t87 = -qJD(5) * qJD(2) + t107 * t61 - t15 * t72 - t5 + t59;
t86 = -g(1) * (-pkin(4) * t134 + pkin(7) * t132) - g(2) * (-pkin(4) * t135 + pkin(7) * t133) - g(3) * (t61 * pkin(4) + t60 * pkin(7));
t53 = qJDD(5) * t82 - t85 * t79;
t52 = qJDD(5) * t79 + t85 * t82;
t45 = qJDD(2) + t113;
t35 = t74 * t70 - 0.2e1 * t109;
t34 = t73 * t70 + 0.2e1 * t109;
t28 = -0.2e1 * t126 * t72 * qJD(5) + 0.2e1 * t79 * t129;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t75 ^ 2 + t77 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t40 * qJD(3) + qJDD(3) * t148, -t39 * qJD(3) - t42 * qJDD(3), 0, t148 * t88 + t25 * t42 - t37 * t40 + t38 * t39 - g(3), 0, 0, 0, 0, 0, 0, t105, -t27 * t70 - t138, 0, -t17 * t10 + t103 * t8 - t110 * t27 + t18 * t9 - g(3), 0, 0, 0, 0, 0, 0, t94 * t79 - t96 * t82, t96 * t79 + t94 * t82, t27 * t112 + t125 * t138, t15 * t10 - t103 * t6 - t104 * t9 + t91 * t27 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, t53, -t52, 0, -t104 * qJD(5) + t2 * t79 + t3 * t82 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t93 + t106, g(3) * t64 + t107 * t65 - t146, 0, 0, 0, 0, 0, 0, 0, t70, t70 * t139 + t21 * t72 + (-t128 + (-pkin(3) * t72 - t33) * t80) * qJD(4) - t111 + t145, t22 * t72 + (-t72 * t123 - t70 * t80) * pkin(3) + t99, 0, t17 * t21 - t18 * t22 + (-t110 * t80 + t8 * t83 + (-t17 * t80 + t18 * t83) * qJD(4) + t93) * pkin(3), t34, t28, t52, t35, t53, 0, t90 * t79 + (-t149 - t144) * t82 + t117, t90 * t82 + (-t100 + t144) * t79 + t119, t147 * t72 * t125 + t62 * t112 + t89, t6 * t63 - t15 * t21 + t104 * t22 + ((-t104 * t83 + t15 * t80) * qJD(4) + t93) * pkin(3) + t91 * t62 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t136 + t8 + t145, t99 + t137, 0, 0, t34, t28, t52, t35, t53, 0, t95 * t79 + (-t149 - t98) * t82 + t117, t95 * t82 + (-t100 + t98) * t79 + t119, pkin(7) * t112 - t125 * t137 + t89, -t6 * pkin(4) + t91 * pkin(7) + t104 * t17 - t15 * t18 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t126 * t69, t79 * t70, t118, t129, qJDD(5), t87 * t79 + t97 * t82 + t120 + t66, t121 + (-qJDD(2) - t97) * t79 + t87 * t82, 0, 0;];
tau_reg = t1;
