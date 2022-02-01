% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:23
% DurationCPUTime: 0.58s
% Computational Cost: add. (1140->132), mult. (1902->178), div. (0->0), fcn. (1145->14), ass. (0->99)
t75 = sin(qJ(2));
t130 = pkin(1) * t75;
t105 = qJD(1) * t130;
t79 = cos(qJ(2));
t110 = qJD(1) * t79;
t67 = qJD(1) + qJD(2);
t43 = pkin(1) * t110 + t67 * pkin(2);
t71 = sin(pkin(9));
t72 = cos(pkin(9));
t23 = -t71 * t105 + t72 * t43;
t21 = t67 * pkin(3) + t23;
t24 = t72 * t105 + t71 * t43;
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t11 = t74 * t21 + t78 * t24;
t107 = qJDD(1) * t75;
t133 = pkin(1) * (qJD(2) * t110 + t107);
t123 = t79 * pkin(1);
t60 = qJDD(1) * t123;
t66 = qJDD(1) + qJDD(2);
t31 = t66 * pkin(2) - qJD(2) * t105 + t60;
t18 = -t71 * t133 + t72 * t31;
t13 = t66 * pkin(3) + t18;
t19 = t72 * t133 + t71 * t31;
t70 = qJ(1) + qJ(2);
t58 = pkin(9) + qJ(4) + t70;
t52 = cos(t58);
t136 = g(2) * t52 + t11 * qJD(4) - t78 * t13 + t74 * t19;
t62 = qJDD(4) + t66;
t125 = t62 * pkin(4);
t135 = -t125 + t136;
t129 = pkin(2) * t71;
t54 = t72 * pkin(2) + pkin(3);
t113 = t78 * t129 + t74 * t54;
t119 = t72 * t75;
t92 = pkin(1) * (-t71 * t79 - t119);
t35 = qJD(1) * t92;
t120 = t71 * t75;
t91 = pkin(1) * (t72 * t79 - t120);
t37 = qJD(1) * t91;
t63 = qJD(4) + t67;
t103 = (-t113 * qJD(4) - t78 * t35 + t74 * t37) * t63;
t93 = -t74 * t129 + t78 * t54;
t33 = -pkin(4) - t93;
t34 = pkin(8) + t113;
t81 = qJD(5) ^ 2;
t134 = -t33 * t62 - t34 * t81 + t103;
t59 = pkin(2) + t123;
t97 = -pkin(1) * t120 + t72 * t59;
t32 = pkin(3) + t97;
t39 = pkin(1) * t119 + t71 * t59;
t116 = t74 * t32 + t78 * t39;
t64 = sin(t70);
t65 = cos(t70);
t132 = g(1) * t64 - g(2) * t65;
t118 = t74 * t24;
t51 = sin(t58);
t83 = g(1) * t52 + g(2) * t51 - (qJD(4) * t21 + t19) * t78 + qJD(4) * t118 - t74 * t13;
t49 = g(1) * t51;
t36 = qJD(2) * t92;
t38 = qJD(2) * t91;
t126 = (t116 * qJD(4) - t78 * t36 + t74 * t38) * t63;
t124 = t63 * pkin(4);
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t10 = t78 * t21 - t118;
t8 = -t10 - t124;
t122 = t8 * qJD(5) * t73 + t77 * t49;
t121 = t11 * t63;
t117 = t77 * t62;
t115 = -t93 * qJD(4) + t74 * t35 + t78 * t37;
t112 = g(1) * t65 + g(2) * t64;
t68 = t73 ^ 2;
t111 = -t77 ^ 2 + t68;
t109 = t77 * qJD(5);
t108 = qJDD(3) - g(3);
t106 = t8 * t109 + t135 * t73;
t99 = qJD(1) * (-qJD(2) + t67);
t98 = qJD(2) * (-qJD(1) - t67);
t96 = t60 + t132;
t95 = t78 * t32 - t74 * t39;
t89 = pkin(8) * t81 - t121 - t125;
t14 = -pkin(4) - t95;
t15 = pkin(8) + t116;
t88 = t14 * t62 + t15 * t81 + t126;
t87 = -t62 * pkin(8) - t8 * t63 + t83;
t86 = -pkin(8) * qJDD(5) + (t10 - t124) * qJD(5);
t4 = t95 * qJD(4) + t74 * t36 + t78 * t38;
t85 = -qJDD(5) * t15 + (t14 * t63 - t4) * qJD(5);
t84 = -qJDD(5) * t34 + (t33 * t63 + t115) * qJD(5);
t82 = t49 - t136;
t80 = cos(qJ(1));
t76 = sin(qJ(1));
t61 = t63 ^ 2;
t45 = qJDD(5) * t77 - t81 * t73;
t44 = qJDD(5) * t73 + t81 * t77;
t26 = 0.2e1 * t73 * t63 * t109 + t68 * t62;
t20 = -0.2e1 * t111 * t63 * qJD(5) + 0.2e1 * t73 * t117;
t1 = [qJDD(1), g(1) * t76 - g(2) * t80, g(1) * t80 + g(2) * t76, t66, (t66 * t79 + t75 * t98) * pkin(1) + t96, ((-qJDD(1) - t66) * t75 + t79 * t98) * pkin(1) + t112, t19 * t39 + t24 * t38 + t18 * t97 + t23 * t36 - g(1) * (-t76 * pkin(1) - pkin(2) * t64) - g(2) * (t80 * pkin(1) + pkin(2) * t65), t62, t95 * t62 - t126 + t82, -t116 * t62 - t4 * t63 + t83, t26, t20, t44, t45, 0, t85 * t73 + (-t135 - t88) * t77 + t122, t85 * t77 + (t88 - t49) * t73 + t106; 0, 0, 0, t66, t99 * t130 + t96, (t79 * t99 - t107) * pkin(1) + t112, -t23 * t35 - t24 * t37 + (t18 * t72 + t19 * t71 + t132) * pkin(2), t62, t93 * t62 + t103 + t82, -t113 * t62 + t115 * t63 + t83, t26, t20, t44, t45, 0, t84 * t73 + (-t135 + t134) * t77 + t122, t84 * t77 + (-t49 - t134) * t73 + t106; 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t44; 0, 0, 0, 0, 0, 0, 0, t62, t82 + t121, t10 * t63 + t83, t26, t20, t44, t45, 0, t86 * t73 + (-t135 - t89) * t77 + t122, t86 * t77 + (t89 - t49) * t73 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t61 * t77, t111 * t61, t73 * t62, t117, qJDD(5), t108 * t77 + t87 * t73, -t108 * t73 + t87 * t77;];
tau_reg = t1;
