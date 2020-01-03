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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:59:13
% EndTime: 2020-01-03 11:59:16
% DurationCPUTime: 0.78s
% Computational Cost: add. (1065->176), mult. (1719->225), div. (0->0), fcn. (978->12), ass. (0->120)
t148 = qJ(5) + pkin(7);
t82 = qJ(1) + qJ(2);
t69 = pkin(8) + t82;
t58 = sin(t69);
t59 = cos(t69);
t105 = g(2) * t59 + g(3) * t58;
t130 = pkin(1) * qJD(2);
t117 = qJD(1) * t130;
t87 = sin(qJ(2));
t124 = qJDD(1) * t87;
t90 = cos(qJ(2));
t147 = pkin(1) * t124 + t90 * t117;
t146 = pkin(1) * t90;
t67 = qJDD(1) * t146;
t78 = qJDD(1) + qJDD(2);
t28 = pkin(2) * t78 - t87 * t117 + t67;
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t13 = -t147 * t83 + t28 * t84;
t100 = pkin(3) * t78 - t105 + t13;
t104 = g(2) * t58 - g(3) * t59;
t145 = pkin(2) * t84;
t89 = cos(qJ(4));
t144 = pkin(4) * t89;
t143 = g(1) * t89;
t131 = pkin(1) * qJD(1);
t121 = t90 * t131;
t79 = qJD(1) + qJD(2);
t41 = pkin(2) * t79 + t121;
t122 = t87 * t131;
t54 = t84 * t122;
t24 = t83 * t41 + t54;
t113 = t148 * t79 + t24;
t86 = sin(qJ(4));
t11 = t89 * qJD(3) - t113 * t86;
t129 = qJD(4) * pkin(4);
t8 = t11 + t129;
t139 = -t11 + t8;
t138 = t83 * t87;
t137 = t84 * t87;
t136 = t86 * t78;
t135 = t89 * t78;
t66 = pkin(2) + t146;
t134 = pkin(1) * t137 + t83 * t66;
t80 = t86 ^ 2;
t81 = t89 ^ 2;
t133 = -t80 - t81;
t132 = t80 - t81;
t31 = pkin(7) + t134;
t128 = -qJ(5) - t31;
t60 = pkin(2) * t83 + pkin(7);
t127 = -qJ(5) - t60;
t126 = qJD(4) * t86;
t125 = qJDD(3) - g(1);
t72 = sin(t82);
t63 = pkin(2) * t72;
t65 = pkin(3) + t144;
t123 = -t148 * t59 + t58 * t65 + t63;
t14 = t147 * t84 + t83 * t28;
t120 = pkin(4) * t126;
t119 = t79 * t126;
t73 = cos(t82);
t116 = g(2) * t72 - g(3) * t73;
t53 = t83 * t122;
t23 = t84 * t41 - t53;
t115 = -pkin(1) * t138 + t66 * t84;
t19 = -pkin(3) * t79 - t23;
t114 = t19 * qJD(4) * t89 - t100 * t86;
t112 = qJD(4) * t128;
t111 = qJD(4) * t127;
t110 = qJD(1) * (-qJD(2) + t79);
t109 = qJD(2) * (-qJD(1) - t79);
t101 = t113 * qJD(4);
t10 = pkin(7) * t78 + t14;
t93 = qJ(5) * t78 + qJD(3) * qJD(4) + qJD(5) * t79 + t10;
t3 = (qJDD(3) - t101) * t86 + t93 * t89;
t108 = t3 * t89 - t104;
t30 = -pkin(3) - t115;
t64 = pkin(2) * t73;
t106 = t148 * t58 + t59 * t65 + t64;
t103 = -g(2) * t73 - g(3) * t72;
t12 = t86 * qJD(3) + t113 * t89;
t102 = t12 * t89 - t8 * t86;
t99 = t103 + t67;
t33 = (t83 * t90 + t137) * t130;
t92 = qJD(4) ^ 2;
t98 = t30 * t78 + t31 * t92 + t33 * t79;
t32 = t83 * t121 + t54;
t61 = -pkin(3) - t145;
t97 = -t32 * t79 + t60 * t92 + t61 * t78;
t96 = -t19 * t79 - t10 + t104;
t35 = (t84 * t90 - t138) * t130;
t95 = -qJDD(4) * t31 + (t30 * t79 - t35) * qJD(4);
t34 = t84 * t121 - t53;
t94 = -qJDD(4) * t60 + (t61 * t79 + t34) * qJD(4);
t4 = pkin(4) * t119 - t65 * t78 + qJDD(5) - t13;
t91 = cos(qJ(1));
t88 = sin(qJ(1));
t77 = t79 ^ 2;
t76 = t91 * pkin(1);
t75 = t88 * pkin(1);
t74 = t89 * qJ(5);
t70 = t89 * qJD(5);
t68 = t89 * qJDD(3);
t43 = qJDD(4) * t89 - t86 * t92;
t42 = qJDD(4) * t86 + t89 * t92;
t38 = t60 * t89 + t74;
t37 = t127 * t86;
t29 = 0.2e1 * t89 * t119 + t78 * t80;
t27 = -t86 * qJD(5) + t89 * t111;
t26 = t86 * t111 + t70;
t22 = t31 * t89 + t74;
t21 = t128 * t86;
t18 = -0.2e1 * t132 * t79 * qJD(4) + 0.2e1 * t86 * t135;
t16 = t19 * t126;
t15 = -t65 * t79 + qJD(5) - t23;
t6 = (-qJD(5) - t35) * t86 + t89 * t112;
t5 = t86 * t112 + t89 * t35 + t70;
t2 = qJDD(4) * pkin(4) - t89 * t101 - t93 * t86 + t68;
t1 = [qJDD(1), -g(2) * t91 - g(3) * t88, g(2) * t88 - g(3) * t91, t78, (t87 * t109 + t78 * t90) * pkin(1) + t99, ((-qJDD(1) - t78) * t87 + t90 * t109) * pkin(1) + t116, t14 * t134 + t24 * t35 + t13 * t115 - t23 * t33 - g(2) * (t64 + t76) - g(3) * (t63 + t75), t29, t18, t42, t43, 0, t16 + t95 * t86 + (t100 - t98) * t89, t98 * t86 + t95 * t89 + t114, (t22 * t78 + t5 * t79 + (-t21 * t79 - t8) * qJD(4)) * t89 + (-t21 * t78 - t6 * t79 - t2 + (-t22 * t79 - t12) * qJD(4)) * t86 + t108, t3 * t22 + t12 * t5 + t2 * t21 + t8 * t6 + t4 * (t30 - t144) + t15 * (t33 + t120) - g(2) * (t106 + t76) - g(3) * (t75 + t123); 0, 0, 0, t78, t87 * pkin(1) * t110 + t99, (t90 * t110 - t124) * pkin(1) + t116, t23 * t32 - t24 * t34 + (t13 * t84 + t14 * t83 + t103) * pkin(2), t29, t18, t42, t43, 0, t16 + t94 * t86 + (t100 - t97) * t89, t97 * t86 + t94 * t89 + t114, (-qJD(4) * t8 + t38 * t78) * t89 + (-qJD(4) * t12 - t37 * t78 - t2) * t86 + (t26 * t89 - t27 * t86 + t133 * t34 + (-t37 * t89 - t38 * t86) * qJD(4)) * t79 + t108, t3 * t38 + t2 * t37 + t4 * (-t65 - t145) - g(2) * t106 - g(3) * t123 + (t34 * t86 + t27) * t8 + (-t32 + t120) * t15 + (-t34 * t89 + t26) * t12; 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, t43, -t42, 0, t102 * qJD(4) + t2 * t89 + t3 * t86 - g(1); 0, 0, 0, 0, 0, 0, 0, -t86 * t77 * t89, t132 * t77, t136, t135, qJDD(4), t96 * t86 - t143 + t68, -t125 * t86 + t96 * t89, -pkin(4) * t136 + (-t129 + t139) * t89 * t79, t139 * t12 + (-t143 + t2 + (-t15 * t79 + t104) * t86) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t77, -t102 * t79 + t105 + t4;];
tau_reg = t1;
