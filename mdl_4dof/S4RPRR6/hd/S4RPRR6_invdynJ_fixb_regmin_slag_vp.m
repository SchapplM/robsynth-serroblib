% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:43
% EndTime: 2019-12-31 16:52:45
% DurationCPUTime: 0.74s
% Computational Cost: add. (903->154), mult. (2190->211), div. (0->0), fcn. (1732->12), ass. (0->99)
t80 = cos(pkin(7));
t85 = cos(qJ(3));
t123 = t85 * t80;
t79 = sin(pkin(7));
t82 = sin(qJ(3));
t124 = t79 * t82;
t48 = -t123 + t124;
t40 = t48 * qJD(1);
t49 = t79 * t85 + t80 * t82;
t41 = t49 * qJD(1);
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t18 = t84 * t40 + t41 * t81;
t78 = qJD(3) + qJD(4);
t126 = t18 * t78;
t117 = qJD(4) * t84;
t118 = qJD(4) * t81;
t111 = qJD(1) * t124;
t114 = t80 * qJDD(1);
t115 = t79 * qJDD(1);
t119 = qJD(3) * t85;
t112 = t80 * qJD(1) * t119 + t82 * t114 + t85 * t115;
t25 = -qJD(3) * t111 + t112;
t103 = -t85 * t114 + t82 * t115;
t43 = t49 * qJD(3);
t26 = qJD(1) * t43 + t103;
t4 = -t40 * t117 - t41 * t118 + t84 * t25 - t81 * t26;
t141 = t4 + t126;
t98 = -t40 * t81 + t84 * t41;
t140 = t98 * t18;
t127 = t98 * t78;
t5 = t98 * qJD(4) + t81 * t25 + t84 * t26;
t139 = -t5 + t127;
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t105 = g(1) * t86 + g(2) * t83;
t138 = -t18 ^ 2 + t98 ^ 2;
t122 = pkin(5) + qJ(2);
t58 = t122 * t79;
t50 = qJD(1) * t58;
t59 = t122 * t80;
t51 = qJD(1) * t59;
t97 = t50 * t82 - t51 * t85;
t13 = -pkin(6) * t40 - t97;
t113 = qJD(1) * qJD(2);
t132 = t122 * qJDD(1) + t113;
t33 = t132 * t79;
t34 = t132 * t80;
t108 = -t85 * t33 - t82 * t34;
t2 = qJDD(3) * pkin(3) - pkin(6) * t25 + t97 * qJD(3) + t108;
t69 = -pkin(2) * t80 - pkin(1);
t53 = t69 * qJD(1) + qJD(2);
t29 = t40 * pkin(3) + t53;
t77 = pkin(7) + qJ(3);
t73 = qJ(4) + t77;
t67 = sin(t73);
t68 = cos(t73);
t137 = t29 * t18 + t13 * t118 + g(3) * t67 + (-t13 * t78 - t2) * t81 + t105 * t68;
t121 = -t82 * t58 + t85 * t59;
t135 = -t85 * t50 - t51 * t82;
t134 = qJ(2) * qJDD(1);
t12 = -pkin(6) * t41 + t135;
t11 = qJD(3) * pkin(3) + t12;
t128 = t13 * t84;
t102 = -t11 * t81 - t128;
t99 = -t82 * t33 + t85 * t34;
t3 = -pkin(6) * t26 + t135 * qJD(3) + t99;
t133 = -g(3) * t68 + t102 * qJD(4) + t105 * t67 + t84 * t2 - t29 * t98 - t81 * t3;
t131 = pkin(3) * t43;
t120 = t79 ^ 2 + t80 ^ 2;
t116 = qJDD(1) * pkin(1);
t110 = t120 * qJD(1) ^ 2;
t107 = -t85 * t58 - t59 * t82;
t106 = 0.2e1 * t120;
t104 = g(1) * t83 - g(2) * t86;
t15 = -pkin(6) * t49 + t107;
t16 = -pkin(6) * t48 + t121;
t101 = t15 * t84 - t16 * t81;
t100 = t15 * t81 + t16 * t84;
t27 = t84 * t48 + t49 * t81;
t28 = -t48 * t81 + t49 * t84;
t52 = t69 * qJDD(1) + qJDD(2);
t96 = -t104 - t116;
t70 = qJDD(2) - t116;
t95 = -t70 - t96;
t93 = -t58 * t119 + qJD(2) * t123 + (-qJD(2) * t79 - qJD(3) * t59) * t82;
t91 = t106 * t113 - t105;
t89 = -t49 * qJD(2) - t121 * qJD(3);
t74 = qJDD(3) + qJDD(4);
t72 = cos(t77);
t71 = sin(t77);
t42 = t48 * qJD(3);
t31 = pkin(3) * t48 + t69;
t14 = t26 * pkin(3) + t52;
t9 = pkin(6) * t42 + t89;
t8 = -pkin(6) * t43 + t93;
t7 = t28 * qJD(4) - t81 * t42 + t84 * t43;
t6 = -t27 * qJD(4) - t84 * t42 - t81 * t43;
t1 = [qJDD(1), t104, t105, t95 * t80, -t95 * t79, t106 * t134 + t91, (t104 - t70) * pkin(1) + (t120 * t134 + t91) * qJ(2), t25 * t49 - t41 * t42, -t25 * t48 - t26 * t49 + t40 * t42 - t41 * t43, -qJD(3) * t42 + qJDD(3) * t49, -qJD(3) * t43 - qJDD(3) * t48, 0, t89 * qJD(3) + t107 * qJDD(3) + t104 * t72 + t69 * t26 + t53 * t43 + t52 * t48, -t93 * qJD(3) - t121 * qJDD(3) - t104 * t71 + t69 * t25 - t53 * t42 + t52 * t49, t28 * t4 + t6 * t98, -t18 * t6 - t27 * t4 - t28 * t5 - t7 * t98, t28 * t74 + t6 * t78, -t27 * t74 - t7 * t78, 0, t18 * t131 + t31 * t5 + t14 * t27 + t29 * t7 + (-t100 * qJD(4) - t81 * t8 + t84 * t9) * t78 + t101 * t74 + t104 * t68, t98 * t131 + t31 * t4 + t14 * t28 + t29 * t6 - (t101 * qJD(4) + t84 * t8 + t81 * t9) * t78 - t100 * t74 - t104 * t67; 0, 0, 0, -t114, t115, -t110, -qJ(2) * t110 + qJDD(2) + t96, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t41 + t103, (-t40 - t111) * qJD(3) + t112, 0, 0, 0, 0, 0, t5 + t127, t4 - t126; 0, 0, 0, 0, 0, 0, 0, t41 * t40, -t40 ^ 2 + t41 ^ 2, (t40 - t111) * qJD(3) + t112, -t103, qJDD(3), -g(3) * t72 + t105 * t71 - t53 * t41 + t108, g(3) * t71 + t105 * t72 + t53 * t40 - t99, t140, t138, t141, t139, t74, -(-t12 * t81 - t128) * t78 + (-t78 * t118 - t18 * t41 + t74 * t84) * pkin(3) + t133, (-qJD(4) * t11 + t12 * t78 - t3) * t84 + (-t78 * t117 - t41 * t98 - t74 * t81) * pkin(3) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t138, t141, t139, t74, -t102 * t78 + t133, (-t3 + (-qJD(4) + t78) * t11) * t84 + t137;];
tau_reg = t1;
