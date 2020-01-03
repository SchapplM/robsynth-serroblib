% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP6
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:12
% EndTime: 2019-12-31 17:19:15
% DurationCPUTime: 1.18s
% Computational Cost: add. (952->224), mult. (2227->319), div. (0->0), fcn. (1380->6), ass. (0->119)
t73 = sin(qJ(2));
t115 = qJD(3) * t73;
t156 = qJD(1) * t115 - qJDD(2);
t111 = t73 * qJDD(1);
t76 = cos(qJ(2));
t122 = qJD(1) * t76;
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t11 = ((qJD(3) + t122) * qJD(2) + t111) * t72 + t156 * t75;
t110 = qJD(1) * qJD(2);
t67 = t76 * qJDD(1);
t155 = -t73 * t110 + t67;
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t96 = g(1) * t77 + g(2) * t74;
t136 = t72 * t76;
t28 = t74 * t136 + t75 * t77;
t133 = t76 * t77;
t30 = -t72 * t133 + t74 * t75;
t154 = -g(1) * t30 + g(2) * t28;
t102 = t76 * t110;
t62 = pkin(5) * t111;
t27 = -qJDD(2) * pkin(2) + pkin(5) * t102 + t62;
t58 = -qJD(3) + t122;
t142 = g(3) * t76;
t81 = -t96 * t73 + t142;
t153 = -qJD(3) * pkin(6) * t58 + t27 + t81;
t119 = qJD(2) * t72;
t123 = qJD(1) * t73;
t44 = t75 * t123 + t119;
t151 = t44 ^ 2;
t48 = -pkin(2) * t76 - pkin(6) * t73 - pkin(1);
t36 = t48 * qJD(1);
t64 = pkin(5) * t122;
t52 = qJD(2) * pkin(6) + t64;
t16 = t75 * t36 - t52 * t72;
t7 = -qJ(4) * t44 + t16;
t6 = -pkin(3) * t58 + t7;
t150 = -t7 + t6;
t149 = pkin(3) * t72;
t148 = pkin(5) * t72;
t143 = g(3) * t73;
t132 = qJ(4) + pkin(6);
t100 = qJD(3) * t132;
t108 = t72 * t123;
t97 = pkin(2) * t73 - pkin(6) * t76;
t46 = t97 * qJD(1);
t128 = pkin(5) * t108 + t75 * t46;
t124 = qJ(4) * t76;
t88 = pkin(3) * t73 - t75 * t124;
t141 = -t88 * qJD(1) - t72 * qJD(4) - t75 * t100 - t128;
t113 = t75 * qJD(2);
t10 = -qJD(3) * t113 + (-t102 - t111) * t75 + t156 * t72;
t140 = t10 * t72;
t42 = t108 - t113;
t139 = t42 * t58;
t138 = t44 * t58;
t137 = t44 * t75;
t135 = t73 * t75;
t134 = t75 * t76;
t112 = t75 * qJD(4);
t32 = t72 * t46;
t131 = -t72 * t100 + t112 - t32 - (-pkin(5) * t135 - t72 * t124) * qJD(1);
t114 = qJD(3) * t75;
t47 = t97 * qJD(2);
t130 = t48 * t114 + t72 * t47;
t118 = qJD(2) * t73;
t129 = t118 * t148 + t75 * t47;
t59 = pkin(5) * t134;
t127 = t72 * t48 + t59;
t69 = t73 ^ 2;
t126 = -t76 ^ 2 + t69;
t125 = qJ(4) * t73;
t121 = qJD(2) * t42;
t120 = qJD(2) * t44;
t117 = qJD(2) * t76;
t116 = qJD(3) * t72;
t107 = t58 * t113;
t106 = t58 * t116;
t105 = t58 * t114;
t104 = pkin(5) + t149;
t51 = -qJD(2) * pkin(2) + pkin(5) * t123;
t26 = pkin(5) * t155 + qJDD(2) * pkin(6);
t99 = -qJD(3) * t36 - t26;
t95 = g(1) * t74 - g(2) * t77;
t18 = qJD(1) * t47 + t48 * qJDD(1);
t14 = t75 * t18;
t94 = t52 * t114 - t14;
t17 = t36 * t72 + t52 * t75;
t8 = -qJ(4) * t42 + t17;
t93 = t6 * t75 + t72 * t8;
t61 = pkin(3) * t75 + pkin(2);
t92 = t132 * t73 + t61 * t76;
t39 = qJDD(3) - t155;
t90 = -pkin(6) * t39 + qJD(3) * t51;
t89 = pkin(1) + t92;
t86 = t72 * t39 - t105;
t85 = t75 * t39 + t106;
t84 = t36 * t114 - t52 * t116 + t72 * t18 + t75 * t26;
t83 = -0.2e1 * pkin(1) * t110 - pkin(5) * qJDD(2);
t79 = qJD(1) ^ 2;
t82 = pkin(1) * t79 + t96;
t5 = pkin(3) * t11 + qJDD(4) + t27;
t78 = qJD(2) ^ 2;
t80 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t78 + t95;
t50 = t132 * t75;
t49 = t132 * t72;
t41 = t75 * t48;
t38 = t42 ^ 2;
t31 = t75 * t133 + t72 * t74;
t29 = -t74 * t134 + t72 * t77;
t20 = pkin(3) * t42 + qJD(4) + t51;
t19 = -t72 * t125 + t127;
t15 = -t75 * t125 + t41 + (-pkin(3) - t148) * t76;
t4 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t135 + (-qJD(4) * t73 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t76) * t72 + t130;
t3 = -t73 * t112 + t88 * qJD(2) + (-t59 + (-t48 + t125) * t72) * qJD(3) + t129;
t2 = -qJ(4) * t11 - qJD(4) * t42 + t84;
t1 = pkin(3) * t39 + qJ(4) * t10 - t17 * qJD(3) - qJD(4) * t44 - t72 * t26 + t14;
t9 = [qJDD(1), t95, t96, qJDD(1) * t69 + 0.2e1 * t73 * t102, -0.2e1 * t126 * t110 + 0.2e1 * t73 * t67, qJDD(2) * t73 + t76 * t78, qJDD(2) * t76 - t73 * t78, 0, t83 * t73 + t80 * t76, -t80 * t73 + t83 * t76, -t10 * t135 + (t113 * t76 - t115 * t72) * t44, (-t42 * t75 - t44 * t72) * t117 + (t140 - t11 * t75 + (t42 * t72 - t137) * qJD(3)) * t73, (t10 - t107) * t76 + (t85 + t120) * t73, (t119 * t58 + t11) * t76 + (-t86 - t121) * t73, -t118 * t58 - t39 * t76, -(-t116 * t48 + t129) * t58 + t41 * t39 - g(1) * t29 - g(2) * t31 + ((t105 + t121) * pkin(5) + (-pkin(5) * t39 + qJD(2) * t51 - t99) * t72 + t94) * t76 + (pkin(5) * t11 + qJD(2) * t16 + t114 * t51 + t27 * t72) * t73, t130 * t58 - t127 * t39 - g(1) * t28 - g(2) * t30 + (t51 * t113 + (-t106 + t120) * pkin(5) + t84) * t76 + (-t51 * t116 - t17 * qJD(2) + t27 * t75 + (-t10 - t107) * pkin(5)) * t73, t10 * t15 - t11 * t19 - t3 * t44 - t4 * t42 - t93 * t117 + (-t1 * t75 - t2 * t72 + (t6 * t72 - t75 * t8) * qJD(3) + t95) * t73, t1 * t15 + t2 * t19 + t6 * t3 + t8 * t4 + t20 * t104 * t117 + (pkin(3) * t114 * t20 + t104 * t5) * t73 + (-g(1) * t104 - g(2) * t89) * t77 + (g(1) * t89 - g(2) * t104) * t74; 0, 0, 0, -t73 * t79 * t76, t126 * t79, t111, t67, qJDD(2), t82 * t73 - t142 - t62, t143 + (-pkin(5) * qJDD(1) + t82) * t76, -t58 * t137 - t140, (-t10 + t139) * t75 + (-t11 + t138) * t72, (t58 * t134 - t44 * t73) * qJD(1) + t86, (-t58 * t136 + t42 * t73) * qJD(1) + t85, t58 * t123, -pkin(2) * t11 + t128 * t58 + t90 * t72 + (-t16 * t73 + (-pkin(5) * t42 - t51 * t72) * t76) * qJD(1) - t153 * t75, pkin(2) * t10 - t32 * t58 + t90 * t75 + (-t51 * t134 + t17 * t73 + (t58 * t135 - t44 * t76) * pkin(5)) * qJD(1) + t153 * t72, -t143 - t1 * t72 - t10 * t49 - t11 * t50 + t2 * t75 - t141 * t44 - t131 * t42 - t93 * qJD(3) + (qJD(1) * t93 - t96) * t76, t2 * t50 - t1 * t49 - t5 * t61 - g(3) * t92 + t131 * t8 + t141 * t6 + (-t58 * t149 - t64) * t20 + t96 * (-t132 * t76 + t61 * t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t42, -t38 + t151, -t10 - t139, -t11 - t138, t39, -t17 * t58 - t44 * t51 + (t99 + t143) * t72 - t94 + t154, g(1) * t31 - g(2) * t29 + g(3) * t135 - t16 * t58 + t42 * t51 - t84, pkin(3) * t10 - t150 * t42, t150 * t8 + (t72 * t143 - t20 * t44 + t1 + t154) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 - t151, t42 * t8 + t44 * t6 + t5 + t81;];
tau_reg = t9;
