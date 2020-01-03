% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRP3
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:18
% EndTime: 2019-12-31 17:14:19
% DurationCPUTime: 0.89s
% Computational Cost: add. (882->191), mult. (1387->216), div. (0->0), fcn. (675->8), ass. (0->123)
t74 = sin(qJ(3));
t71 = t74 ^ 2;
t77 = cos(qJ(3));
t72 = t77 ^ 2;
t132 = t71 + t72;
t75 = sin(qJ(2));
t118 = qJDD(1) * t75;
t78 = cos(qJ(2));
t126 = qJD(2) * t78;
t69 = qJDD(1) + qJDD(2);
t20 = t69 * pkin(6) + (qJD(1) * t126 + t118) * pkin(1);
t165 = t132 * t20;
t96 = t77 * pkin(3) + t74 * qJ(4);
t141 = t72 * t69;
t142 = t71 * t69;
t164 = t141 + t142;
t108 = t132 * t78;
t131 = pkin(1) * qJD(1);
t113 = t75 * t131;
t70 = qJD(1) + qJD(2);
t36 = t70 * pkin(6) + t113;
t128 = qJD(1) * t78;
t112 = pkin(1) * t128;
t152 = t70 * pkin(2);
t37 = -t112 - t152;
t163 = t36 * t108 + t37 * t75;
t116 = qJDD(3) * qJ(4);
t140 = t74 * t36;
t16 = t77 * t20;
t4 = t116 + t16 + (qJD(4) - t140) * qJD(3);
t14 = t74 * t20;
t122 = qJDD(3) * pkin(3);
t123 = qJD(3) * t77;
t159 = t36 * t123 - t122;
t5 = qJDD(4) + t14 + t159;
t162 = t4 * t77 + t5 * t74;
t73 = qJ(1) + qJ(2);
t65 = sin(t73);
t66 = cos(t73);
t160 = g(1) * t66 + g(2) * t65;
t119 = qJD(3) * qJ(4);
t139 = t77 * t36;
t23 = t119 + t139;
t121 = t23 * qJD(3);
t102 = -qJD(3) * pkin(3) + qJD(4);
t18 = t102 + t140;
t158 = -t74 * t121 + t18 * t123 + t162;
t125 = qJD(3) * t70;
t133 = -t71 + t72;
t138 = t77 * t69;
t157 = 0.2e1 * t133 * t125 + 0.2e1 * t74 * t138;
t80 = qJD(3) ^ 2;
t156 = pkin(6) * t80;
t58 = g(1) * t65;
t76 = sin(qJ(1));
t155 = g(1) * t76;
t154 = g(2) * t66;
t153 = t69 * pkin(2);
t150 = t78 * pkin(1);
t38 = -pkin(2) - t96;
t148 = t38 * t69;
t147 = t38 * t70;
t60 = t75 * pkin(1) + pkin(6);
t146 = t60 * t80;
t145 = t65 * t74;
t144 = t66 * t74;
t143 = t70 * t74;
t137 = g(1) * t145 - g(2) * t144;
t136 = t66 * pkin(2) + t65 * pkin(6);
t134 = -qJD(2) * t113 + qJDD(1) * t150;
t129 = pkin(6) * qJDD(3);
t127 = qJD(2) * t75;
t124 = qJD(3) * t74;
t120 = t74 * qJD(4);
t117 = qJDD(3) * t60;
t103 = t70 * t113;
t48 = t77 * t58;
t115 = t77 * t103 + t112 * t124 + t48;
t114 = pkin(1) * t126;
t111 = t70 * t127;
t19 = -t134 - t153;
t109 = -t19 - t154;
t106 = t37 * t123 + t19 * t74 - t137;
t105 = t96 * t66 + t136;
t101 = t123 * t143;
t55 = t66 * pkin(6);
t100 = g(1) * (-t65 * pkin(2) + t55);
t99 = -t153 + t156;
t79 = cos(qJ(1));
t97 = -g(2) * t79 + t155;
t95 = pkin(3) * t74 - qJ(4) * t77;
t94 = t18 * t74 + t23 * t77;
t93 = t134 + t58 - t154;
t92 = g(1) * t144 + g(2) * t145 - g(3) * t77 - t14;
t91 = t132 * t70 * t114 + t164 * t60 - t160;
t1 = (t95 * qJD(3) - t120) * t70 + t148 - t134;
t90 = -t1 - t148 - t156;
t29 = t38 - t150;
t89 = t29 * t70 - t114;
t88 = -qJDD(4) + t92;
t24 = pkin(3) * t124 - t77 * t119 - t120;
t12 = pkin(1) * t127 + t24;
t87 = -t12 * t70 - t29 * t69 - t1 - t146;
t61 = -pkin(2) - t150;
t86 = pkin(1) * t111 + t61 * t69 + t146;
t85 = -t70 * t108 * t131 + t164 * pkin(6) - t160;
t84 = -g(1) * t55 - t38 * t58;
t83 = -t117 + (t61 * t70 - t114) * qJD(3);
t82 = (t18 * t77 - t23 * t74) * qJD(3) + t162;
t68 = t70 ^ 2;
t67 = t79 * pkin(1);
t52 = t74 * t69;
t42 = t74 * t68 * t77;
t40 = qJDD(3) * t77 - t80 * t74;
t39 = qJDD(3) * t74 + t80 * t77;
t30 = t133 * t68;
t27 = t37 * t124;
t25 = t95 * t70;
t22 = -0.2e1 * t101 + t141;
t21 = 0.2e1 * t101 + t142;
t9 = -t112 + t147;
t6 = t9 * t124;
t2 = [0, 0, 0, 0, 0, qJDD(1), t97, g(1) * t79 + g(2) * t76, 0, 0, 0, 0, 0, 0, 0, t69, (t69 * t78 - t111) * pkin(1) + t93, ((-qJDD(1) - t69) * t75 + (-qJD(1) - t70) * t126) * pkin(1) + t160, 0, (t97 + (t75 ^ 2 + t78 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t21, t157, t39, t22, t40, 0, t27 + t48 + t83 * t74 + (t109 - t86) * t77, t74 * t86 + t77 * t83 + t106, t91 + t165, t19 * t61 - t100 - g(2) * (t67 + t136) + t60 * t165 + (t163 * qJD(2) + t155) * pkin(1), t21, t39, -t157, 0, -t40, t22, t48 + t6 + (t89 * qJD(3) - t117) * t74 + (t87 - t154) * t77, t91 + t158, (t117 + (-t89 - t9) * qJD(3)) * t77 + t87 * t74 + t137, t1 * t29 + t9 * t12 - g(2) * (t67 + t105) + (t94 * t126 + t155) * pkin(1) + t82 * t60 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t93 + t103, (-t118 + (-qJD(2) + t70) * t128) * pkin(1) + t160, 0, 0, t21, t157, t39, t22, t40, 0, t27 + (-pkin(2) * t125 - t129) * t74 + (t109 - t99) * t77 + t115, (-t129 + (t112 - t152) * qJD(3)) * t77 + (t99 - t103) * t74 + t106, t85 + t165, -t19 * pkin(2) + pkin(6) * t165 - g(2) * t136 - t163 * t131 - t100, t21, t39, -t157, 0, -t40, t22, t6 + (t38 * t125 - t129) * t74 + (-t24 * t70 - t154 + t90) * t77 + t115, t85 + t158, (t129 + (-t112 - t9 - t147) * qJD(3)) * t77 + ((-t24 + t113) * t70 + t90) * t74 + t137, t1 * t38 + t9 * t24 - g(2) * t105 + (-t75 * t9 - t78 * t94) * t131 + t82 * pkin(6) + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t30, t52, t42, t138, qJDD(3), -t37 * t143 + t92, g(3) * t74 - t16 + (-t37 * t70 + t160) * t77, 0, 0, -t42, t52, t30, qJDD(3), -t138, t42, 0.2e1 * t122 + (t25 * t77 - t74 * t9) * t70 + t88, -t95 * t69 + ((t23 - t119) * t74 + (t102 - t18) * t77) * t70, 0.2e1 * t116 + 0.2e1 * qJD(3) * qJD(4) + t16 + (t25 * t70 - g(3)) * t74 + (t70 * t9 - t160) * t77, t4 * qJ(4) - t5 * pkin(3) - t9 * t25 - t18 * t139 - g(3) * t96 + (qJD(4) + t140) * t23 + t160 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t42, t52, -t71 * t68 - t80, t9 * t143 - t121 + t159 - t88;];
tau_reg = t2;
