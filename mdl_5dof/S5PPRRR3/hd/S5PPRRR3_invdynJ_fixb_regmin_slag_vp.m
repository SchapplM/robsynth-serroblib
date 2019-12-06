% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRR3
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
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:02
% EndTime: 2019-12-05 15:17:07
% DurationCPUTime: 1.34s
% Computational Cost: add. (702->187), mult. (1618->281), div. (0->0), fcn. (1325->12), ass. (0->120)
t69 = sin(pkin(9));
t134 = qJD(1) * t69;
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t154 = t78 * qJD(2) - t75 * t134;
t77 = cos(qJ(4));
t61 = -t77 * pkin(4) - pkin(3);
t24 = t61 * qJD(3) - t154;
t158 = t24 + t154;
t65 = qJD(4) + qJD(5);
t151 = pkin(6) + pkin(7);
t39 = t75 * qJD(2) + t78 * t134;
t157 = t39 * qJD(3);
t73 = sin(qJ(5));
t74 = sin(qJ(4));
t76 = cos(qJ(5));
t41 = t73 * t77 + t76 * t74;
t155 = t41 * t65;
t156 = qJD(3) * t155;
t115 = t77 * qJDD(3);
t117 = t74 * qJDD(3);
t7 = -t76 * t115 + t73 * t117 + t156;
t40 = t73 * t74 - t76 * t77;
t90 = t40 * t65;
t153 = t65 * t77;
t71 = cos(pkin(9));
t118 = t71 * qJDD(1);
t70 = sin(pkin(8));
t72 = cos(pkin(8));
t152 = (g(1) * t72 + g(2) * t70) * t69 + t118;
t131 = qJD(3) * t77;
t110 = t76 * t131;
t132 = qJD(3) * t74;
t111 = t73 * t132;
t35 = -t110 + t111;
t37 = -t73 * t131 - t76 * t132;
t150 = t37 * t35;
t64 = qJDD(4) + qJDD(5);
t149 = t40 * t64;
t148 = t41 * t64;
t147 = t69 * t70;
t146 = t69 * t72;
t145 = t69 * t75;
t144 = t69 * t78;
t143 = t70 * t75;
t142 = t70 * t78;
t141 = t72 * t75;
t140 = t72 * t78;
t105 = t151 * qJD(3) + t39;
t133 = qJD(1) * t71;
t51 = t74 * t133;
t13 = t105 * t77 - t51;
t139 = t76 * t13;
t80 = qJD(3) ^ 2;
t138 = t80 * t78;
t66 = t74 ^ 2;
t137 = -t77 ^ 2 + t66;
t79 = qJD(4) ^ 2;
t136 = t79 + t80;
t135 = qJD(3) * pkin(3);
t130 = qJD(3) * t78;
t129 = qJD(4) * t74;
t127 = qJD(5) * t73;
t126 = qJD(5) * t76;
t124 = qJDD(3) * pkin(3);
t121 = qJDD(1) - g(3);
t120 = qJDD(1) * t69;
t119 = qJDD(4) * t74;
t116 = t75 * qJDD(2);
t114 = t78 * qJDD(3);
t113 = qJD(3) * qJD(4);
t112 = qJD(3) * t145;
t108 = qJD(4) * t151;
t107 = t74 * t113;
t106 = t77 * t113;
t19 = qJDD(3) * pkin(6) + qJD(3) * t154 + t78 * t120 + t116;
t104 = pkin(7) * qJDD(3) + t19;
t103 = pkin(4) * t129 - t39;
t28 = t71 * t142 - t141;
t30 = t71 * t140 + t143;
t102 = g(1) * t30 + g(2) * t28;
t12 = -t105 * t74 - t77 * t133;
t11 = qJD(4) * pkin(4) + t12;
t99 = -t73 * t11 - t139;
t33 = -t74 * t144 - t71 * t77;
t34 = t77 * t144 - t71 * t74;
t98 = t76 * t33 - t73 * t34;
t97 = t73 * t33 + t76 * t34;
t44 = t151 * t74;
t45 = t151 * t77;
t96 = -t76 * t44 - t73 * t45;
t95 = -t73 * t44 + t76 * t45;
t94 = -t80 * t75 + t114;
t93 = qJDD(3) * t75 + t138;
t92 = -t78 * qJDD(2) + t75 * t120 + t157;
t88 = t107 - t115;
t87 = -g(1) * (-t71 * t141 + t142) - g(2) * (-t71 * t143 - t140) + g(3) * t145;
t31 = -t154 - t135;
t86 = -t31 * qJD(3) + t102 - t19;
t18 = t92 - t124;
t85 = -pkin(6) * qJDD(4) + (t31 + t154 - t135) * qJD(4);
t6 = qJD(5) * t110 - t65 * t111 + t73 * t115 + (t106 + t117) * t76;
t84 = t87 + t157;
t2 = qJDD(4) * pkin(4) + t51 * qJD(4) - t104 * t74 + (-t105 * qJD(4) - t118) * t77;
t68 = qJ(4) + qJ(5);
t62 = sin(t68);
t63 = cos(t68);
t83 = -g(1) * (-t62 * t146 - t30 * t63) - g(2) * (-t62 * t147 - t28 * t63) - g(3) * (-t63 * t144 + t71 * t62) + t24 * t35 + t13 * t127 + (-t13 * t65 - t2) * t73;
t82 = -pkin(6) * t79 + t124 - t18 + t84;
t3 = t12 * qJD(4) + t104 * t77 - t74 * t118;
t81 = -g(1) * (t63 * t146 - t30 * t62) - g(2) * (t63 * t147 - t28 * t62) - g(3) * (-t62 * t144 - t71 * t63) + t99 * qJD(5) + t76 * t2 + t24 * t37 - t73 * t3;
t43 = t77 * t108;
t42 = t74 * t108;
t23 = t33 * qJD(4) - t77 * t112;
t22 = -t34 * qJD(4) + t74 * t112;
t9 = t88 * pkin(4) + t18;
t8 = -t35 ^ 2 + t37 ^ 2;
t5 = -t37 * t65 - t7;
t4 = t35 * t65 + t6;
t1 = [t121, -g(3) + (t69 ^ 2 + t71 ^ 2) * qJDD(1), 0, -t93 * t69, -t94 * t69, 0, 0, 0, 0, 0, t22 * qJD(4) + t33 * qJDD(4) + (-t77 * t138 + t88 * t75) * t69, -t23 * qJD(4) - t34 * qJDD(4) + (t75 * t106 + t93 * t74) * t69, 0, 0, 0, 0, 0, (-t97 * qJD(5) + t76 * t22 - t73 * t23) * t65 + t98 * t64 + (t35 * t130 + t75 * t7) * t69, -(t98 * qJD(5) + t73 * t22 + t76 * t23) * t65 - t97 * t64 + (-t37 * t130 + t75 * t6) * t69; 0, -g(1) * t70 + g(2) * t72 + qJDD(2), 0, t94, -t93, 0, 0, 0, 0, 0, (-0.2e1 * t107 + t115) * t78 + (-t136 * t77 - t119) * t75, (-qJDD(4) * t75 - 0.2e1 * t78 * t113) * t77 + (t136 * t75 - t114) * t74, 0, 0, 0, 0, 0, (-t7 - t156) * t78 + ((t74 * t127 + t73 * t129 - t153 * t76) * t65 - t148 + qJD(3) * t35) * t75, (qJD(3) * t90 - t6) * t78 + (-(-t74 * t126 - t76 * t129 - t153 * t73) * t65 + t149 - qJD(3) * t37) * t75; 0, 0, qJDD(3), t84 - t92, -t121 * t144 + t102 - t116, t66 * qJDD(3) + 0.2e1 * t74 * t106, -0.2e1 * t137 * t113 + 0.2e1 * t74 * t115, t79 * t77 + t119, qJDD(4) * t77 - t79 * t74, 0, t85 * t74 + t82 * t77, -t82 * t74 + t85 * t77, t37 * t90 + t6 * t41, t155 * t37 + t35 * t90 - t6 * t40 - t41 * t7, -t65 * t90 + t148, -t155 * t65 - t149, 0, (-t95 * qJD(5) + t73 * t42 - t76 * t43) * t65 + t96 * t64 + t61 * t7 + t9 * t40 + t103 * t35 + t87 * t63 + t158 * t155, -(t96 * qJD(5) - t76 * t42 - t73 * t43) * t65 - t95 * t64 + t61 * t6 + t9 * t41 - t103 * t37 - t87 * t62 - t158 * t90; 0, 0, 0, 0, 0, -t74 * t80 * t77, t137 * t80, t117, t115, qJDD(4), -g(3) * t33 - t152 * t77 + t86 * t74, g(3) * t34 + t152 * t74 + t86 * t77, -t150, t8, t4, t5, t64, -(-t73 * t12 - t139) * t65 + (-t65 * t127 - t35 * t132 + t76 * t64) * pkin(4) + t81, (-qJD(5) * t11 + t12 * t65 - t3) * t76 + (-t65 * t126 + t37 * t132 - t73 * t64) * pkin(4) + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t8, t4, t5, t64, -t99 * t65 + t81, (-t3 + (-qJD(5) + t65) * t11) * t76 + t83;];
tau_reg = t1;
