% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:49
% DurationCPUTime: 1.63s
% Computational Cost: add. (4449->264), mult. (10303->367), div. (0->0), fcn. (6909->8), ass. (0->156)
t128 = sin(pkin(7));
t129 = cos(pkin(7));
t135 = cos(qJ(2));
t163 = qJD(1) * t135;
t132 = sin(qJ(2));
t164 = qJD(1) * t132;
t102 = t128 * t164 - t129 * t163;
t104 = t128 * t163 + t129 * t164;
t81 = t104 * t102;
t186 = qJDD(2) - t81;
t191 = t128 * t186;
t190 = t129 * t186;
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t87 = -t134 * qJD(2) + t131 * t104;
t89 = t131 * qJD(2) + t134 * t104;
t67 = t89 * t87;
t118 = t132 * qJDD(1);
t158 = qJD(1) * qJD(2);
t153 = t135 * t158;
t109 = t118 + t153;
t120 = t135 * qJDD(1);
t154 = t132 * t158;
t110 = t120 - t154;
t147 = t128 * t109 - t129 * t110;
t80 = qJDD(4) + t147;
t187 = -t67 + t80;
t189 = t131 * t187;
t188 = t134 * t187;
t162 = qJD(2) * t104;
t68 = t147 + t162;
t83 = t129 * t109 + t128 * t110;
t97 = qJD(2) * t102;
t70 = t83 - t97;
t137 = qJD(1) ^ 2;
t133 = sin(qJ(1));
t184 = cos(qJ(1));
t144 = t184 * g(1) + t133 * g(2);
t170 = qJDD(1) * pkin(5);
t140 = -t137 * pkin(1) - t144 + t170;
t139 = t132 * t140;
t168 = t132 * t137;
t138 = -t139 - t109 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t168 + qJ(3) * t158 - g(3)) * t135;
t126 = t135 ^ 2;
t123 = t126 * t137;
t142 = qJD(2) * pkin(2) - qJ(3) * t164;
t91 = -t132 * g(3) + t135 * t140;
t63 = -pkin(2) * t123 + t110 * qJ(3) - qJD(2) * t142 + t91;
t34 = -0.2e1 * qJD(3) * t102 + t128 * t138 + t129 * t63;
t148 = -t134 * qJDD(2) + t131 * t83;
t99 = qJD(4) + t102;
t39 = (qJD(4) - t99) * t89 + t148;
t151 = t133 * g(1) - t184 * g(2);
t143 = qJDD(1) * pkin(1) + t151;
t65 = t110 * pkin(2) + (qJ(3) * t126 + pkin(5)) * t137 - t142 * t164 - qJDD(3) + t143;
t85 = t87 ^ 2;
t86 = t89 ^ 2;
t98 = t99 ^ 2;
t100 = t102 ^ 2;
t101 = t104 ^ 2;
t185 = 0.2e1 * qJD(3);
t183 = pkin(3) * t128;
t181 = t128 * t65;
t78 = qJDD(2) + t81;
t180 = t128 * t78;
t179 = t129 * t65;
t178 = t129 * t78;
t136 = qJD(2) ^ 2;
t150 = t128 * t63 - t129 * t138;
t75 = t102 * pkin(3) - t104 * pkin(6);
t24 = -qJDD(2) * pkin(3) - t136 * pkin(6) + (t185 + t75) * t104 + t150;
t177 = t131 * t24;
t49 = t67 + t80;
t176 = t131 * t49;
t175 = t131 * t99;
t33 = t104 * t185 + t150;
t15 = t128 * t34 - t129 * t33;
t174 = t132 * t15;
t173 = t134 * t24;
t172 = t134 * t49;
t171 = t134 * t99;
t117 = t135 * t168;
t169 = t132 * (qJDD(2) + t117);
t167 = t135 * (qJDD(2) - t117);
t165 = qJD(4) + t99;
t161 = qJD(2) * t128;
t160 = qJD(2) * t129;
t157 = t128 * t67;
t156 = t129 * t67;
t155 = -pkin(3) * t129 - pkin(2);
t25 = -t136 * pkin(3) + qJDD(2) * pkin(6) - t102 * t75 + t34;
t29 = t68 * pkin(3) - t70 * pkin(6) - t65;
t10 = t131 * t25 - t134 * t29;
t11 = t131 * t29 + t134 * t25;
t5 = t131 * t10 + t134 * t11;
t16 = t128 * t33 + t129 * t34;
t90 = t135 * g(3) + t139;
t149 = t132 * t90 + t135 * t91;
t4 = -t134 * t10 + t131 * t11;
t145 = -t131 * qJDD(2) - t134 * t83;
t69 = -t147 + t162;
t56 = -t87 * qJD(4) - t145;
t125 = t132 ^ 2;
t121 = t125 * t137;
t111 = t120 - 0.2e1 * t154;
t108 = t118 + 0.2e1 * t153;
t106 = t137 * pkin(5) + t143;
t94 = -t101 - t136;
t93 = -t101 + t136;
t92 = t100 - t136;
t76 = -t136 - t100;
t74 = t99 * t87;
t73 = -t86 + t98;
t72 = t85 - t98;
t71 = t83 + t97;
t66 = -t100 - t101;
t64 = t86 - t85;
t59 = -t86 - t98;
t58 = -t128 * t94 - t178;
t57 = t129 * t94 - t180;
t55 = -t89 * qJD(4) - t148;
t54 = -t98 - t85;
t53 = t85 + t86;
t52 = t129 * t76 - t191;
t51 = t128 * t76 + t190;
t47 = (t131 * t89 - t134 * t87) * t99;
t46 = t128 * t71 + t129 * t69;
t45 = t128 * t69 - t129 * t71;
t44 = t165 * t87 + t145;
t43 = t56 + t74;
t42 = t56 - t74;
t40 = -t165 * t89 - t148;
t38 = t134 * t56 - t89 * t175;
t37 = -t131 * t55 + t87 * t171;
t36 = t134 * t72 - t176;
t35 = -t131 * t73 + t188;
t31 = -t131 * t59 - t172;
t30 = t134 * t59 - t176;
t27 = t134 * t54 - t189;
t26 = t131 * t54 + t188;
t23 = t131 * t43 - t134 * t39;
t22 = -t131 * t42 + t134 * t40;
t21 = -t131 * t39 - t134 * t43;
t20 = -t128 * t44 + t129 * t31;
t19 = t128 * t31 + t129 * t44;
t18 = -t128 * t40 + t129 * t27;
t17 = t128 * t27 + t129 * t40;
t14 = -t128 * t53 + t129 * t23;
t13 = t128 * t23 + t129 * t53;
t12 = -pkin(6) * t30 + t173;
t8 = -pkin(6) * t26 + t177;
t7 = -pkin(3) * t30 + t11;
t6 = -pkin(3) * t26 + t10;
t2 = t128 * t5 - t129 * t24;
t1 = -pkin(6) * t21 - t4;
t3 = [0, 0, 0, 0, 0, qJDD(1), t151, t144, 0, 0, (t109 + t153) * t132, t135 * t108 + t132 * t111, t169 + t135 * (-t121 + t136), (t110 - t154) * t135, t132 * (t123 - t136) + t167, 0, t135 * t106 + pkin(1) * t111 + pkin(5) * (t135 * (-t123 - t136) - t169), -t132 * t106 - pkin(1) * t108 + pkin(5) * (-t167 - t132 * (-t121 - t136)), pkin(1) * (t121 + t123) + (t125 + t126) * t170 + t149, pkin(1) * t106 + pkin(5) * t149, t132 * (-t104 * t161 + t129 * t83) + t135 * (t104 * t160 + t128 * t83), t132 * (-t128 * t70 - t129 * t68) + t135 * (-t128 * t68 + t129 * t70), t132 * (-t128 * t93 + t190) + t135 * (t129 * t93 + t191), t132 * (t102 * t160 + t128 * t147) + t135 * (t102 * t161 - t129 * t147), t132 * (t129 * t92 - t180) + t135 * (t128 * t92 + t178), (t132 * (-t102 * t129 + t104 * t128) + t135 * (-t102 * t128 - t104 * t129)) * qJD(2), t132 * (-qJ(3) * t51 - t181) + t135 * (-pkin(2) * t68 + qJ(3) * t52 + t179) - pkin(1) * t68 + pkin(5) * (-t132 * t51 + t135 * t52), t132 * (-qJ(3) * t57 - t179) + t135 * (-pkin(2) * t70 + qJ(3) * t58 - t181) - pkin(1) * t70 + pkin(5) * (-t132 * t57 + t135 * t58), t132 * (-qJ(3) * t45 - t15) + t135 * (-pkin(2) * t66 + qJ(3) * t46 + t16) - pkin(1) * t66 + pkin(5) * (-t132 * t45 + t135 * t46), -qJ(3) * t174 + t135 * (pkin(2) * t65 + qJ(3) * t16) + pkin(1) * t65 + pkin(5) * (t135 * t16 - t174), t132 * (t129 * t38 + t157) + t135 * (t128 * t38 - t156), t132 * (t128 * t64 + t129 * t22) + t135 * (t128 * t22 - t129 * t64), t132 * (t128 * t43 + t129 * t35) + t135 * (t128 * t35 - t129 * t43), t132 * (t129 * t37 - t157) + t135 * (t128 * t37 + t156), t132 * (-t128 * t39 + t129 * t36) + t135 * (t128 * t36 + t129 * t39), t132 * (t128 * t80 + t129 * t47) + t135 * (t128 * t47 - t129 * t80), t132 * (-qJ(3) * t17 - t128 * t6 + t129 * t8) + t135 * (-pkin(2) * t26 + qJ(3) * t18 + t128 * t8 + t129 * t6) - pkin(1) * t26 + pkin(5) * (-t132 * t17 + t135 * t18), t132 * (-qJ(3) * t19 + t129 * t12 - t128 * t7) + t135 * (-pkin(2) * t30 + qJ(3) * t20 + t128 * t12 + t129 * t7) - pkin(1) * t30 + pkin(5) * (-t132 * t19 + t135 * t20), t132 * (-qJ(3) * t13 + t129 * t1) + t135 * (qJ(3) * t14 + t128 * t1) + pkin(5) * (-t132 * t13 + t135 * t14) + (t132 * t183 + t135 * t155 - pkin(1)) * t21, (t132 * (-pkin(6) * t129 + t183) + t135 * (-pkin(6) * t128 + t155) - pkin(1)) * t4 + (pkin(5) + qJ(3)) * (-t132 * t2 + t135 * (t128 * t24 + t129 * t5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t121 - t123, t118, t117, t120, qJDD(2), -t90, -t91, 0, 0, t81, t101 - t100, t71, -t81, t69, qJDD(2), pkin(2) * t51 - t33, pkin(2) * t57 - t34, pkin(2) * t45, pkin(2) * t15, t131 * t56 + t89 * t171, t131 * t40 + t134 * t42, t134 * t73 + t189, t134 * t55 + t87 * t175, t131 * t72 + t172, (-t131 * t87 - t134 * t89) * t99, pkin(2) * t17 + pkin(3) * t40 + pkin(6) * t27 - t173, pkin(2) * t19 + pkin(3) * t44 + pkin(6) * t31 + t177, pkin(2) * t13 + pkin(3) * t53 + pkin(6) * t23 + t5, pkin(2) * t2 - pkin(3) * t24 + pkin(6) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t70, t66, -t65, 0, 0, 0, 0, 0, 0, t26, t30, t21, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t64, t43, -t67, -t39, t80, -t10, -t11, 0, 0;];
tauJ_reg = t3;
