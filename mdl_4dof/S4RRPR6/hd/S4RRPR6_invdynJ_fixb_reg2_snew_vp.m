% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:59
% DurationCPUTime: 1.75s
% Computational Cost: add. (5430->269), mult. (12871->377), div. (0->0), fcn. (8790->8), ass. (0->163)
t194 = -2 * qJD(3);
t134 = sin(pkin(7));
t135 = cos(pkin(7));
t140 = cos(qJ(2));
t137 = sin(qJ(2));
t164 = qJD(1) * t137;
t109 = -qJD(1) * t135 * t140 + t134 * t164;
t111 = (t134 * t140 + t135 * t137) * qJD(1);
t94 = t111 * t109;
t187 = qJDD(2) - t94;
t193 = t134 * t187;
t192 = t135 * t187;
t136 = sin(qJ(4));
t130 = qJDD(2) + qJDD(4);
t139 = cos(qJ(4));
t85 = t109 * t139 + t111 * t136;
t87 = -t109 * t136 + t111 * t139;
t63 = t87 * t85;
t189 = -t63 + t130;
t191 = t136 * t189;
t190 = t139 * t189;
t106 = qJD(2) * t109;
t124 = t137 * qJDD(1);
t157 = qJD(1) * qJD(2);
t155 = t140 * t157;
t117 = t124 + t155;
t125 = t140 * qJDD(1);
t156 = t137 * t157;
t118 = t125 - t156;
t96 = t117 * t135 + t118 * t134;
t188 = -t96 - t106;
t142 = qJD(1) ^ 2;
t166 = t137 * t142;
t138 = sin(qJ(1));
t185 = cos(qJ(1));
t147 = g(1) * t185 + g(2) * t138;
t171 = qJDD(1) * pkin(5);
t114 = -pkin(1) * t142 - t147 + t171;
t168 = t137 * t114;
t71 = qJDD(2) * pkin(2) - t117 * qJ(3) - t168 + (pkin(2) * t166 + qJ(3) * t157 - g(3)) * t140;
t133 = t140 ^ 2;
t127 = t133 * t142;
t145 = qJD(2) * pkin(2) - qJ(3) * t164;
t99 = -t137 * g(3) + t140 * t114;
t72 = -pkin(2) * t127 + t118 * qJ(3) - qJD(2) * t145 + t99;
t150 = t111 * t194 - t134 * t72 + t135 * t71;
t46 = t109 * t194 + t134 * t71 + t135 * t72;
t186 = pkin(6) * t188 + t150;
t83 = t85 ^ 2;
t84 = t87 ^ 2;
t107 = t109 ^ 2;
t108 = t111 ^ 2;
t131 = qJD(2) + qJD(4);
t129 = t131 ^ 2;
t143 = pkin(3) * t187 + t186;
t149 = qJD(2) * pkin(3) - pkin(6) * t111;
t95 = -t117 * t134 + t118 * t135;
t28 = -t107 * pkin(3) + t95 * pkin(6) - qJD(2) * t149 + t46;
t12 = t136 * t28 - t139 * t143;
t174 = t139 * t28;
t13 = t136 * t143 + t174;
t6 = -t12 * t139 + t13 * t136;
t184 = t134 * t6;
t183 = t135 * t6;
t154 = t138 * g(1) - t185 * g(2);
t146 = qJDD(1) * pkin(1) + t154;
t73 = t118 * pkin(2) - qJDD(3) - t145 * t164 + (qJ(3) * t133 + pkin(5)) * t142 + t146;
t181 = t134 * t73;
t91 = qJDD(2) + t94;
t180 = t134 * t91;
t179 = t135 * t73;
t178 = t135 * t91;
t44 = t95 * pkin(3) + t107 * pkin(6) - t111 * t149 + t73;
t177 = t136 * t44;
t60 = t63 + t130;
t176 = t136 * t60;
t25 = t134 * t46 + t135 * t150;
t175 = t137 * t25;
t173 = t139 * t44;
t172 = t139 * t60;
t170 = t131 * t136;
t169 = t131 * t139;
t123 = t140 * t166;
t167 = t137 * (qJDD(2) + t123);
t165 = t140 * (qJDD(2) - t123);
t163 = qJD(2) * t111;
t162 = qJD(2) * t134;
t161 = qJD(2) * t135;
t158 = qJD(4) + t131;
t7 = t12 * t136 + t13 * t139;
t26 = -t134 * t150 + t135 * t46;
t152 = t136 * t96 - t139 * t95;
t98 = t140 * g(3) + t168;
t151 = t137 * t98 + t140 * t99;
t148 = t136 * t95 + t139 * t96;
t76 = t95 + t163;
t144 = (-qJD(4) + t131) * t87 - t152;
t53 = -qJD(4) * t85 + t148;
t141 = qJD(2) ^ 2;
t132 = t137 ^ 2;
t126 = t132 * t142;
t119 = t125 - 0.2e1 * t156;
t116 = t124 + 0.2e1 * t155;
t113 = pkin(5) * t142 + t146;
t102 = -t108 - t141;
t101 = -t108 + t141;
t100 = t107 - t141;
t89 = -t141 - t107;
t82 = t131 * t85;
t81 = -t84 + t129;
t80 = t83 - t129;
t79 = -t84 - t129;
t77 = -t106 + t96;
t75 = -t95 + t163;
t74 = -t107 - t108;
t67 = -t102 * t134 - t178;
t66 = t102 * t135 - t180;
t65 = t135 * t89 - t193;
t64 = t134 * t89 + t192;
t62 = t84 - t83;
t58 = -t129 - t83;
t57 = (t136 * t87 - t139 * t85) * t131;
t56 = (-t136 * t85 - t139 * t87) * t131;
t55 = -t134 * t188 + t135 * t76;
t54 = t134 * t76 + t135 * t188;
t52 = -qJD(4) * t87 - t152;
t51 = -t83 - t84;
t50 = t139 * t80 - t176;
t49 = -t136 * t81 + t190;
t48 = t136 * t80 + t172;
t47 = t139 * t81 + t191;
t43 = -t136 * t79 - t172;
t42 = t139 * t79 - t176;
t40 = t53 + t82;
t39 = t53 - t82;
t38 = -t158 * t85 + t148;
t35 = t158 * t87 + t152;
t34 = t139 * t53 - t170 * t87;
t33 = t136 * t53 + t169 * t87;
t32 = -t136 * t52 + t169 * t85;
t31 = t139 * t52 + t170 * t85;
t30 = t139 * t58 - t191;
t29 = t136 * t58 + t190;
t24 = -pkin(6) * t42 - t173;
t23 = -t134 * t42 + t135 * t43;
t22 = t134 * t43 + t135 * t42;
t21 = t136 * t40 + t139 * t144;
t20 = -t136 * t39 - t139 * t35;
t19 = t136 * t144 - t139 * t40;
t18 = -t136 * t35 + t139 * t39;
t17 = -pkin(6) * t29 - t177;
t16 = -t134 * t29 + t135 * t30;
t15 = t134 * t30 + t135 * t29;
t14 = -pkin(3) * t38 + pkin(6) * t43 - t177;
t10 = -pkin(3) * t35 + pkin(6) * t30 + t173;
t9 = -t134 * t19 + t135 * t21;
t8 = t134 * t21 + t135 * t19;
t5 = pkin(3) * t44 + pkin(6) * t7;
t4 = -pkin(6) * t19 - t6;
t3 = -pkin(3) * t51 + pkin(6) * t21 + t7;
t2 = t135 * t7 - t184;
t1 = t134 * t7 + t183;
t11 = [0, 0, 0, 0, 0, qJDD(1), t154, t147, 0, 0, (t117 + t155) * t137, t116 * t140 + t119 * t137, t167 + t140 * (-t126 + t141), (t118 - t156) * t140, t137 * (t127 - t141) + t165, 0, t140 * t113 + pkin(1) * t119 + pkin(5) * (t140 * (-t127 - t141) - t167), -t137 * t113 - pkin(1) * t116 + pkin(5) * (-t165 - t137 * (-t126 - t141)), pkin(1) * (t126 + t127) + (t132 + t133) * t171 + t151, pkin(1) * t113 + pkin(5) * t151, t137 * (-t111 * t162 + t135 * t96) + t140 * (t111 * t161 + t134 * t96), t137 * (-t134 * t77 - t135 * t75) + t140 * (-t134 * t75 + t135 * t77), t137 * (-t101 * t134 + t192) + t140 * (t101 * t135 + t193), t137 * (t109 * t161 - t134 * t95) + t140 * (t109 * t162 + t135 * t95), t137 * (t100 * t135 - t180) + t140 * (t100 * t134 + t178), (t137 * (-t109 * t135 + t111 * t134) + t140 * (-t109 * t134 - t111 * t135)) * qJD(2), t137 * (-qJ(3) * t64 - t181) + t140 * (-pkin(2) * t75 + qJ(3) * t65 + t179) - pkin(1) * t75 + pkin(5) * (-t137 * t64 + t140 * t65), t137 * (-qJ(3) * t66 - t179) + t140 * (-pkin(2) * t77 + qJ(3) * t67 - t181) - pkin(1) * t77 + pkin(5) * (-t137 * t66 + t140 * t67), t137 * (-qJ(3) * t54 - t25) + t140 * (-pkin(2) * t74 + qJ(3) * t55 + t26) - pkin(1) * t74 + pkin(5) * (-t137 * t54 + t140 * t55), -qJ(3) * t175 + t140 * (pkin(2) * t73 + qJ(3) * t26) + pkin(1) * t73 + pkin(5) * (t140 * t26 - t175), t137 * (-t134 * t33 + t135 * t34) + t140 * (t134 * t34 + t135 * t33), t137 * (-t134 * t18 + t135 * t20) + t140 * (t134 * t20 + t135 * t18), t137 * (-t134 * t47 + t135 * t49) + t140 * (t134 * t49 + t135 * t47), t137 * (-t134 * t31 + t135 * t32) + t140 * (t134 * t32 + t135 * t31), t137 * (-t134 * t48 + t135 * t50) + t140 * (t134 * t50 + t135 * t48), t137 * (-t134 * t56 + t135 * t57) + t140 * (t134 * t57 + t135 * t56), t137 * (-qJ(3) * t15 - t10 * t134 + t135 * t17) + t140 * (-pkin(2) * t35 + qJ(3) * t16 + t10 * t135 + t134 * t17) - pkin(1) * t35 + pkin(5) * (-t137 * t15 + t140 * t16), t137 * (-qJ(3) * t22 - t134 * t14 + t135 * t24) + t140 * (-pkin(2) * t38 + qJ(3) * t23 + t134 * t24 + t135 * t14) - pkin(1) * t38 + pkin(5) * (-t137 * t22 + t140 * t23), t137 * (-qJ(3) * t8 - t134 * t3 + t135 * t4) + t140 * (-pkin(2) * t51 + qJ(3) * t9 + t134 * t4 + t135 * t3) - pkin(1) * t51 + pkin(5) * (-t137 * t8 + t140 * t9), t137 * (-pkin(6) * t183 - qJ(3) * t1 - t134 * t5) + t140 * (pkin(2) * t44 - pkin(6) * t184 + qJ(3) * t2 + t135 * t5) + pkin(1) * t44 + pkin(5) * (-t1 * t137 + t140 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t126 - t127, t124, t123, t125, qJDD(2), -t98, -t99, 0, 0, t94, t108 - t107, -t188, -t94, t76, qJDD(2), pkin(2) * t64 + t150, pkin(2) * t66 - t46, pkin(2) * t54, pkin(2) * t25, t63, t62, t40, -t63, t144, t130, pkin(2) * t15 + pkin(3) * t29 - t12, -t174 - t136 * t186 + pkin(2) * t22 + (-t136 * t187 + t42) * pkin(3), pkin(2) * t8 + pkin(3) * t19, pkin(2) * t1 + pkin(3) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t77, t74, -t73, 0, 0, 0, 0, 0, 0, t35, t38, t51, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t62, t40, -t63, t144, t130, -t12, -t13, 0, 0;];
tauJ_reg = t11;
