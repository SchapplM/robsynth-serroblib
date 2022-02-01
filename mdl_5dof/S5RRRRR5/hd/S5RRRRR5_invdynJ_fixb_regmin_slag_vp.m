% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:12
% DurationCPUTime: 1.56s
% Computational Cost: add. (2163->243), mult. (2960->323), div. (0->0), fcn. (1907->16), ass. (0->175)
t128 = sin(qJ(3));
t133 = cos(qJ(3));
t134 = cos(qJ(2));
t191 = qJD(1) * t134;
t172 = qJD(2) * t191;
t129 = sin(qJ(2));
t185 = qJDD(1) * t129;
t143 = (t172 + t185) * pkin(1);
t217 = t134 * pkin(1);
t107 = qJDD(1) * t217;
t119 = qJDD(1) + qJDD(2);
t210 = pkin(1) * qJD(1);
t183 = t129 * t210;
t50 = t119 * pkin(2) - qJD(2) * t183 + t107;
t121 = qJD(1) + qJD(2);
t180 = pkin(1) * t191;
t71 = t121 * pkin(2) + t180;
t166 = qJD(3) * t183;
t84 = t128 * t166;
t230 = -t128 * t50 - (qJD(3) * t71 + t143) * t133 + t84;
t111 = qJD(3) + t121;
t126 = sin(qJ(5));
t132 = cos(qJ(4));
t127 = sin(qJ(4));
t131 = cos(qJ(5));
t195 = t131 * t127;
t59 = t126 * t132 + t195;
t46 = t59 * t111;
t190 = qJD(3) * t128;
t197 = t128 * t129;
t94 = pkin(1) * t197;
t164 = t128 * pkin(1) * t172 + qJDD(1) * t94 + t71 * t190 + (t166 - t50) * t133;
t110 = qJDD(3) + t119;
t221 = t110 * pkin(3);
t13 = t164 - t221;
t125 = qJ(1) + qJ(2);
t116 = qJ(3) + t125;
t102 = cos(t116);
t222 = g(2) * t102;
t229 = t13 + t222;
t189 = qJD(3) * t133;
t90 = t128 * t183;
t54 = t133 * t180 - t90;
t228 = pkin(2) * t189 - t54;
t101 = sin(t116);
t226 = g(1) * t102 + g(2) * t101;
t196 = t129 * t133;
t156 = t128 * t134 + t196;
t53 = t156 * t210;
t163 = pkin(2) * t190 - t53;
t120 = qJD(4) + qJD(5);
t103 = t128 * pkin(2) + pkin(8);
t218 = t133 * pkin(2);
t104 = -pkin(3) - t218;
t136 = qJD(4) ^ 2;
t225 = -t103 * t136 - t104 * t110 - t111 * t163;
t224 = -pkin(8) - pkin(9);
t105 = pkin(2) + t217;
t212 = pkin(1) * t196 + t128 * t105;
t52 = pkin(8) + t212;
t223 = -pkin(9) - t52;
t92 = g(1) * t101;
t220 = t111 * pkin(3);
t219 = t132 * pkin(4);
t194 = t131 * t132;
t178 = t111 * t194;
t200 = t126 * t127;
t179 = t111 * t200;
t44 = -t178 + t179;
t216 = t46 * t44;
t215 = -pkin(9) - t103;
t188 = qJD(4) * t127;
t42 = t133 * t71 - t90;
t36 = -t42 - t220;
t214 = t132 * t92 + t36 * t188;
t108 = pkin(4) * t188;
t213 = t108 + t163;
t113 = sin(t125);
t115 = cos(t125);
t211 = g(1) * t115 + g(2) * t113;
t43 = t128 * t71 + t133 * t183;
t37 = t111 * pkin(8) + t43;
t175 = pkin(9) * t111 + t37;
t26 = t175 * t132;
t209 = t131 * t26;
t30 = t105 * t190 + (t156 * qJD(2) + t129 * t189) * pkin(1);
t208 = t30 * t111;
t207 = t43 * t111;
t124 = qJ(4) + qJ(5);
t112 = sin(t124);
t205 = t101 * t112;
t114 = cos(t124);
t204 = t101 * t114;
t203 = t102 * t112;
t202 = t102 * t114;
t201 = t111 * t127;
t198 = t127 * t110;
t193 = t132 * t110;
t122 = t127 ^ 2;
t192 = -t132 ^ 2 + t122;
t187 = qJD(4) * t132;
t186 = qJD(5) * t126;
t184 = t229 * t127 + t36 * t187;
t106 = -pkin(3) - t219;
t177 = qJD(4) * t224;
t176 = t111 * t187;
t173 = qJD(4) * t223;
t170 = qJD(4) * t215;
t169 = t133 * t105 - t94;
t168 = qJD(1) * (-qJD(2) + t121);
t167 = qJD(2) * (-qJD(1) - t121);
t51 = -pkin(3) - t169;
t165 = g(1) * t113 - g(2) * t115 + t107;
t162 = -t43 + t108;
t161 = t126 * t198 - t131 * t193;
t25 = t175 * t127;
t24 = qJD(4) * pkin(4) - t25;
t160 = -t126 * t24 - t209;
t38 = t223 * t127;
t117 = t132 * pkin(9);
t39 = t132 * t52 + t117;
t159 = -t126 * t39 + t131 * t38;
t158 = t126 * t38 + t131 * t39;
t58 = -t194 + t200;
t149 = t111 * t188 - t193;
t10 = t149 * pkin(4) + t13;
t28 = t106 * t111 - t42;
t34 = t120 * t58;
t155 = -g(1) * t205 + g(2) * t203 + t10 * t59 - t28 * t34;
t35 = t120 * t59;
t154 = g(1) * t204 - g(2) * t202 + t10 * t58 + t28 * t35;
t57 = t132 * t103 + t117;
t153 = qJD(5) * t57 + t228 * t127 - t132 * t170;
t87 = t132 * pkin(8) + t117;
t152 = qJD(5) * t87 - t127 * t42 - t132 * t177;
t56 = t215 * t127;
t151 = -qJD(5) * t56 - t127 * t170 - t228 * t132;
t86 = t224 * t127;
t150 = -qJD(5) * t86 - t127 * t177 + t132 * t42;
t147 = pkin(8) * t136 - t207 - t221;
t146 = t110 * t51 + t136 * t52 + t208;
t145 = -t164 + t92 - t222;
t12 = t110 * pkin(8) - t230;
t144 = -t111 * t36 - t12 + t226;
t142 = -pkin(8) * qJDD(4) + (t42 - t220) * qJD(4);
t29 = t105 * t189 + (-t129 * t190 + (t133 * t134 - t197) * qJD(2)) * pkin(1);
t141 = -qJDD(4) * t52 + (t111 * t51 - t29) * qJD(4);
t14 = qJD(5) * t178 + t110 * t195 - t120 * t179 + t126 * t193 + t131 * t176;
t4 = -t37 * t187 + qJDD(4) * pkin(4) - t127 * t12 + (-t176 - t198) * pkin(9);
t140 = t28 * t44 + t26 * t186 + g(2) * t204 + g(1) * t202 + g(3) * t112 + (-t26 * t120 - t4) * t126;
t139 = -qJDD(4) * t103 + (t104 * t111 - t228) * qJD(4);
t5 = -t149 * pkin(9) + t132 * t12 - t37 * t188;
t138 = g(1) * t203 + g(2) * t205 - g(3) * t114 + t160 * qJD(5) - t126 * t5 + t131 * t4 - t28 * t46;
t137 = t226 + t230;
t135 = cos(qJ(1));
t130 = sin(qJ(1));
t118 = qJDD(4) + qJDD(5);
t109 = t111 ^ 2;
t80 = t106 - t218;
t79 = qJDD(4) * t132 - t136 * t127;
t78 = qJDD(4) * t127 + t136 * t132;
t49 = t122 * t110 + 0.2e1 * t127 * t176;
t48 = t51 - t219;
t33 = -0.2e1 * t192 * t111 * qJD(4) + 0.2e1 * t127 * t193;
t27 = t108 + t30;
t22 = -t58 * t118 - t35 * t120;
t21 = t59 * t118 - t34 * t120;
t20 = -t127 * t29 + t132 * t173;
t19 = t127 * t173 + t132 * t29;
t18 = -t44 ^ 2 + t46 ^ 2;
t15 = t35 * t111 + t161;
t8 = t44 * t120 + t14;
t2 = t14 * t59 - t46 * t34;
t1 = -t14 * t58 - t59 * t15 + t34 * t44 - t46 * t35;
t3 = [qJDD(1), g(1) * t130 - g(2) * t135, g(1) * t135 + g(2) * t130, t119, (t119 * t134 + t129 * t167) * pkin(1) + t165, ((-qJDD(1) - t119) * t129 + t134 * t167) * pkin(1) + t211, t110, t110 * t169 + t145 - t208, -t212 * t110 - t29 * t111 + t137, t49, t33, t78, t79, 0, t141 * t127 + (-t146 - t229) * t132 + t214, t141 * t132 + (t146 - t92) * t127 + t184, t2, t1, t21, t22, 0, t27 * t44 + t48 * t15 + (-qJD(5) * t158 - t126 * t19 + t131 * t20) * t120 + t159 * t118 + t154, t27 * t46 + t48 * t14 - (qJD(5) * t159 + t126 * t20 + t131 * t19) * t120 - t158 * t118 + t155; 0, 0, 0, t119, pkin(1) * t129 * t168 + t165, (t134 * t168 - t185) * pkin(1) + t211, t110, t53 * t111 + (t110 * t133 - t111 * t190) * pkin(2) + t145, t54 * t111 + t84 + (-pkin(2) * t110 - t50) * t128 + ((-pkin(2) * t111 - t71) * qJD(3) - t143) * t133 + t226, t49, t33, t78, t79, 0, t139 * t127 + (-t229 + t225) * t132 + t214, t139 * t132 + (-t225 - t92) * t127 + t184, t2, t1, t21, t22, 0, t80 * t15 + (-t126 * t57 + t131 * t56) * t118 + t213 * t44 + (t126 * t151 - t131 * t153) * t120 + t154, t80 * t14 - (t126 * t56 + t131 * t57) * t118 + t213 * t46 + (t126 * t153 + t131 * t151) * t120 + t155; 0, 0, 0, 0, 0, 0, t110, t145 + t207, t42 * t111 + t137, t49, t33, t78, t79, 0, t142 * t127 + (-t147 - t229) * t132 + t214, t142 * t132 + (t147 - t92) * t127 + t184, t2, t1, t21, t22, 0, t106 * t15 + (-t126 * t87 + t131 * t86) * t118 + t162 * t44 + (t126 * t150 - t131 * t152) * t120 + t154, t106 * t14 - (t126 * t86 + t131 * t87) * t118 + t162 * t46 + (t126 * t152 + t131 * t150) * t120 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t109 * t132, t192 * t109, t198, t193, qJDD(4), -g(3) * t132 + t127 * t144, g(3) * t127 + t132 * t144, t216, t18, t8, -t161, t118, -(t126 * t25 - t209) * t120 + (t131 * t118 - t120 * t186 - t44 * t201) * pkin(4) + t138, (-qJD(5) * t24 - t25 * t120 - t5) * t131 + (-qJD(5) * t131 * t120 - t126 * t118 - t46 * t201) * pkin(4) + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t18, t8, -t161, t118, -t120 * t160 + t138, (-t5 + (-qJD(5) + t120) * t24) * t131 + t140;];
tau_reg = t3;
