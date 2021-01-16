% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:45
% EndTime: 2021-01-15 21:02:55
% DurationCPUTime: 2.34s
% Computational Cost: add. (2151->345), mult. (4460->426), div. (0->0), fcn. (2592->6), ass. (0->198)
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t116 = sin(qJ(2));
t180 = qJD(1) * qJD(2);
t166 = t116 * t180;
t119 = cos(qJ(2));
t178 = t119 * qJDD(1);
t248 = -t166 + t178;
t191 = qJD(2) * t115;
t192 = qJD(1) * t119;
t63 = t118 * t192 + t191;
t21 = t63 * qJD(4) - t118 * qJDD(2) + t248 * t115;
t193 = qJD(1) * t116;
t90 = qJD(4) + t193;
t225 = t63 * t90;
t250 = t21 - t225;
t249 = t21 + t225;
t235 = pkin(3) + pkin(6);
t121 = -pkin(2) - pkin(7);
t204 = t116 * qJ(3);
t59 = t121 * t119 - pkin(1) - t204;
t34 = t59 * qJD(1);
t100 = pkin(6) * t193;
t241 = qJD(3) + t100;
t195 = pkin(3) * t193 + t241;
t38 = t121 * qJD(2) + t195;
t15 = t115 * t38 + t118 * t34;
t211 = qJ(3) * t119;
t147 = pkin(7) * t116 - t211;
t182 = t116 * qJD(3);
t131 = t147 * qJD(2) - t182;
t88 = pkin(2) * t166;
t13 = t131 * qJD(1) + t59 * qJDD(1) + t88;
t165 = t119 * t180;
t179 = t116 * qJDD(1);
t136 = t165 + t179;
t87 = pkin(6) * t165;
t97 = pkin(6) * t179;
t171 = qJDD(3) + t87 + t97;
t24 = t136 * pkin(3) + t121 * qJDD(2) + t171;
t163 = -t115 * t13 + t118 * t24;
t129 = -t15 * qJD(4) + t163;
t213 = t21 * qJ(5);
t62 = qJDD(4) + t136;
t230 = t62 * pkin(4);
t167 = t115 * t192;
t189 = qJD(2) * t118;
t65 = -t167 + t189;
t1 = -t65 * qJD(5) + t129 + t213 + t230;
t120 = cos(qJ(1));
t201 = t116 * t120;
t117 = sin(qJ(1));
t203 = t116 * t117;
t177 = g(1) * t201 + g(2) * t203 - g(3) * t119;
t186 = qJD(4) * t118;
t187 = qJD(4) * t115;
t161 = -t115 * t24 - t118 * t13 - t38 * t186 + t34 * t187;
t22 = -qJD(4) * t167 + t115 * qJDD(2) + (qJD(2) * qJD(4) + t248) * t118;
t218 = qJ(5) * t22;
t2 = -qJD(5) * t63 - t161 - t218;
t8 = -qJ(5) * t63 + t15;
t229 = t8 * t90;
t14 = -t115 * t34 + t118 * t38;
t7 = -qJ(5) * t65 + t14;
t6 = pkin(4) * t90 + t7;
t247 = -(t6 * t90 - t2) * t115 + (t1 + t229) * t118 - t177;
t224 = t65 * t90;
t245 = -t22 + t224;
t244 = t22 + t224;
t47 = t118 * t62;
t243 = -t90 * t187 + t47;
t242 = qJ(5) - t121;
t240 = t22 * pkin(4) + qJDD(5);
t198 = t118 * t119;
t197 = t118 * t120;
t51 = -t115 * t117 + t116 * t197;
t199 = t117 * t118;
t53 = t115 * t120 + t116 * t199;
t239 = -g(1) * t51 - g(2) * t53 + g(3) * t198;
t112 = qJD(2) * qJ(3);
t101 = pkin(6) * t192;
t72 = pkin(3) * t192 + t101;
t49 = t112 + t72;
t237 = t121 * t62 + t49 * t90;
t236 = t65 ^ 2;
t234 = -t7 + t6;
t228 = g(1) * t117;
t227 = g(2) * t120;
t226 = g(3) * t116;
t104 = pkin(2) * t193;
t43 = t147 * qJD(1) + t104;
t223 = t115 * t72 + t118 * t43;
t181 = t118 * qJD(5);
t206 = t115 * t116;
t57 = t118 * t72;
t222 = t242 * t187 - t181 + t115 * t43 - t57 - (pkin(4) * t119 - qJ(5) * t206) * qJD(1);
t75 = t242 * t118;
t221 = -qJ(5) * t118 * t193 - qJD(4) * t75 - t115 * qJD(5) - t223;
t80 = t235 * t116;
t220 = t115 * t80 + t118 * t59;
t170 = -pkin(4) * t118 - pkin(3);
t219 = pkin(4) * t186 - t170 * t193 + t241;
t217 = t115 * t62;
t216 = t118 * t21;
t215 = t118 * t65;
t212 = pkin(6) * qJDD(2);
t149 = pkin(2) * t119 + t204;
t77 = pkin(1) + t149;
t210 = qJD(1) * t77;
t209 = qJD(2) * t63;
t208 = qJD(2) * t65;
t207 = qJDD(2) * pkin(2);
t205 = t115 * t119;
t202 = t116 * t118;
t123 = qJD(1) ^ 2;
t200 = t116 * t123;
t81 = t235 * t119;
t113 = t116 ^ 2;
t114 = t119 ^ 2;
t194 = t113 - t114;
t190 = qJD(2) * t116;
t188 = qJD(2) * t119;
t185 = qJD(4) * t119;
t184 = qJD(4) * t121;
t183 = qJDD(1) * t77;
t110 = qJDD(2) * qJ(3);
t111 = qJD(2) * qJD(3);
t98 = pkin(6) * t178;
t176 = t110 + t111 + t98;
t164 = -pkin(4) * t63 - qJD(5);
t26 = -t164 + t49;
t175 = t26 * t187;
t173 = t26 * t186;
t172 = t119 * t200;
t169 = t235 * qJD(2);
t168 = t115 * t185;
t96 = pkin(4) * t115 + qJ(3);
t162 = t116 * t96 + t119 * t242;
t160 = qJ(5) * t119 - t59;
t159 = -t97 + t177;
t158 = -qJD(2) * pkin(2) + qJD(3);
t157 = pkin(3) * t178 + t176;
t95 = pkin(6) - t170;
t156 = qJD(1) * t169;
t155 = g(1) * t53 - g(2) * t51;
t52 = t115 * t201 + t199;
t54 = -t115 * t203 + t197;
t154 = -g(1) * t54 - g(2) * t52;
t122 = qJD(2) ^ 2;
t153 = pkin(6) * t122 + t227;
t152 = g(1) * t120 + g(2) * t117;
t151 = -t227 + t228;
t148 = -pkin(2) * t116 + t211;
t76 = t100 + t158;
t78 = -t101 - t112;
t146 = t116 * t78 + t119 * t76;
t145 = t115 * t90;
t144 = -t193 * t210 + qJDD(3) - t159;
t141 = -t90 * t186 - t217;
t140 = t152 * t119;
t139 = -0.2e1 * pkin(1) * t180 - t212;
t103 = pkin(2) * t190;
t31 = t103 + t131;
t73 = t235 * t188;
t138 = t115 * t73 + t118 * t31 + t80 * t186 - t59 * t187;
t137 = -qJ(3) * t188 - t182;
t135 = 0.2e1 * qJDD(1) * pkin(1) - t153;
t134 = 0.2e1 * t210 * qJD(2) + t212;
t133 = -t140 - t226;
t25 = -t116 * t156 + t157;
t5 = t25 + t240;
t132 = t5 + t133;
t130 = t25 + t133;
t23 = t137 * qJD(1) - t183 + t88;
t46 = t103 + t137;
t128 = qJD(1) * t46 + t153 - t183 + t23;
t127 = g(1) * t52 - g(2) * t54 - g(3) * t205 + t161;
t125 = t129 + t239;
t35 = pkin(6) * t166 - t176;
t41 = t171 - t207;
t124 = t146 * qJD(2) + t41 * t116 - t35 * t119 - t152;
t93 = t119 * t228;
t74 = t242 * t115;
t71 = t116 * t169;
t69 = -qJ(3) * t192 + t104;
t68 = t118 * t80;
t61 = t63 ^ 2;
t58 = t118 * t73;
t48 = pkin(4) * t198 + t81;
t42 = pkin(1) + t162;
t27 = -pkin(4) * t168 - t95 * t190;
t18 = -qJ(5) * t198 + t220;
t17 = t116 * pkin(4) + t160 * t115 + t68;
t10 = -t90 ^ 2 * t118 - t208 - t217;
t9 = -t90 * t145 - t209 + t47;
t4 = -t119 * t181 + (t116 * t189 + t168) * qJ(5) + t138;
t3 = pkin(4) * t188 + t58 + t160 * t186 + (-qJ(5) * t190 - qJD(4) * t80 + qJD(5) * t119 - t31) * t115;
t11 = [qJDD(1), t151, t152, qJDD(1) * t113 + 0.2e1 * t116 * t165, 0.2e1 * t116 * t178 - 0.2e1 * t194 * t180, qJDD(2) * t116 + t119 * t122, qJDD(2) * t119 - t116 * t122, 0, t139 * t116 + t135 * t119 + t93, t139 * t119 + (-t135 - t228) * t116, (t113 + t114) * qJDD(1) * pkin(6) + t124, t116 * t134 + t128 * t119 - t93, t134 * t119 + (-t128 + t228) * t116, -t210 * t46 + (t151 - t23) * t77 + t124 * pkin(6), -t185 * t215 + (t119 * t21 + t65 * t190) * t115, (-t115 * t63 + t215) * t190 + (t115 * t22 + t216 + (t115 * t65 + t118 * t63) * qJD(4)) * t119, (t90 * t191 - t21) * t116 + (t141 + t208) * t119, (t90 * t189 - t22) * t116 + (-t209 - t243) * t119, t116 * t62 + t90 * t188, (-t115 * t31 + t58) * t90 + (-t115 * t59 + t68) * t62 + t163 * t116 - t71 * t63 + t81 * t22 + t25 * t198 + (t119 * t14 - t202 * t49) * qJD(2) + (-t15 * t116 - t49 * t205 - t220 * t90) * qJD(4) + t154, -t138 * t90 - t220 * t62 - t71 * t65 - t81 * t21 + (t191 * t49 + t161) * t116 + (-qJD(2) * t15 - t25 * t115 - t186 * t49) * t119 + t155, t17 * t62 + t48 * t22 + t27 * t63 + t3 * t90 + (-t189 * t26 + t1) * t116 + (qJD(2) * t6 + t118 * t5 - t175) * t119 + t154, -t18 * t62 - t48 * t21 + t27 * t65 - t4 * t90 + (t191 * t26 - t2) * t116 + (-qJD(2) * t8 - t115 * t5 - t173) * t119 + t155, t17 * t21 - t18 * t22 - t3 * t65 - t4 * t63 + t93 + (-t115 * t6 + t118 * t8) * t190 + (-t227 + t1 * t115 - t118 * t2 + (t115 * t8 + t118 * t6) * qJD(4)) * t119, t2 * t18 + t8 * t4 + t1 * t17 + t6 * t3 + t5 * t48 + t26 * t27 - g(1) * (-t117 * t42 + t120 * t95) - g(2) * (t117 * t95 + t120 * t42); 0, 0, 0, -t172, t194 * t123, t179, t178, qJDD(2), pkin(1) * t200 + t159, t226 - t98 + (pkin(1) * t123 + t152) * t119, t148 * qJDD(1) + ((-t78 - t112) * t116 + (t158 - t76) * t119) * qJD(1), -t69 * t192 + t144 - 0.2e1 * t207, 0.2e1 * t110 + 0.2e1 * t111 + t98 + (qJD(1) * t69 - g(3)) * t116 + (-qJD(1) * t210 - t152) * t119, -t146 * qJD(1) * pkin(6) - t41 * pkin(2) - g(3) * t149 - t35 * qJ(3) - t78 * qJD(3) - t152 * t148 + t210 * t69, -t65 * t145 - t216, t249 * t115 - t244 * t118, (-t119 * t65 - t90 * t206) * qJD(1) + t243, (t119 * t63 - t90 * t202) * qJD(1) + t141, -t90 * t192, -t14 * t192 + qJ(3) * t22 - t57 * t90 + t195 * t63 + t237 * t118 + ((t43 - t184) * t90 + t130) * t115, -qJ(3) * t21 + t223 * t90 + t15 * t192 + t195 * t65 - t237 * t115 + (-t184 * t90 + t130) * t118, t173 + t22 * t96 - t62 * t75 + t222 * t90 + t219 * t63 + (-t119 * t6 + t202 * t26) * qJD(1) + t132 * t115, -t175 - t21 * t96 + t62 * t74 - t221 * t90 + t219 * t65 + (t119 * t8 - t206 * t26) * qJD(1) + t132 * t118, -t21 * t75 + t22 * t74 - t221 * t63 - t222 * t65 - t247, -t2 * t74 - t1 * t75 + t5 * t96 - g(3) * t162 + t221 * t8 + t222 * t6 - t152 * (-t116 * t242 + t119 * t96) + t219 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, qJDD(2) + t172, -t113 * t123 - t122, qJD(2) * t78 + t144 - t207 + t87, 0, 0, 0, 0, 0, t9, t10, t9, t10, t245 * t115 + t250 * t118, -t26 * qJD(2) + t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t61 + t236, -t250, t245, t62, t15 * t90 - t49 * t65 + t125, t14 * t90 + t49 * t63 + t127, 0.2e1 * t230 + t213 + t229 + (t164 - t26) * t65 + t125, -pkin(4) * t236 + t218 + t7 * t90 + (qJD(5) + t26) * t63 + t127, t21 * pkin(4) - t234 * t63, t234 * t8 + (-t26 * t65 + t1 + t239) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t249, -t61 - t236, t6 * t65 + t63 * t8 - t140 + (-g(3) - t156) * t116 + t157 + t240;];
tau_reg = t11;
