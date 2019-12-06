% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:37
% EndTime: 2019-12-05 16:37:49
% DurationCPUTime: 3.35s
% Computational Cost: add. (2100->373), mult. (5169->556), div. (0->0), fcn. (4157->12), ass. (0->179)
t131 = sin(qJ(3));
t195 = qJD(2) * qJD(3);
t178 = t131 * t195;
t134 = cos(qJ(3));
t191 = t134 * qJDD(2);
t147 = t178 - t191;
t126 = sin(pkin(10));
t127 = sin(pkin(5));
t132 = sin(qJ(2));
t128 = cos(pkin(10));
t135 = cos(qJ(2));
t214 = t128 * t135;
t71 = (t126 * t132 + t134 * t214) * t127;
t162 = pkin(3) * t131 - qJ(4) * t134;
t80 = t162 * qJD(3) - qJD(4) * t131;
t237 = qJD(1) * t71 - t126 * t80;
t206 = qJD(1) * t127;
t186 = t135 * t206;
t187 = t132 * t206;
t218 = t126 * t134;
t236 = t186 * t218 + (-t187 + t80) * t128;
t100 = qJD(2) * pkin(7) + t187;
t129 = cos(pkin(5));
t213 = t129 * t134;
t235 = qJD(1) * t213 - t131 * t100;
t234 = -qJDD(3) * pkin(3) + qJDD(4);
t136 = qJD(3) ^ 2;
t196 = qJD(1) * qJD(2);
t180 = t132 * t196;
t216 = t127 * t135;
t161 = -qJDD(1) * t216 + t127 * t180;
t222 = sin(pkin(9));
t171 = t222 * t132;
t223 = cos(pkin(9));
t172 = t223 * t135;
t81 = -t129 * t172 + t171;
t231 = g(2) * t81;
t170 = t222 * t135;
t173 = t223 * t132;
t83 = t129 * t170 + t173;
t232 = g(1) * t83;
t169 = t231 + t232;
t233 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t136 + t127 * (-g(3) * t135 + t180) - t161 + t169;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t202 = qJD(2) * t134;
t182 = t130 * t202;
t179 = t134 * t195;
t192 = t131 * qJDD(2);
t148 = t179 + t192;
t193 = qJDD(3) * t126;
t70 = t128 * t148 + t193;
t201 = qJD(3) * t126;
t203 = qJD(2) * t131;
t97 = t128 * t203 + t201;
t12 = -qJD(5) * t182 + t130 * t70 + (qJD(5) * t97 - t147) * t133;
t205 = qJD(1) * t131;
t115 = t129 * t205;
t194 = qJDD(1) * t129;
t199 = qJD(3) * t134;
t75 = qJDD(2) * pkin(7) + (qJDD(1) * t132 + t135 * t196) * t127;
t149 = qJD(3) * t115 + t100 * t199 + t131 * t75 - t134 * t194;
t17 = t149 + t234;
t122 = t128 * qJDD(3);
t69 = t126 * t148 - t122;
t10 = pkin(4) * t69 - pkin(8) * t70 + t17;
t67 = t134 * t100 + t115;
t54 = qJD(3) * qJ(4) + t67;
t155 = pkin(3) * t134 + qJ(4) * t131 + pkin(2);
t68 = -t155 * qJD(2) - t186;
t20 = t126 * t68 + t128 * t54;
t16 = -pkin(8) * t202 + t20;
t51 = -qJD(3) * pkin(3) + qJD(4) - t235;
t123 = t128 * qJD(3);
t95 = t126 * t203 - t123;
t18 = pkin(4) * t95 - pkin(8) * t97 + t51;
t160 = t130 * t16 - t133 * t18;
t177 = t131 * t194;
t14 = t177 + qJDD(3) * qJ(4) + t134 * t75 + (qJD(4) + t235) * qJD(3);
t24 = t80 * qJD(2) - t155 * qJDD(2) + t161;
t8 = t126 * t24 + t128 * t14;
t4 = pkin(8) * t147 + t8;
t1 = -t160 * qJD(5) + t130 * t10 + t133 * t4;
t188 = pkin(7) * t126 + pkin(4);
t200 = qJD(3) * t131;
t230 = -t188 * t200 - t236;
t190 = pkin(7) * t200;
t229 = t126 * t190 + t236;
t228 = -t128 * t190 - t237;
t99 = t162 * qJD(2);
t32 = t126 * t99 + t128 * t235;
t227 = qJD(2) * pkin(2);
t63 = qJDD(5) + t69;
t225 = t130 * t63;
t224 = t133 * t63;
t215 = t128 * t134;
t73 = pkin(7) * t215 - t126 * t155;
t219 = t155 * t128;
t217 = t127 * t132;
t212 = t130 * t131;
t211 = t131 * t133;
t210 = t133 * t134;
t209 = -qJD(5) - t95;
t208 = qJDD(1) - g(3);
t124 = t131 ^ 2;
t207 = -t134 ^ 2 + t124;
t204 = qJD(2) * t127;
t198 = qJD(5) * t130;
t197 = qJD(5) * t133;
t189 = t131 * t216;
t184 = t132 * t204;
t183 = t135 * t204;
t181 = qJ(4) * t191;
t176 = t127 * t223;
t175 = t127 * t222;
t174 = t209 ^ 2;
t82 = t129 * t173 + t170;
t84 = -t129 * t171 + t172;
t168 = g(1) * t84 + g(2) * t82;
t165 = pkin(4) * t126 - pkin(8) * t128;
t158 = pkin(7) + t165;
t77 = t158 * t131;
t167 = -qJD(5) * t77 - (-pkin(7) * t128 + pkin(8)) * t200 + t237;
t102 = -pkin(4) * t128 - pkin(8) * t126 - pkin(3);
t163 = t102 * t63 + (t165 * t202 + t67) * t209;
t7 = -t126 * t14 + t128 * t24;
t19 = -t126 * t54 + t128 * t68;
t31 = -t126 * t235 + t128 * t99;
t6 = t130 * t18 + t133 * t16;
t87 = t129 * t131 + t134 * t217;
t44 = -t126 * t216 + t128 * t87;
t86 = t131 * t217 - t213;
t26 = t130 * t86 + t133 * t44;
t25 = -t130 * t44 + t133 * t86;
t159 = -qJ(4) * t63 + qJD(4) * t209;
t137 = qJD(2) ^ 2;
t156 = qJDD(2) * t135 - t132 * t137;
t58 = t130 * t97 + t133 * t202;
t153 = t128 * t210 + t212;
t89 = t128 * t212 + t210;
t45 = t82 * t131 + t134 * t176;
t47 = t131 * t84 - t134 * t175;
t151 = g(1) * t47 + g(2) * t45 + g(3) * t86;
t46 = -t131 * t176 + t134 * t82;
t48 = t131 * t175 + t84 * t134;
t150 = g(1) * t48 + g(2) * t46 + g(3) * t87;
t146 = t151 - t17;
t145 = g(3) * t216 - t169;
t144 = -qJ(4) * t200 + (qJD(4) - t51) * t134;
t57 = -pkin(8) * t134 + t73;
t143 = qJD(5) * t57 + t131 * t186 - t158 * t199;
t142 = -qJ(4) * qJD(5) * t209 - t151;
t2 = -t6 * qJD(5) + t133 * t10 - t130 * t4;
t101 = -t186 - t227;
t140 = -pkin(7) * qJDD(3) + (t101 + t186 - t227) * qJD(3);
t139 = -(pkin(8) * t203 - qJD(5) * t102 + t32) * t209 - t150;
t138 = -t149 + t151;
t11 = -qJD(5) * t58 + t147 * t130 + t133 * t70;
t90 = t128 * t211 - t130 * t134;
t79 = t153 * qJD(2);
t78 = t128 * t182 - t133 * t203;
t72 = -pkin(7) * t218 - t219;
t60 = t133 * t97 - t182;
t56 = t188 * t134 + t219;
t50 = qJD(3) * t87 + t131 * t183;
t49 = -qJD(3) * t86 + t134 * t183;
t43 = t126 * t87 + t127 * t214;
t37 = -t133 * t200 - t134 * t198 + (t130 * t199 + t131 * t197) * t128;
t36 = qJD(3) * t153 - qJD(5) * t89;
t34 = t126 * t84 - t83 * t215;
t33 = t126 * t82 - t81 * t215;
t30 = t126 * t184 + t128 * t49;
t29 = t126 * t49 - t128 * t184;
t27 = -pkin(4) * t203 - t31;
t23 = t126 * t83 + t128 * t48;
t22 = t126 * t81 + t128 * t46;
t15 = pkin(4) * t202 - t19;
t3 = -pkin(4) * t147 - t7;
t5 = [t208, 0, t156 * t127, (-qJDD(2) * t132 - t135 * t137) * t127, 0, 0, 0, 0, 0, -qJD(3) * t50 - qJDD(3) * t86 + (t156 * t134 - t135 * t178) * t127, -qJD(3) * t49 - qJDD(3) * t87 + (-t156 * t131 - t135 * t179) * t127, t43 * t191 + t50 * t95 + t69 * t86 + (t134 * t29 - t43 * t200) * qJD(2), t44 * t191 + t50 * t97 + t70 * t86 + (t134 * t30 - t44 * t200) * qJD(2), t29 * t97 - t30 * t95 + t43 * t70 - t44 * t69, t17 * t86 - t19 * t29 + t20 * t30 - t43 * t7 + t44 * t8 + t50 * t51 - g(3), 0, 0, 0, 0, 0, -(-qJD(5) * t26 - t130 * t30 + t133 * t50) * t209 + t25 * t63 + t29 * t58 + t43 * t12, (qJD(5) * t25 + t130 * t50 + t133 * t30) * t209 - t26 * t63 + t29 * t60 + t43 * t11; 0, qJDD(2), t208 * t216 + t169, -t208 * t217 + t168, qJDD(2) * t124 + 0.2e1 * t134 * t178, 0.2e1 * t131 * t191 - 0.2e1 * t207 * t195, qJDD(3) * t131 + t134 * t136, qJDD(3) * t134 - t131 * t136, 0, t140 * t131 + t233 * t134, -t233 * t131 + t140 * t134, -g(1) * t34 - g(2) * t33 - g(3) * t71 + (-t95 * t186 + pkin(7) * t69 + t126 * t17 + (qJD(2) * t72 + t19) * qJD(3)) * t131 + (-qJDD(2) * t72 - t7 + (pkin(7) * t95 + t126 * t51) * qJD(3) - t229 * qJD(2)) * t134, (-g(3) * t217 - t168) * t128 + (-t97 * t186 + pkin(7) * t70 + t17 * t128 + (-qJD(2) * t73 - t20) * qJD(3)) * t131 + (t73 * qJDD(2) + t8 + (pkin(7) * t97 + t128 * t51) * qJD(3) + t228 * qJD(2) + t145 * t126) * t134, -t69 * t73 - t70 * t72 - t229 * t97 - t228 * t95 + (-t126 * t20 - t128 * t19) * t199 + (-t126 * t8 - t128 * t7 - t145) * t131, t7 * t72 + t8 * t73 + t228 * t20 + t229 * t19 + t155 * t232 + t155 * t231 + (t17 * t131 + t51 * t199 - t168) * pkin(7) + (-g(3) * pkin(7) * t132 + (-g(3) * t155 - t51 * t205) * t135) * t127, t11 * t90 + t36 * t60, -t11 * t89 - t12 * t90 - t36 * t58 - t37 * t60, -t36 * t209 + t63 * t90 + (t11 * t131 + t60 * t199) * t126, t37 * t209 - t63 * t89 + (-t12 * t131 - t58 * t199) * t126, (t131 * t63 - t199 * t209) * t126, (-t130 * t57 + t133 * t77) * t63 + t56 * t12 + t3 * t89 + t15 * t37 - g(1) * (t133 * t34 - t83 * t212) - g(2) * (t133 * t33 - t81 * t212) - g(3) * (t130 * t189 + t133 * t71) - (t130 * t167 - t133 * t143) * t209 + t230 * t58 + (t2 * t131 - t160 * t199) * t126, -(t130 * t77 + t133 * t57) * t63 + t56 * t11 + t3 * t90 + t15 * t36 - g(1) * (-t130 * t34 - t83 * t211) - g(2) * (-t130 * t33 - t81 * t211) - g(3) * (-t130 * t71 + t133 * t189) - (t130 * t143 + t133 * t167) * t209 + t230 * t60 + (-t1 * t131 - t6 * t199) * t126; 0, 0, 0, 0, -t131 * t137 * t134, t207 * t137, t192, t191, qJDD(3), qJD(3) * t67 - t101 * t203 + t138, -t177 + (-qJD(2) * t101 - t75) * t134 + t150, t126 * t181 - pkin(3) * t69 - t67 * t95 + t146 * t128 + (t126 * t144 - t131 * t19 + t134 * t31) * qJD(2), t128 * t181 - pkin(3) * t70 - t67 * t97 - t146 * t126 + (t128 * t144 + t20 * t131 - t134 * t32) * qJD(2), t31 * t97 + t32 * t95 + (-qJ(4) * t69 - qJD(4) * t95 + t19 * t202 + t8) * t128 + (qJ(4) * t70 + qJD(4) * t97 + t20 * t202 - t7) * t126 - t150, -t19 * t31 - t20 * t32 - t51 * t67 + (-t19 * t126 + t20 * t128) * qJD(4) + t146 * pkin(3) + (-t7 * t126 + t8 * t128 - t150) * qJ(4), -t60 * t79 + (t11 * t133 - t60 * t198) * t126, t58 * t79 + t60 * t78 + (-t11 * t130 - t12 * t133 + (t130 * t58 - t133 * t60) * qJD(5)) * t126, -t11 * t128 + t79 * t209 + (t198 * t209 - t60 * t202 + t224) * t126, t12 * t128 - t78 * t209 + (t197 * t209 + t58 * t202 - t225) * t126, t126 * t202 * t209 - t128 * t63, -t15 * t78 - t27 * t58 + t163 * t133 + t139 * t130 + (t130 * t159 - t133 * t142 - t2) * t128 + (qJ(4) * t12 + qJD(4) * t58 + t3 * t130 + t15 * t197 + t160 * t202) * t126, -t15 * t79 - t27 * t60 - t163 * t130 + t139 * t133 + (t130 * t142 + t133 * t159 + t1) * t128 + (qJ(4) * t11 + qJD(4) * t60 + t3 * t133 - t15 * t198 + t6 * t202) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126 * t192 - t122 + (-t97 + t201) * t202, t128 * t192 + t193 + (t95 + t123) * t202, -t95 ^ 2 - t97 ^ 2, t19 * t97 + t20 * t95 - t138 + t234, 0, 0, 0, 0, 0, -t130 * t174 - t58 * t97 + t224, -t133 * t174 - t60 * t97 - t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, -t209 * t58 + t11, -t209 * t60 - t12, t63, -t6 * t209 - t15 * t60 - g(1) * (-t130 * t23 + t133 * t47) - g(2) * (-t130 * t22 + t133 * t45) - g(3) * t25 + t2, t160 * t209 + t15 * t58 - g(1) * (-t130 * t47 - t133 * t23) - g(2) * (-t130 * t45 - t133 * t22) + g(3) * t26 - t1;];
tau_reg = t5;
