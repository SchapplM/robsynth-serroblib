% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:40
% EndTime: 2019-12-31 20:30:47
% DurationCPUTime: 2.40s
% Computational Cost: add. (2181->324), mult. (4711->412), div. (0->0), fcn. (3259->8), ass. (0->189)
t128 = sin(qJ(4));
t129 = sin(qJ(2));
t132 = cos(qJ(4));
t133 = cos(qJ(2));
t65 = t128 * t129 + t132 * t133;
t245 = t65 * qJD(1);
t252 = qJD(5) + t245;
t119 = qJD(2) - qJD(4);
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t202 = qJD(1) * t133;
t203 = qJD(1) * t129;
t61 = -t128 * t202 + t132 * t203;
t37 = t131 * t119 + t127 * t61;
t258 = t252 * t37;
t160 = t119 * t127 - t131 * t61;
t257 = t160 * t252;
t198 = qJD(4) * t132;
t199 = qJD(4) * t128;
t200 = qJD(2) * t133;
t256 = t128 * t200 + t129 * t198 - t133 * t199;
t108 = pkin(6) * t203;
t255 = -pkin(7) * t203 + qJD(3) + t108;
t112 = t129 * qJDD(1);
t195 = qJD(1) * qJD(2);
t184 = t133 * t195;
t254 = t184 + t112;
t185 = t129 * t195;
t194 = t133 * qJDD(1);
t253 = t185 - t194;
t17 = t256 * qJD(1) + qJDD(1) * t65 - t132 * t185;
t14 = qJDD(5) + t17;
t221 = t127 * t14;
t246 = t131 * t252;
t251 = t160 * t61 + t246 * t252 + t221;
t118 = qJDD(2) - qJDD(4);
t150 = t65 * qJD(4);
t151 = t253 * t128 + t254 * t132;
t16 = -qJD(1) * t150 + t151;
t196 = qJD(5) * t131;
t197 = qJD(5) * t127;
t7 = -t127 * t118 - t119 * t196 + t131 * t16 - t61 * t197;
t225 = t7 * t127;
t250 = t160 * t246 - t225;
t8 = -qJD(5) * t160 + t131 * t118 + t127 * t16;
t249 = (-t8 + t257) * t127 - (-t7 + t258) * t131;
t205 = t133 * pkin(2) + t129 * qJ(3);
t75 = -pkin(1) - t205;
t135 = -pkin(2) - pkin(3);
t188 = t135 * qJD(2);
t48 = t188 + t255;
t122 = qJD(2) * qJ(3);
t109 = pkin(6) * t202;
t72 = -pkin(7) * t202 + t109;
t62 = t122 + t72;
t22 = -t128 * t62 + t132 * t48;
t19 = pkin(4) * t119 - t22;
t248 = t252 * t19;
t247 = t61 * t119 + t17;
t206 = t132 * qJ(3) + t128 * t135;
t134 = cos(qJ(1));
t229 = g(2) * t134;
t130 = sin(qJ(1));
t232 = g(1) * t130;
t244 = -t229 + t232;
t230 = g(2) * t130;
t231 = g(1) * t134;
t167 = t230 + t231;
t209 = t130 * t133;
t212 = t129 * t132;
t53 = t128 * t209 - t130 * t212;
t211 = t129 * t134;
t213 = t128 * t133;
t55 = -t132 * t211 + t134 * t213;
t152 = g(1) * t55 + g(2) * t53 + g(3) * t65;
t103 = pkin(6) * t112;
t93 = pkin(6) * t184;
t186 = qJDD(3) + t103 + t93;
t29 = -t254 * pkin(7) + t135 * qJDD(2) + t186;
t104 = pkin(6) * t194;
t120 = qJDD(2) * qJ(3);
t121 = qJD(2) * qJD(3);
t46 = -pkin(6) * t185 + t104 + t120 + t121;
t30 = t253 * pkin(7) + t46;
t158 = t128 * t30 - t132 * t29 + t62 * t198 + t48 * t199;
t4 = pkin(4) * t118 + t158;
t147 = t152 - t4;
t28 = t61 * pkin(4) + pkin(8) * t245;
t101 = qJ(3) * t202;
t52 = t135 * t203 + t101;
t69 = -pkin(8) + t206;
t241 = (qJD(5) * t69 - t28 + t52) * t252 + t147;
t240 = (pkin(8) * qJD(5) + t28) * t252 - t147;
t219 = t131 * t14;
t239 = t127 * t252 ^ 2 - t37 * t61 - t219;
t216 = pkin(6) * qJDD(2);
t63 = -qJD(1) * pkin(1) - pkin(2) * t202 - qJ(3) * t203;
t238 = (qJD(1) * t75 + t63) * qJD(2) - t216;
t236 = pkin(6) - pkin(7);
t77 = t236 * t129;
t78 = t236 * t133;
t40 = t128 * t78 - t132 * t77;
t201 = qJD(2) * t129;
t71 = t236 * t201;
t73 = qJD(2) * t78;
t10 = -qJD(4) * t40 + t128 * t73 - t132 * t71;
t45 = pkin(3) * t202 - t63;
t12 = pkin(4) * t245 - pkin(8) * t61 + t45;
t154 = t128 * t29 + t132 * t30 + t48 * t198 - t62 * t199;
t182 = -pkin(8) * t118 + qJD(5) * t12 + t154;
t64 = t133 * pkin(3) - t75;
t66 = t212 - t213;
t21 = pkin(4) * t65 - pkin(8) * t66 + t64;
t34 = qJD(2) * t65 - t150;
t41 = t128 * t77 + t132 * t78;
t237 = -t41 * t14 - (qJD(5) * t21 + t10) * t252 - t182 * t65 + t19 * t34 + t4 * t66 + t231;
t54 = t65 * t130;
t235 = g(1) * t54;
t23 = t128 * t48 + t132 * t62;
t20 = -pkin(8) * t119 + t23;
t5 = t12 * t131 - t127 * t20;
t234 = t5 * t61;
t6 = t12 * t127 + t131 * t20;
t233 = t6 * t61;
t228 = t19 * t66;
t227 = t252 * t61;
t226 = t61 * t245;
t157 = -qJ(3) * t128 + t132 * t135;
t224 = qJD(4) * t157 - t128 * t72 + t255 * t132;
t223 = t206 * qJD(4) + t255 * t128 + t132 * t72;
t218 = t245 * t119;
t215 = qJD(5) * t20;
t125 = qJDD(1) * pkin(1);
t214 = qJDD(2) * pkin(2);
t137 = qJD(1) ^ 2;
t210 = t129 * t137;
t113 = t129 * qJD(3);
t207 = qJ(3) * t200 + t113;
t123 = t129 ^ 2;
t124 = t133 ^ 2;
t204 = t123 - t124;
t193 = t245 ^ 2 - t61 ^ 2;
t192 = t252 * t203;
t191 = t66 * t197;
t190 = t252 * t196;
t189 = t133 * t210;
t175 = -qJD(2) * pkin(2) + qJD(3);
t174 = t252 * t119;
t173 = t119 ^ 2;
t172 = t129 * t188;
t171 = g(2) * t54 + g(3) * t66;
t170 = t21 * t14 + t235;
t169 = pkin(2) * t194 + t254 * qJ(3) + qJD(1) * t113 + t125;
t136 = qJD(2) ^ 2;
t168 = pkin(6) * t136 + t229;
t166 = t14 * t66 + t252 * t34;
t163 = pkin(2) * t129 - qJ(3) * t133;
t162 = -t215 - t229;
t74 = t108 + t175;
t76 = t109 + t122;
t161 = t129 * t76 - t133 * t74;
t159 = g(1) * t211 - g(3) * t133 + t129 * t230 - t103;
t155 = -0.2e1 * pkin(1) * t195 - t216;
t153 = -qJDD(3) + t159;
t44 = t172 + t207;
t148 = -t168 + 0.2e1 * t125;
t145 = t171 - t182;
t144 = -pkin(8) * t14 + t22 * t252 + t248;
t142 = -t69 * t14 - t224 * t252 - t248;
t27 = pkin(2) * t185 - t169;
t57 = pkin(2) * t201 - t207;
t141 = -qJD(1) * t57 - qJDD(1) * t75 - t168 - t27;
t18 = pkin(3) * t194 + qJD(1) * t172 + t169;
t140 = t45 * t61 - t152 + t158;
t56 = t65 * t134;
t139 = -g(1) * t56 - t245 * t45 + t154 - t171;
t51 = t186 - t214;
t138 = -qJD(2) * t161 + t51 * t129 + t46 * t133 - t167;
t95 = g(1) * t209;
t68 = pkin(4) - t157;
t67 = pkin(2) * t203 - t101;
t36 = -t127 * t130 + t131 * t56;
t35 = -t127 * t56 - t130 * t131;
t33 = -t132 * t201 + t256;
t11 = qJD(4) * t41 - t128 * t71 - t132 * t73;
t9 = t33 * pkin(4) - t34 * pkin(8) + t44;
t2 = pkin(4) * t17 - pkin(8) * t16 + t18;
t1 = t131 * t2;
t3 = [qJDD(1), t244, t167, qJDD(1) * t123 + 0.2e1 * t129 * t184, 0.2e1 * t129 * t194 - 0.2e1 * t204 * t195, qJDD(2) * t129 + t133 * t136, qJDD(2) * t133 - t129 * t136, 0, t129 * t155 + t133 * t148 + t95, t155 * t133 + (-t148 - t232) * t129, t238 * t129 + t141 * t133 + t95, (t123 + t124) * qJDD(1) * pkin(6) + t138, -t238 * t133 + (t141 + t232) * t129, pkin(6) * t138 + t63 * t57 + (-t244 + t27) * t75, t16 * t66 + t34 * t61, -t16 * t65 - t17 * t66 - t245 * t34 - t33 * t61, -t118 * t66 - t119 * t34, t118 * t65 + t119 * t33, 0, -g(2) * t56 + t11 * t119 + t118 * t40 + t17 * t64 + t18 * t65 + t245 * t44 + t33 * t45 + t235, -g(1) * t53 + g(2) * t55 + t10 * t119 + t118 * t41 + t16 * t64 + t18 * t66 + t34 * t45 + t44 * t61, t160 * t191 + (-t160 * t34 + t66 * t7) * t131, (t127 * t160 - t131 * t37) * t34 + (-t225 - t131 * t8 + (t127 * t37 + t131 * t160) * qJD(5)) * t66, t131 * t166 - t160 * t33 - t191 * t252 + t7 * t65, -t127 * t166 - t190 * t66 - t37 * t33 - t8 * t65, t14 * t65 + t252 * t33, -g(2) * t36 + t1 * t65 + t11 * t37 + t5 * t33 + t40 * t8 + (t9 * t252 + (-t20 * t65 - t252 * t41 + t228) * qJD(5) + t170) * t131 + t237 * t127, -g(2) * t35 - t11 * t160 - t6 * t33 + t40 * t7 + (-(-qJD(5) * t41 + t9) * t252 - (t2 - t215) * t65 - qJD(5) * t228 - t170) * t127 + t237 * t131; 0, 0, 0, -t189, t204 * t137, t112, t194, qJDD(2), pkin(1) * t210 + t159, g(3) * t129 - t104 + (pkin(1) * t137 + t167) * t133, 0.2e1 * t214 + (-t129 * t63 + t133 * t67) * qJD(1) + t153, -t163 * qJDD(1) + ((t76 - t122) * t129 + (t175 - t74) * t133) * qJD(1), t104 + 0.2e1 * t120 + 0.2e1 * t121 + (qJD(1) * t67 - g(3)) * t129 + (qJD(1) * t63 - t167) * t133, t161 * qJD(1) * pkin(6) - t51 * pkin(2) - g(3) * t205 + t46 * qJ(3) + t76 * qJD(3) + t167 * t163 - t63 * t67, -t226, t193, qJD(4) * t245 - t151 + t218, t247, t118, -t157 * t118 + t223 * t119 - t245 * t52 + t140, t206 * t118 + t224 * t119 - t52 * t61 + t139, t250, -t249, -t251, t239, t227, t142 * t127 - t241 * t131 + t223 * t37 + t68 * t8 + t234, t241 * t127 + t142 * t131 - t160 * t223 + t68 * t7 - t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t189, t112, -t123 * t137 - t136, -qJD(2) * t76 + t63 * t203 - t153 - t214 + t93, 0, 0, 0, 0, 0, -t132 * t118 - t128 * t173 - t203 * t245, t128 * t118 - t132 * t173 - t203 * t61, 0, 0, 0, 0, 0, -t131 * t192 + (t127 * t174 - t8) * t132 + (-t119 * t37 - t190 - t221) * t128, t127 * t192 + (t131 * t174 - t7) * t132 + (t119 * t160 + t197 * t252 - t219) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t193, t16 - t218, -t247, -t118, -t119 * t23 - t140, -t119 * t22 - t139, -t250, t249, t251, -t239, -t227, -pkin(4) * t8 + t144 * t127 - t240 * t131 - t23 * t37 - t234, -pkin(4) * t7 + t240 * t127 + t144 * t131 + t160 * t23 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t37, t160 ^ 2 - t37 ^ 2, t7 + t258, -t8 - t257, t14, -g(1) * t35 + t127 * t145 + t131 * t162 + t160 * t19 + t252 * t6 + t1, g(1) * t36 + t19 * t37 + t5 * t252 + (-t162 - t2) * t127 + t145 * t131;];
tau_reg = t3;
