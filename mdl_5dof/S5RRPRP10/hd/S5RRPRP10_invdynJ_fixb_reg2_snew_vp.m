% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:06
% EndTime: 2019-12-31 20:11:11
% DurationCPUTime: 1.44s
% Computational Cost: add. (4194->248), mult. (8799->290), div. (0->0), fcn. (4982->6), ass. (0->165)
t227 = -pkin(2) - pkin(7);
t154 = cos(qJ(2));
t195 = qJD(1) * qJD(2);
t188 = t154 * t195;
t151 = sin(qJ(2));
t193 = t151 * qJDD(1);
t118 = t188 + t193;
t108 = qJDD(4) + t118;
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t199 = qJD(1) * t154;
t112 = t150 * qJD(2) + t153 * t199;
t114 = t153 * qJD(2) - t150 * t199;
t83 = t114 * t112;
t235 = t108 - t83;
t240 = pkin(4) * t235;
t146 = t151 ^ 2;
t157 = qJD(1) ^ 2;
t139 = t146 * t157;
t156 = qJD(2) ^ 2;
t128 = -t139 - t156;
t202 = t151 * t157;
t190 = t154 * t202;
t125 = qJDD(2) - t190;
t201 = t154 * t125;
t239 = pkin(6) * (t151 * t128 + t201);
t238 = t150 * t235;
t237 = t153 * t235;
t137 = t151 * t195;
t192 = t154 * qJDD(1);
t119 = -t137 + t192;
t74 = -t112 * qJD(4) + t153 * qJDD(2) - t150 * t119;
t196 = t151 * qJD(1);
t134 = qJD(4) + t196;
t93 = t134 * t112;
t236 = t74 - t93;
t106 = t112 ^ 2;
t131 = t134 ^ 2;
t75 = -t131 - t106;
t38 = t150 * t75 + t237;
t39 = t153 * t75 - t238;
t234 = pkin(3) * t38 - qJ(3) * t39;
t69 = t108 + t83;
t217 = t150 * t69;
t107 = t114 ^ 2;
t79 = -t107 - t131;
t45 = t153 * t79 - t217;
t215 = t153 * t69;
t46 = -t150 * t79 - t215;
t233 = pkin(3) * t45 - qJ(3) * t46;
t175 = t118 + t188;
t232 = t175 * qJ(3);
t126 = pkin(3) * t196 - qJD(2) * pkin(7);
t133 = pkin(2) * t137;
t189 = qJD(3) * t196;
t136 = -0.2e1 * t189;
t147 = t154 ^ 2;
t152 = sin(qJ(1));
t155 = cos(qJ(1));
t187 = t152 * g(1) - t155 * g(2);
t168 = -qJDD(1) * pkin(1) - t187;
t30 = -t126 * t196 + t133 + t136 + (-pkin(3) * t147 - pkin(6)) * t157 + t227 * t119 - t232 + t168;
t208 = t151 * qJ(3);
t225 = t154 * pkin(2);
t176 = -t208 - t225;
t115 = t176 * qJD(1);
t179 = t155 * g(1) + t152 * g(2);
t212 = qJDD(1) * pkin(6);
t101 = -t157 * pkin(1) - t179 + t212;
t207 = t151 * t101;
t165 = -qJDD(2) * pkin(2) - t156 * qJ(3) + t115 * t196 + qJDD(3) + t207;
t41 = t118 * pkin(3) - qJDD(2) * pkin(7) + (-pkin(3) * t195 - pkin(7) * t202 + g(3)) * t154 + t165;
t13 = t150 * t30 - t153 * t41;
t231 = qJ(5) * t93 + 0.2e1 * qJD(5) * t114 + t13 - t240;
t14 = t150 * t41 + t153 * t30;
t184 = t150 * qJDD(2) + t153 * t119;
t73 = -t114 * qJD(4) - t184;
t86 = t134 * pkin(4) - t114 * qJ(5);
t10 = -t106 * pkin(4) + t73 * qJ(5) - 0.2e1 * qJD(5) * t112 - t134 * t86 + t14;
t209 = t147 * t157;
t230 = t201 + t151 * (-t156 + t209);
t120 = -0.2e1 * t137 + t192;
t130 = -t156 - t209;
t124 = qJDD(2) + t190;
t205 = t151 * t124;
t228 = pkin(6) * (-t154 * t130 + t205) - pkin(1) * t120;
t224 = t154 * g(3);
t163 = (-qJD(4) + t134) * t114 - t184;
t61 = t74 + t93;
t25 = t150 * t163 - t153 * t61;
t27 = t150 * t61 + t153 * t163;
t67 = -t106 - t107;
t223 = pkin(6) * (t151 * t25 + t154 * t67) - pkin(1) * t27;
t56 = (qJD(4) + t134) * t114 + t184;
t222 = pkin(6) * (t151 * t38 + t154 * t56) - pkin(1) * t39;
t221 = pkin(6) * (t151 * t45 + t154 * t236) - pkin(1) * t46;
t194 = qJD(3) * qJD(2);
t140 = 0.2e1 * t194;
t88 = -t151 * g(3) + t154 * t101;
t178 = -t156 * pkin(2) + qJDD(2) * qJ(3) + t115 * t199 + t88;
t162 = t119 * pkin(3) - pkin(7) * t209 + qJD(2) * t126 + t178;
t40 = t140 + t162;
t218 = t150 * t40;
t216 = t153 * t40;
t214 = qJ(5) * t150;
t213 = qJ(5) * t153;
t211 = t134 * t150;
t210 = t134 * t153;
t206 = t151 * t120;
t122 = t139 + t209;
t200 = (t146 + t147) * t212 + pkin(1) * t122;
t191 = t151 * t83;
t186 = pkin(3) * t25 - qJ(3) * t27;
t87 = t207 + t224;
t185 = t151 * t87 + t154 * t88;
t183 = qJ(3) * t67 + t227 * t25;
t182 = qJ(3) * t56 + t227 * t38;
t181 = qJ(3) * t236 + t227 * t45;
t5 = -t153 * t13 + t150 * t14;
t174 = t150 * t13 + t153 * t14;
t173 = t154 * (-t139 + t156) + t205;
t171 = pkin(3) * t56 + t227 * t39;
t170 = pkin(3) * t236 + t227 * t46;
t169 = pkin(3) * t67 + t227 * t27;
t62 = t140 + t178;
t166 = pkin(4) * t79 - t10;
t164 = t154 * t227 - pkin(1) - t208;
t100 = t157 * pkin(6) - t168;
t63 = t165 + t224;
t9 = -t74 * qJ(5) - t231;
t161 = t9 + t240;
t160 = t119 * pkin(2) + t100 - t133;
t159 = -t73 * pkin(4) - t106 * qJ(5) + t114 * t86 + qJDD(5) + t162;
t158 = t160 + 0.2e1 * t189;
t15 = t140 + t159;
t123 = t139 - t209;
t117 = 0.2e1 * t188 + t193;
t90 = -t107 + t131;
t89 = t106 - t131;
t85 = t175 * t151;
t84 = (t119 - t137) * t154;
t80 = t107 - t106;
t77 = t154 * t117 + t206;
t64 = (-t112 * t153 + t114 * t150) * t134;
t55 = pkin(4) * t61;
t52 = -t114 * t211 + t153 * t74;
t51 = t112 * t210 - t150 * t73;
t49 = t151 * t108 + t154 * (t112 * t150 + t114 * t153) * t134;
t48 = t153 * t89 - t217;
t47 = -t150 * t90 + t237;
t32 = t191 + t154 * (-t114 * t210 - t150 * t74);
t31 = -t191 + t154 * (-t112 * t211 - t153 * t73);
t28 = -pkin(4) * t236 - qJ(5) * t69;
t26 = -t150 * t236 - t153 * t56;
t21 = t151 * t61 + t154 * (-t153 * t90 - t238);
t20 = t151 * t163 + t154 * (-t150 * t89 - t215);
t17 = t151 * t80 + t154 * (t150 * t56 - t153 * t236);
t12 = -qJ(5) * t79 + t15;
t11 = -pkin(4) * t56 + qJ(5) * t75 - t159 - 0.2e1 * t194;
t8 = pkin(4) * t9;
t7 = (t61 + t74) * qJ(5) + t231;
t4 = -pkin(4) * t67 + qJ(5) * t163 + t10;
t3 = -pkin(4) * t15 + qJ(5) * t10;
t1 = t150 * t10 + t153 * t9;
t2 = [0, 0, 0, 0, 0, qJDD(1), t187, t179, 0, 0, t85, t77, t173, t84, t230, 0, t154 * t100 - t228, -pkin(1) * t117 - t151 * t100 - t239, t185 + t200, pkin(1) * t100 + pkin(6) * t185, 0, -t173, -t230, t85, t77, t84, t151 * (qJ(3) * t122 + t63) + t154 * (pkin(2) * t122 + t62) + t200, t154 * (-pkin(2) * t120 + t136 - t160) + (-t154 * t175 - t206) * qJ(3) + t228, t151 * t158 + t239 + (pkin(1) + t225) * t117 + (t117 + t175) * t208, pkin(6) * (t151 * t63 + t154 * t62) + (pkin(1) - t176) * (t158 + t232), t32, t17, t21, t31, t20, t49, t151 * (-t13 + t234) + t154 * (t171 + t216) + t222, t151 * (-t14 + t233) + t154 * (t170 - t218) + t221, t151 * t186 + t154 * (t169 - t174) + t223, t164 * t174 + (pkin(3) + pkin(6)) * (t151 * t5 + t154 * t40), t32, t17, t21, t31, t20, t49, t151 * (t161 + t234) + t154 * (-t153 * t11 + t214 * t235 + t171) + t222, t151 * (t166 + t233) + t154 * (-t150 * t12 - t153 * t28 + t170) + t221, t151 * (t186 - t55) + t154 * (-t150 * t7 - t153 * t4 + t169) + t223, t151 * (pkin(3) * t1 + t8) + t154 * (pkin(3) * t15 - t153 * t3 + t9 * t214) + pkin(6) * (t151 * t1 + t154 * t15) + t164 * (t153 * t10 - t150 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t123, t193, t190, t192, qJDD(2), -t87, -t88, 0, 0, qJDD(2), -t193, -t192, -t190, t123, t190, (-pkin(2) * t151 + qJ(3) * t154) * qJDD(1), -pkin(2) * t124 - qJ(3) * t130 + t63, -pkin(2) * t128 + qJ(3) * t125 + t62, -pkin(2) * t63 + qJ(3) * t62, t52, t26, t47, t51, t48, t64, t182 + t218, t181 + t216, t183 - t5, qJ(3) * t40 + t227 * t5, t52, t26, t47, t51, t48, t64, -t150 * t11 - t213 * t235 + t182, t153 * t12 - t150 * t28 + t181, -t150 * t4 + t153 * t7 + t183, qJ(3) * t15 + t227 * t1 - t150 * t3 - t9 * t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t124, t128, t63, 0, 0, 0, 0, 0, 0, t38, t45, t25, t5, 0, 0, 0, 0, 0, 0, t38, t45, t25, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t80, t61, -t83, t163, t108, -t13, -t14, 0, 0, t83, t80, t61, -t83, t163, t108, t161, t166, -t55, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t236, t67, t15;];
tauJ_reg = t2;
