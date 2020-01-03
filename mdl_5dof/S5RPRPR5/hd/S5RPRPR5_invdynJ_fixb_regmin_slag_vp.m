% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:07
% EndTime: 2020-01-03 11:43:17
% DurationCPUTime: 2.28s
% Computational Cost: add. (2337->303), mult. (5876->425), div. (0->0), fcn. (4375->12), ass. (0->192)
t145 = sin(pkin(8));
t140 = t145 ^ 2;
t252 = 0.2e1 * t140;
t147 = cos(pkin(8));
t206 = t147 * qJD(1);
t122 = -qJD(3) + t206;
t116 = -qJD(5) + t122;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t144 = sin(pkin(9));
t146 = cos(pkin(9));
t150 = sin(qJ(3));
t153 = cos(qJ(3));
t103 = t144 * t153 + t146 * t150;
t167 = qJD(1) * t103;
t77 = t145 * t167;
t69 = t152 * t77;
t211 = qJD(1) * t153;
t192 = t145 * t211;
t212 = qJD(1) * t145;
t193 = t150 * t212;
t80 = -t144 * t193 + t146 * t192;
t37 = t149 * t80 + t69;
t239 = t37 * t116;
t207 = qJD(5) * t149;
t228 = t146 * t153;
t174 = t144 * t150 - t228;
t165 = t174 * qJD(3);
t49 = (qJD(1) * t165 - qJDD(1) * t103) * t145;
t166 = t103 * qJD(3);
t199 = t145 * qJDD(1);
t201 = qJDD(1) * t150;
t50 = -t199 * t228 + (qJD(1) * t166 + t144 * t201) * t145;
t6 = -qJD(5) * t69 + t149 * t49 - t152 * t50 - t80 * t207;
t265 = t6 - t239;
t175 = -t149 * t77 + t152 * t80;
t264 = t175 * t37;
t240 = t175 * t116;
t7 = qJD(5) * t175 - t149 * t50 - t152 * t49;
t263 = -t7 - t240;
t234 = qJDD(1) * pkin(1);
t131 = qJDD(2) - t234;
t151 = sin(qJ(1));
t154 = cos(qJ(1));
t257 = -g(2) * t154 - g(3) * t151;
t262 = t257 - t131;
t261 = t175 ^ 2 - t37 ^ 2;
t249 = t77 * pkin(7);
t236 = qJ(2) * t153;
t121 = t147 * t236;
t105 = -t147 * pkin(2) - t145 * pkin(6) - pkin(1);
t91 = qJD(1) * t105 + qJD(2);
t59 = -qJ(4) * t193 + qJD(1) * t121 + t150 * t91;
t241 = t146 * t59;
t229 = t145 * t153;
t195 = qJ(4) * t229;
t237 = qJ(2) * t150;
t196 = t147 * t237;
t164 = -t195 - t196;
t87 = t153 * t91;
t58 = qJD(1) * t164 + t87;
t48 = -t122 * pkin(3) + t58;
t26 = t144 * t48 + t241;
t12 = t26 - t249;
t11 = t12 * t207;
t134 = qJ(3) + pkin(9) + qJ(5);
t128 = cos(t134);
t246 = g(1) * t145;
t92 = pkin(3) * t193 + qJ(2) * t212 + qJD(4);
t57 = t77 * pkin(4) + t92;
t127 = sin(t134);
t223 = t154 * t127;
t226 = t151 * t128;
t73 = t147 * t226 - t223;
t222 = t154 * t128;
t227 = t151 * t127;
t75 = t147 * t222 + t227;
t260 = g(2) * t73 - g(3) * t75 + t128 * t246 + t57 * t37 + t11;
t202 = qJDD(1) * qJ(2);
t204 = qJD(1) * qJD(2);
t171 = t202 + t204;
t258 = t171 * t147;
t256 = -t145 * (-qJ(4) - pkin(6)) + (t153 * pkin(3) + pkin(2)) * t147;
t220 = t154 * t153;
t225 = t151 * t150;
t95 = -t147 * t225 - t220;
t221 = t154 * t150;
t224 = t151 * t153;
t97 = t147 * t221 - t224;
t255 = -g(2) * t95 - g(3) * t97;
t254 = qJD(5) + t116;
t198 = t147 * qJDD(1);
t120 = -qJDD(3) + t198;
t208 = qJD(4) * t145;
t210 = qJD(2) * t147;
t163 = -t150 * t210 - t153 * t208;
t230 = t145 * t150;
t197 = qJ(4) * t230;
t235 = qJD(3) * t91;
t90 = qJDD(1) * t105 + qJDD(2);
t86 = t153 * t90;
t15 = -t150 * t235 - t120 * pkin(3) + t86 + t164 * qJDD(1) + ((-t121 + t197) * qJD(3) + t163) * qJD(1);
t157 = qJD(3) * t164 - t150 * t208;
t191 = t153 * t204;
t209 = qJD(3) * t153;
t186 = qJDD(1) * t121 + t147 * t191 + t150 * t90 + t91 * t209;
t189 = t150 * t199;
t21 = -qJ(4) * t189 + t157 * qJD(1) + t186;
t4 = -t144 * t21 + t146 * t15;
t2 = -t120 * pkin(4) + t50 * pkin(7) + t4;
t5 = t144 * t15 + t146 * t21;
t3 = t49 * pkin(7) + t5;
t194 = -t149 * t3 + t152 * t2;
t72 = -t147 * t227 - t222;
t74 = t147 * t223 - t226;
t253 = -g(2) * t72 - g(3) * t74 + t127 * t246 - t57 * t175 + t194;
t141 = t147 ^ 2;
t248 = t80 * pkin(7);
t247 = pkin(3) * t144;
t245 = g(2) * t151;
t244 = g(3) * t154;
t238 = t105 * t209 + t153 * t210;
t46 = t157 + t238;
t47 = (-t121 + (qJ(4) * t145 - t105) * t150) * qJD(3) + t163;
t19 = t144 * t47 + t146 * t46;
t53 = t144 * t59;
t29 = t146 * t58 - t53;
t101 = t153 * t105;
t63 = -t195 + t101 + (-pkin(3) - t237) * t147;
t219 = t150 * t105 + t121;
t67 = -t197 + t219;
t32 = t144 * t63 + t146 * t67;
t243 = -t147 * t167 + t166;
t242 = t174 * t206 - t165;
t155 = qJD(1) ^ 2;
t232 = t140 * t155;
t218 = (pkin(3) * t209 + qJD(2)) * t145;
t217 = pkin(3) * t230 + t145 * qJ(2);
t216 = t154 * pkin(1) + t151 * qJ(2);
t214 = t140 + t141;
t143 = t153 ^ 2;
t213 = t150 ^ 2 - t143;
t205 = qJD(3) + t122;
t203 = qJD(1) * qJD(3);
t200 = qJDD(1) * t153;
t190 = t150 * t203;
t18 = -t144 * t46 + t146 * t47;
t25 = t146 * t48 - t53;
t28 = -t144 * t58 - t241;
t31 = -t144 * t67 + t146 * t63;
t187 = t214 * t155;
t185 = qJD(1) * t205;
t184 = pkin(3) * t192;
t183 = t120 + t198;
t182 = 0.2e1 * t214;
t181 = qJD(3) * t196;
t180 = -qJD(5) * t174 + t242;
t179 = qJD(5) * t103 + t243;
t177 = -t244 + t245;
t10 = -t122 * pkin(4) - t248 + t25;
t176 = -t149 * t10 - t152 * t12;
t88 = t103 * t145;
t89 = t174 * t145;
t51 = -t149 * t89 + t152 * t88;
t52 = -t149 * t88 - t152 * t89;
t173 = pkin(3) * t189 + qJ(2) * t199 + qJD(3) * t184 + t145 * t204 + qJDD(4);
t172 = qJD(3) * (t122 + t206);
t129 = t146 * pkin(3) + pkin(4);
t170 = t149 * t129 + t152 * t247;
t169 = t152 * t129 - t149 * t247;
t162 = t234 + t262;
t160 = -t122 ^ 2 - t232;
t159 = t182 * t204 + t244;
t136 = t151 * pkin(1);
t111 = -qJDD(5) + t120;
t98 = t147 * t220 + t225;
t96 = t147 * t224 - t221;
t83 = t145 * t166;
t79 = t145 * t165;
t66 = t80 * pkin(4) + t184;
t64 = t88 * pkin(4) + t217;
t60 = -t79 * pkin(4) + t218;
t30 = -t49 * pkin(4) + t173;
t27 = -t88 * pkin(7) + t32;
t24 = -t147 * pkin(4) + t89 * pkin(7) + t31;
t23 = qJD(5) * t52 - t149 * t83 - t152 * t79;
t22 = -qJD(5) * t51 + t149 * t79 - t152 * t83;
t17 = t29 - t248;
t16 = t28 + t249;
t9 = t79 * pkin(7) + t19;
t8 = t83 * pkin(7) + t18;
t1 = [qJDD(1), t257, t177, t162 * t147, -t162 * t145, t182 * t202 + t159 - t245, -t131 * pkin(1) - g(2) * t216 - g(3) * t136 + (t214 * t202 + t159) * qJ(2), (qJDD(1) * t143 - 0.2e1 * t153 * t190) * t140, (-t150 * t200 + t213 * t203) * t252, (t150 * t172 - t153 * t183) * t145, (t150 * t183 + t153 * t172) * t145, t120 * t147, -g(2) * t98 - g(3) * t96 - t101 * t120 - t86 * t147 + (t122 * t147 + (t252 + t141) * qJD(1)) * qJ(2) * t209 + (qJD(3) * t105 * t122 + t171 * t252 + (qJ(2) * t120 + qJD(2) * t122 + t235 + t258) * t147) * t150, (-t181 + t238) * t122 + t219 * t120 + (-qJD(1) * t181 + t186) * t147 + g(2) * t97 - g(3) * t95 + (t191 + (-t190 + t200) * qJ(2)) * t252, t145 * t257 - t18 * t80 - t19 * t77 + t25 * t83 + t26 * t79 + t31 * t50 + t32 * t49 + t4 * t89 - t5 * t88, t5 * t32 + t26 * t19 + t4 * t31 + t25 * t18 + t173 * t217 + t92 * t218 - g(2) * (pkin(3) * t225 + t216) - g(3) * (t256 * t151 + t136) + (-g(2) * t256 - g(3) * (-pkin(3) * t150 - qJ(2))) * t154, t175 * t22 + t6 * t52, -t175 * t23 - t22 * t37 - t6 * t51 - t52 * t7, -t52 * t111 - t22 * t116 - t6 * t147, t51 * t111 + t23 * t116 + t7 * t147, t111 * t147, -(-t149 * t9 + t152 * t8) * t116 - (-t149 * t27 + t152 * t24) * t111 - t194 * t147 + t60 * t37 + t64 * t7 + t30 * t51 + t57 * t23 - g(2) * t75 - g(3) * t73 + (-(-t149 * t24 - t152 * t27) * t116 - t176 * t147) * qJD(5), g(2) * t74 - g(3) * t72 - t11 * t147 + t57 * t22 + t30 * t52 + t60 * t175 + t64 * t6 + ((-qJD(5) * t27 + t8) * t116 + t24 * t111 + t2 * t147) * t149 + ((qJD(5) * t24 + t9) * t116 + t27 * t111 + (qJD(5) * t10 + t3) * t147) * t152; 0, 0, 0, -t198, t199, -t187, -qJ(2) * t187 - t262, 0, 0, 0, 0, 0, -t153 * t120 + t150 * t160, t150 * t120 + t153 * t160, t103 * t49 - t174 * t50 - t242 * t77 + t243 * t80, t5 * t103 - t174 * t4 - t92 * t212 + t242 * t26 - t243 * t25 - t257, 0, 0, 0, 0, 0, -(-t149 * t103 - t152 * t174) * t111 - t37 * t212 + (t149 * t180 + t152 * t179) * t116, (t152 * t103 - t149 * t174) * t111 - t175 * t212 + (-t149 * t179 + t152 * t180) * t116; 0, 0, 0, 0, 0, 0, 0, t153 * t150 * t232, -t213 * t232, (-t150 * t185 + t200) * t145, (-t153 * t185 - t201) * t145, -t120, t86 + (-t147 * t185 - t232) * t236 + (-t205 * t91 + t246 - t258) * t150 + t255, g(1) * t229 + g(2) * t96 - g(3) * t98 - t87 * t122 + (t205 * t206 + t232) * t237 - t186, (t26 + t28) * t80 - (t25 - t29) * t77 + (t144 * t49 + t146 * t50) * pkin(3), -t25 * t28 - t26 * t29 + (t5 * t144 + t4 * t146 + (g(1) * t150 - t92 * t211) * t145 + t255) * pkin(3), t264, t261, t265, t263, -t111, -t169 * t111 + (-t149 * t17 + t152 * t16) * t116 - t66 * t37 + (t116 * t170 + t176) * qJD(5) + t253, t170 * t111 - t152 * t3 - t149 * t2 - (t149 * t16 + t152 * t17) * t116 - t66 * t175 + (-t152 * t10 + t116 * t169) * qJD(5) + t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 ^ 2 - t80 ^ 2, g(1) * t147 - t145 * t177 + t25 * t80 + t26 * t77 + t173, 0, 0, 0, 0, 0, t7 - t240, t6 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t261, t265, t263, -t111, t176 * t254 + t253, (t12 * t116 - t2) * t149 + (-t10 * t254 - t3) * t152 + t260;];
tau_reg = t1;
