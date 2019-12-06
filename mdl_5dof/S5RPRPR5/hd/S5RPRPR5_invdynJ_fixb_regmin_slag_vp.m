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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:57:16
% EndTime: 2019-12-05 17:57:25
% DurationCPUTime: 2.00s
% Computational Cost: add. (2337->304), mult. (5876->420), div. (0->0), fcn. (4375->12), ass. (0->192)
t141 = sin(pkin(8));
t136 = t141 ^ 2;
t248 = 0.2e1 * t136;
t143 = cos(pkin(8));
t204 = t143 * qJD(1);
t122 = -qJD(3) + t204;
t116 = -qJD(5) + t122;
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t140 = sin(pkin(9));
t142 = cos(pkin(9));
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t103 = t140 * t149 + t142 * t146;
t164 = qJD(1) * t103;
t77 = t141 * t164;
t69 = t148 * t77;
t209 = qJD(1) * t149;
t190 = t141 * t209;
t210 = qJD(1) * t141;
t191 = t146 * t210;
t80 = -t140 * t191 + t142 * t190;
t37 = t145 * t80 + t69;
t233 = t37 * t116;
t205 = qJD(5) * t145;
t224 = t142 * t149;
t172 = t140 * t146 - t224;
t162 = t172 * qJD(3);
t49 = (qJD(1) * t162 - t103 * qJDD(1)) * t141;
t163 = t103 * qJD(3);
t197 = t141 * qJDD(1);
t199 = qJDD(1) * t146;
t50 = -t197 * t224 + (qJD(1) * t163 + t140 * t199) * t141;
t6 = -qJD(5) * t69 + t145 * t49 - t148 * t50 - t80 * t205;
t258 = t6 - t233;
t173 = -t145 * t77 + t148 * t80;
t257 = t173 * t37;
t234 = t173 * t116;
t7 = t173 * qJD(5) - t145 * t50 - t148 * t49;
t256 = -t7 - t234;
t255 = t173 ^ 2 - t37 ^ 2;
t245 = t77 * pkin(7);
t230 = qJ(2) * t149;
t121 = t143 * t230;
t105 = -t143 * pkin(2) - t141 * pkin(6) - pkin(1);
t91 = t105 * qJD(1) + qJD(2);
t59 = -qJ(4) * t191 + qJD(1) * t121 + t146 * t91;
t235 = t142 * t59;
t225 = t141 * t149;
t193 = qJ(4) * t225;
t231 = qJ(2) * t146;
t194 = t143 * t231;
t161 = -t193 - t194;
t87 = t149 * t91;
t58 = t161 * qJD(1) + t87;
t48 = -t122 * pkin(3) + t58;
t26 = t140 * t48 + t235;
t12 = t26 - t245;
t11 = t12 * t205;
t134 = qJ(3) + pkin(9) + qJ(5);
t128 = cos(t134);
t241 = g(1) * t141;
t92 = pkin(3) * t191 + qJ(2) * t210 + qJD(4);
t57 = t77 * pkin(4) + t92;
t127 = sin(t134);
t150 = cos(qJ(1));
t219 = t150 * t127;
t147 = sin(qJ(1));
t222 = t147 * t128;
t73 = t143 * t222 - t219;
t218 = t150 * t128;
t223 = t147 * t127;
t75 = -t143 * t218 - t223;
t254 = -g(2) * t73 - g(3) * t75 + t128 * t241 + t57 * t37 + t11;
t238 = g(3) * t150;
t200 = qJDD(1) * qJ(2);
t202 = qJD(1) * qJD(2);
t169 = t200 + t202;
t252 = t143 * t169;
t216 = t150 * t149;
t221 = t147 * t146;
t95 = t143 * t221 + t216;
t217 = t150 * t146;
t220 = t147 * t149;
t97 = t143 * t217 - t220;
t251 = -g(2) * t95 + g(3) * t97;
t250 = qJD(5) + t116;
t196 = t143 * qJDD(1);
t120 = -qJDD(3) + t196;
t206 = qJD(4) * t141;
t208 = qJD(2) * t143;
t159 = -t146 * t208 - t149 * t206;
t226 = t141 * t146;
t195 = qJ(4) * t226;
t229 = qJD(3) * t91;
t90 = t105 * qJDD(1) + qJDD(2);
t86 = t149 * t90;
t15 = -t146 * t229 - t120 * pkin(3) + t86 + t161 * qJDD(1) + ((-t121 + t195) * qJD(3) + t159) * qJD(1);
t153 = t161 * qJD(3) - t146 * t206;
t189 = t149 * t202;
t207 = qJD(3) * t149;
t184 = qJDD(1) * t121 + t143 * t189 + t146 * t90 + t91 * t207;
t187 = t146 * t197;
t21 = -qJ(4) * t187 + t153 * qJD(1) + t184;
t4 = -t140 * t21 + t142 * t15;
t2 = -t120 * pkin(4) + t50 * pkin(7) + t4;
t5 = t140 * t15 + t142 * t21;
t3 = t49 * pkin(7) + t5;
t192 = -t145 * t3 + t148 * t2;
t72 = t143 * t223 + t218;
t74 = t143 * t219 - t222;
t249 = -g(2) * t72 + g(3) * t74 + t127 * t241 - t57 * t173 + t192;
t137 = t143 ^ 2;
t244 = t80 * pkin(7);
t243 = pkin(3) * t140;
t242 = pkin(3) * t146;
t240 = g(2) * t147;
t239 = qJ(2) * t238;
t232 = t105 * t207 + t149 * t208;
t46 = t153 + t232;
t47 = (-t121 + (qJ(4) * t141 - t105) * t146) * qJD(3) + t159;
t19 = t140 * t47 + t142 * t46;
t53 = t140 * t59;
t29 = t142 * t58 - t53;
t101 = t149 * t105;
t63 = -t193 + t101 + (-pkin(3) - t231) * t143;
t215 = t146 * t105 + t121;
t67 = -t195 + t215;
t32 = t140 * t63 + t142 * t67;
t237 = -t143 * t164 + t163;
t236 = t172 * t204 - t162;
t228 = qJDD(1) * pkin(1);
t151 = qJD(1) ^ 2;
t227 = t136 * t151;
t214 = (pkin(3) * t207 + qJD(2)) * t141;
t213 = pkin(3) * t226 + t141 * qJ(2);
t212 = t136 + t137;
t139 = t149 ^ 2;
t211 = t146 ^ 2 - t139;
t203 = qJD(3) + t122;
t201 = qJD(1) * qJD(3);
t198 = qJDD(1) * t149;
t188 = t146 * t201;
t18 = -t140 * t46 + t142 * t47;
t25 = t142 * t48 - t53;
t28 = -t140 * t58 - t235;
t31 = -t140 * t67 + t142 * t63;
t185 = t212 * t151;
t183 = qJD(1) * t203;
t182 = pkin(3) * t190;
t181 = t120 + t196;
t180 = 0.2e1 * t212;
t179 = qJD(3) * t194;
t178 = -qJD(5) * t172 + t236;
t177 = qJD(5) * t103 + t237;
t176 = g(2) * t150 + g(3) * t147;
t175 = -t238 + t240;
t10 = -t122 * pkin(4) - t244 + t25;
t174 = -t145 * t10 - t148 * t12;
t88 = t103 * t141;
t89 = t172 * t141;
t51 = -t145 * t89 + t148 * t88;
t52 = -t145 * t88 - t148 * t89;
t171 = pkin(3) * t187 + qJ(2) * t197 + qJD(3) * t182 + t141 * t202 + qJDD(4);
t170 = qJD(3) * (t122 + t204);
t129 = t142 * pkin(3) + pkin(4);
t168 = t145 * t129 + t148 * t243;
t167 = t148 * t129 - t145 * t243;
t166 = (t149 * pkin(3) + pkin(2)) * t143 - t141 * (-qJ(4) - pkin(6)) + pkin(1);
t160 = -t176 - t228;
t131 = qJDD(2) - t228;
t158 = -t131 - t160;
t156 = -t122 ^ 2 - t227;
t155 = t180 * t202 + t240;
t111 = -qJDD(5) + t120;
t98 = -t143 * t216 - t221;
t96 = t143 * t220 - t217;
t83 = t141 * t163;
t79 = t141 * t162;
t66 = t80 * pkin(4) + t182;
t64 = t88 * pkin(4) + t213;
t60 = -t79 * pkin(4) + t214;
t30 = -t49 * pkin(4) + t171;
t27 = -t88 * pkin(7) + t32;
t24 = -t143 * pkin(4) + t89 * pkin(7) + t31;
t23 = t52 * qJD(5) - t145 * t83 - t148 * t79;
t22 = -t51 * qJD(5) + t145 * t79 - t148 * t83;
t17 = t29 - t244;
t16 = t28 + t245;
t9 = t79 * pkin(7) + t19;
t8 = t83 * pkin(7) + t18;
t1 = [qJDD(1), t176, -t175, t158 * t143, -t158 * t141, t180 * t200 + t155 - t238, -t239 + (-t131 + t176) * pkin(1) + (t212 * t200 + t155) * qJ(2), (qJDD(1) * t139 - 0.2e1 * t149 * t188) * t136, (-t146 * t198 + t211 * t201) * t248, (t146 * t170 - t181 * t149) * t141, (t181 * t146 + t149 * t170) * t141, t120 * t143, -g(2) * t98 + g(3) * t96 - t101 * t120 - t86 * t143 + (t122 * t143 + (t248 + t137) * qJD(1)) * qJ(2) * t207 + (qJD(3) * t105 * t122 + t169 * t248 + (qJ(2) * t120 + qJD(2) * t122 + t229 + t252) * t143) * t146, (-t179 + t232) * t122 + t215 * t120 + (-qJD(1) * t179 + t184) * t143 - g(2) * t97 - g(3) * t95 + (t189 + (-t188 + t198) * qJ(2)) * t248, t141 * t176 - t18 * t80 - t19 * t77 + t25 * t83 + t26 * t79 + t31 * t50 + t32 * t49 + t4 * t89 - t5 * t88, t5 * t32 + t26 * t19 + t4 * t31 + t25 * t18 + t171 * t213 + t92 * t214 - t239 + (g(2) * t166 - g(3) * t242) * t150 + (-g(2) * (-qJ(2) - t242) + g(3) * t166) * t147, t173 * t22 + t6 * t52, -t173 * t23 - t22 * t37 - t6 * t51 - t52 * t7, -t52 * t111 - t22 * t116 - t6 * t143, t51 * t111 + t23 * t116 + t7 * t143, t111 * t143, -(-t145 * t9 + t148 * t8) * t116 - (-t145 * t27 + t148 * t24) * t111 - t192 * t143 + t60 * t37 + t64 * t7 + t30 * t51 + t57 * t23 - g(2) * t75 + g(3) * t73 + (-(-t145 * t24 - t148 * t27) * t116 - t174 * t143) * qJD(5), -g(2) * t74 - g(3) * t72 - t11 * t143 + t57 * t22 + t30 * t52 + t60 * t173 + t64 * t6 + ((-qJD(5) * t27 + t8) * t116 + t24 * t111 + t2 * t143) * t145 + ((qJD(5) * t24 + t9) * t116 + t27 * t111 + (qJD(5) * t10 + t3) * t143) * t148; 0, 0, 0, -t196, t197, -t185, -qJ(2) * t185 + qJDD(2) + t160, 0, 0, 0, 0, 0, -t149 * t120 + t146 * t156, t146 * t120 + t149 * t156, t103 * t49 - t172 * t50 - t236 * t77 + t237 * t80, t5 * t103 - t172 * t4 - t92 * t210 + t236 * t26 - t237 * t25 - t176, 0, 0, 0, 0, 0, -(-t145 * t103 - t148 * t172) * t111 - t37 * t210 + (t145 * t178 + t148 * t177) * t116, (t148 * t103 - t145 * t172) * t111 - t173 * t210 + (-t145 * t177 + t148 * t178) * t116; 0, 0, 0, 0, 0, 0, 0, t149 * t146 * t227, -t211 * t227, (-t146 * t183 + t198) * t141, (-t149 * t183 - t199) * t141, -t120, t86 + (-t143 * t183 - t227) * t230 + (-t203 * t91 + t241 - t252) * t146 + t251, g(1) * t225 - g(2) * t96 - g(3) * t98 - t87 * t122 + (t203 * t204 + t227) * t231 - t184, (t26 + t28) * t80 - (t25 - t29) * t77 + (t140 * t49 + t142 * t50) * pkin(3), -t25 * t28 - t26 * t29 + (t5 * t140 + t4 * t142 + (g(1) * t146 - t209 * t92) * t141 + t251) * pkin(3), t257, t255, t258, t256, -t111, -t167 * t111 + (-t145 * t17 + t148 * t16) * t116 - t66 * t37 + (t116 * t168 + t174) * qJD(5) + t249, t168 * t111 - t148 * t3 - t145 * t2 - (t145 * t16 + t148 * t17) * t116 - t66 * t173 + (-t148 * t10 + t116 * t167) * qJD(5) + t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 ^ 2 - t80 ^ 2, g(1) * t143 + t141 * t175 + t25 * t80 + t26 * t77 + t171, 0, 0, 0, 0, 0, t7 - t234, t6 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, t255, t258, t256, -t111, t174 * t250 + t249, (t12 * t116 - t2) * t145 + (-t10 * t250 - t3) * t148 + t254;];
tau_reg = t1;
