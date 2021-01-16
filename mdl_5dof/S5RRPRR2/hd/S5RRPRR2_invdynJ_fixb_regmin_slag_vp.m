% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:35
% EndTime: 2021-01-15 21:23:46
% DurationCPUTime: 2.80s
% Computational Cost: add. (3837->305), mult. (9353->408), div. (0->0), fcn. (7361->16), ass. (0->174)
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t182 = cos(qJ(2));
t233 = t174 * t182;
t222 = qJD(1) * t233;
t178 = sin(qJ(2));
t231 = qJD(1) * t178;
t117 = t173 * t231 - t222;
t129 = t173 * t182 + t174 * t178;
t120 = t129 * qJD(1);
t200 = t177 * t117 - t181 * t120;
t119 = t129 * qJD(2);
t224 = t182 * qJDD(1);
t225 = t178 * qJDD(1);
t206 = t173 * t225 - t174 * t224;
t87 = qJD(1) * t119 + t206;
t226 = qJD(1) * qJD(2);
t221 = t178 * t226;
t196 = t129 * qJDD(1) - t173 * t221;
t220 = t182 * t226;
t88 = t174 * t220 + t196;
t189 = t200 * qJD(4) - t177 * t88 - t181 * t87;
t229 = qJD(4) * t181;
t230 = qJD(4) * t177;
t23 = -t117 * t229 - t120 * t230 - t177 * t87 + t181 * t88;
t77 = -t181 * t117 - t177 * t120;
t29 = t176 * t77 - t180 * t200;
t190 = -qJD(5) * t29 - t176 * t23 + t180 * t189;
t169 = qJD(2) + qJD(4);
t166 = qJD(5) + t169;
t240 = t29 * t166;
t255 = t190 + t240;
t33 = t176 * t200 + t180 * t77;
t263 = t29 * t33;
t239 = t33 * t166;
t227 = qJD(5) * t180;
t228 = qJD(5) * t176;
t4 = t176 * t189 + t180 * t23 + t200 * t228 + t227 * t77;
t260 = t4 - t239;
t261 = t29 ^ 2 - t33 ^ 2;
t170 = qJ(2) + pkin(9);
t167 = qJ(4) + t170;
t159 = qJ(5) + t167;
t152 = sin(t159);
t153 = cos(t159);
t168 = qJDD(2) + qJDD(4);
t245 = t120 * pkin(7);
t175 = -qJ(3) - pkin(6);
t145 = t175 * t182;
t134 = qJD(1) * t145;
t123 = t173 * t134;
t144 = t175 * t178;
t133 = qJD(1) * t144;
t242 = qJD(2) * pkin(2);
t127 = t133 + t242;
t85 = t174 * t127 + t123;
t51 = qJD(2) * pkin(3) - t245 + t85;
t246 = t117 * pkin(7);
t234 = t174 * t134;
t86 = t173 * t127 - t234;
t57 = t86 - t246;
t201 = -t177 * t51 - t181 * t57;
t212 = qJD(2) * t175;
t115 = -t178 * qJD(3) + t182 * t212;
t84 = qJDD(2) * pkin(2) + t115 * qJD(1) + qJDD(1) * t144;
t114 = t182 * qJD(3) + t178 * t212;
t91 = t114 * qJD(1) - qJDD(1) * t145;
t36 = -t173 * t91 + t174 * t84;
t25 = qJDD(2) * pkin(3) - t88 * pkin(7) + t36;
t37 = t173 * t84 + t174 * t91;
t26 = -t87 * pkin(7) + t37;
t192 = t201 * qJD(4) - t177 * t26 + t181 * t25;
t2 = t168 * pkin(4) - t23 * pkin(8) + t192;
t179 = sin(qJ(1));
t183 = cos(qJ(1));
t208 = g(1) * t183 + g(2) * t179;
t252 = (qJD(4) * t51 + t26) * t181 + t177 * t25 - t57 * t230;
t3 = pkin(8) * t189 + t252;
t160 = t182 * pkin(2) + pkin(1);
t139 = -t160 * qJD(1) + qJD(3);
t94 = t117 * pkin(3) + t139;
t42 = -pkin(4) * t77 + t94;
t270 = -g(3) * t153 + t208 * t152 - t176 * t3 + t180 * t2 - t42 * t29;
t266 = t77 * pkin(8);
t14 = -t201 + t266;
t269 = g(3) * t152 + t14 * t228 + t208 * t153 - t42 * t33;
t237 = t77 * t169;
t268 = t23 - t237;
t238 = t200 * t169;
t267 = t189 - t238;
t265 = pkin(8) * t200;
t264 = t200 * t77;
t262 = t200 ^ 2 - t77 ^ 2;
t156 = sin(t167);
t157 = cos(t167);
t259 = g(3) * t156 + t208 * t157 - t94 * t77 - t252;
t258 = (-t14 * t166 - t2) * t176 + t269;
t257 = -g(3) * t157 + t208 * t156 + t200 * t94 + t192;
t216 = -t177 * t57 + t181 * t51;
t13 = t216 + t265;
t11 = t169 * pkin(4) + t13;
t241 = t180 * t14;
t205 = -t176 * t11 - t241;
t256 = t205 * qJD(5) + t270;
t158 = t174 * pkin(2) + pkin(3);
t251 = pkin(2) * t173;
t209 = t181 * t158 - t177 * t251;
t92 = -t173 * t133 + t234;
t58 = t92 + t246;
t93 = t174 * t133 + t123;
t59 = t93 - t245;
t254 = -t209 * qJD(4) + t177 * t58 + t181 * t59;
t113 = t177 * t158 + t181 * t251;
t253 = -t113 * qJD(4) + t177 * t59 - t181 * t58;
t250 = pkin(2) * t178;
t247 = g(3) * t182;
t95 = t174 * t144 + t173 * t145;
t70 = -t129 * pkin(7) + t95;
t128 = t173 * t178 - t233;
t96 = t173 * t144 - t174 * t145;
t71 = -t128 * pkin(7) + t96;
t243 = t177 * t70 + t181 * t71;
t236 = t254 + t265;
t235 = t253 + t266;
t69 = t174 * t114 + t173 * t115;
t171 = t178 ^ 2;
t232 = -t182 ^ 2 + t171;
t162 = t178 * t242;
t98 = t119 * pkin(3) + t162;
t97 = pkin(2) * t231 + t120 * pkin(3);
t219 = -qJD(5) * t11 - t3;
t214 = -t177 * t71 + t181 * t70;
t68 = -t173 * t114 + t174 * t115;
t207 = g(1) * t179 - g(2) * t183;
t90 = -t177 * t128 + t181 * t129;
t17 = -t90 * pkin(8) + t214;
t89 = t181 * t128 + t177 * t129;
t18 = -t89 * pkin(8) + t243;
t204 = t180 * t17 - t176 * t18;
t203 = t176 * t17 + t180 * t18;
t40 = t176 * t90 + t180 * t89;
t41 = -t176 * t89 + t180 * t90;
t100 = t128 * pkin(3) - t160;
t199 = -0.2e1 * pkin(1) * t226 - pkin(6) * qJDD(2);
t122 = t128 * qJD(2);
t47 = t122 * pkin(7) + t68;
t48 = -t119 * pkin(7) + t69;
t198 = t177 * t47 + t181 * t48 + t70 * t229 - t71 * t230;
t111 = pkin(2) * t221 - t160 * qJDD(1) + qJDD(3);
t184 = qJD(2) ^ 2;
t194 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t184 + t207;
t185 = qJD(1) ^ 2;
t193 = pkin(1) * t185 - pkin(6) * qJDD(1) + t208;
t53 = t87 * pkin(3) + t111;
t191 = -t243 * qJD(4) - t177 * t48 + t181 * t47;
t165 = cos(t170);
t164 = sin(t170);
t163 = qJDD(5) + t168;
t112 = pkin(4) + t209;
t52 = t89 * pkin(4) + t100;
t43 = -pkin(4) * t200 + t97;
t39 = t90 * qJD(4) + t181 * t119 - t177 * t122;
t38 = -t89 * qJD(4) - t177 * t119 - t181 * t122;
t27 = t39 * pkin(4) + t98;
t10 = -pkin(4) * t189 + t53;
t9 = t41 * qJD(5) + t176 * t38 + t180 * t39;
t8 = -t40 * qJD(5) - t176 * t39 + t180 * t38;
t7 = -t38 * pkin(8) + t191;
t6 = -t39 * pkin(8) + t198;
t1 = [qJDD(1), t207, t208, t171 * qJDD(1) + 0.2e1 * t178 * t220, 0.2e1 * t178 * t224 - 0.2e1 * t232 * t226, qJDD(2) * t178 + t184 * t182, qJDD(2) * t182 - t184 * t178, 0, t178 * t199 + t182 * t194, -t178 * t194 + t182 * t199, t95 * qJDD(2) + t111 * t128 + t139 * t119 - t160 * t87 + t207 * t165 + (t117 * t250 + t68) * qJD(2), -t96 * qJDD(2) + t111 * t129 - t139 * t122 - t160 * t88 - t207 * t164 + (t120 * t250 - t69) * qJD(2), -t69 * t117 - t86 * t119 - t68 * t120 + t85 * t122 - t37 * t128 - t36 * t129 - t96 * t87 - t95 * t88 - t208, t37 * t96 + t86 * t69 + t36 * t95 + t85 * t68 - t111 * t160 + t139 * t162 - g(1) * (-t179 * t160 - t175 * t183) - g(2) * (t183 * t160 - t179 * t175), -t200 * t38 + t23 * t90, t189 * t90 + t200 * t39 - t23 * t89 + t38 * t77, t90 * t168 + t38 * t169, -t89 * t168 - t39 * t169, 0, -t100 * t189 + t207 * t157 + t214 * t168 + t191 * t169 + t94 * t39 + t53 * t89 - t77 * t98, t100 * t23 - t207 * t156 - t243 * t168 - t198 * t169 - t200 * t98 + t94 * t38 + t53 * t90, t29 * t8 + t4 * t41, t190 * t41 - t29 * t9 + t33 * t8 - t4 * t40, t41 * t163 + t8 * t166, -t40 * t163 - t9 * t166, 0, -t27 * t33 - t52 * t190 + t10 * t40 + t42 * t9 + (-qJD(5) * t203 - t176 * t6 + t180 * t7) * t166 + t204 * t163 + t207 * t153, t27 * t29 + t52 * t4 + t10 * t41 + t42 * t8 - (qJD(5) * t204 + t176 * t7 + t180 * t6) * t166 - t203 * t163 - t207 * t152; 0, 0, 0, -t178 * t185 * t182, t232 * t185, t225, t224, qJDD(2), t178 * t193 - t247, g(3) * t178 + t182 * t193, -g(3) * t165 - t92 * qJD(2) - t139 * t120 + t208 * t164 + (qJDD(2) * t174 - t117 * t231) * pkin(2) + t36, g(3) * t164 + t93 * qJD(2) + t139 * t117 + t208 * t165 + (-qJDD(2) * t173 - t120 * t231) * pkin(2) - t37, (t86 + t92) * t120 + (-t85 + t93) * t117 + (-t173 * t87 - t174 * t88) * pkin(2), -t85 * t92 - t86 * t93 + (-t247 + t173 * t37 + t174 * t36 + (-qJD(1) * t139 + t208) * t178) * pkin(2), t264, t262, t268, t267, t168, t209 * t168 + t253 * t169 + t77 * t97 + t257, -t113 * t168 + t254 * t169 + t200 * t97 + t259, -t263, t261, t260, t255, t163, (t180 * t112 - t176 * t113) * t163 + t43 * t33 + (t176 * t236 + t180 * t235) * t166 + ((-t176 * t112 - t180 * t113) * t166 + t205) * qJD(5) + t270, -t43 * t29 + (-t112 * t163 - t2 + (qJD(5) * t113 - t235) * t166) * t176 + (-t113 * t163 + (-qJD(5) * t112 + t236) * t166 + t219) * t180 + t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t120 * qJD(2) + t206, (-t117 + t222) * qJD(2) + t196, -t117 ^ 2 - t120 ^ 2, t86 * t117 + t85 * t120 + t111 - t207, 0, 0, 0, 0, 0, -t189 - t238, t23 + t237, 0, 0, 0, 0, 0, -t190 + t240, t4 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t262, t268, t267, t168, -t201 * t169 + t257, t216 * t169 + t259, -t263, t261, t260, t255, t163, -(-t176 * t13 - t241) * t166 + (t180 * t163 - t166 * t228 - t200 * t33) * pkin(4) + t256, (t13 * t166 + t219) * t180 + (-t176 * t163 - t166 * t227 + t200 * t29) * pkin(4) + t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, t261, t260, t255, t163, -t166 * t205 + t256, (-t3 + (-qJD(5) + t166) * t11) * t180 + t258;];
tau_reg = t1;
