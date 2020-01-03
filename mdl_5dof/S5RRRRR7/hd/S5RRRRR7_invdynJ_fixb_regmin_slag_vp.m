% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR7
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
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:30
% EndTime: 2019-12-31 22:22:39
% DurationCPUTime: 3.56s
% Computational Cost: add. (5557->348), mult. (13097->466), div. (0->0), fcn. (10007->14), ass. (0->204)
t184 = cos(qJ(5));
t253 = qJD(5) * t184;
t181 = sin(qJ(3));
t182 = sin(qJ(2));
t260 = qJD(1) * t182;
t242 = t181 * t260;
t186 = cos(qJ(3));
t187 = cos(qJ(2));
t259 = qJD(1) * t187;
t243 = t186 * t259;
t110 = -t242 + t243;
t111 = -t181 * t259 - t186 * t260;
t180 = sin(qJ(4));
t185 = cos(qJ(4));
t80 = t185 * t110 + t111 * t180;
t308 = t184 * t80;
t313 = t253 - t308;
t179 = sin(qJ(5));
t174 = qJDD(2) + qJDD(3);
t169 = qJDD(4) + t174;
t175 = qJD(2) + qJD(3);
t170 = qJD(4) + t175;
t215 = t110 * t180 - t185 * t111;
t254 = qJD(5) * t179;
t255 = qJD(4) * t185;
t256 = qJD(4) * t180;
t250 = t187 * qJDD(1);
t252 = qJD(1) * qJD(2);
t240 = t187 * t252;
t251 = t182 * qJDD(1);
t305 = t240 + t251;
t65 = qJD(3) * t243 - t175 * t242 + t181 * t250 + t305 * t186;
t217 = t181 * t251 - t186 * t250;
t121 = t181 * t187 + t182 * t186;
t90 = t175 * t121;
t66 = qJD(1) * t90 + t217;
t30 = t110 * t255 + t111 * t256 - t180 * t66 + t185 * t65;
t18 = t179 * t169 + t170 * t253 + t184 * t30 - t215 * t254;
t69 = t170 * t179 + t184 * t215;
t19 = qJD(5) * t69 - t184 * t169 + t179 * t30;
t307 = qJD(5) - t80;
t312 = t179 * t307;
t67 = -t184 * t170 + t179 * t215;
t1 = -t179 * t19 + t18 * t184 - t312 * t69 - t313 * t67;
t31 = qJD(4) * t215 + t180 * t65 + t185 * t66;
t29 = qJDD(5) + t31;
t7 = t184 * t29 + t215 * t67 - t307 * t312;
t178 = qJ(2) + qJ(3);
t173 = qJ(4) + t178;
t161 = sin(t173);
t183 = sin(qJ(1));
t188 = cos(qJ(1));
t220 = g(1) * t188 + g(2) * t183;
t311 = t220 * t161;
t16 = t18 * t179;
t9 = t313 * t69 + t16;
t73 = t307 * t253;
t8 = t179 * t29 - t215 * t69 - t307 * t308 + t73;
t293 = pkin(6) + pkin(7);
t139 = t293 * t187;
t129 = qJD(1) * t139;
t116 = t186 * t129;
t138 = t293 * t182;
t127 = qJD(1) * t138;
t278 = qJD(2) * pkin(2);
t118 = -t127 + t278;
t214 = -t118 * t181 - t116;
t291 = pkin(8) * t110;
t64 = -t214 + t291;
t277 = t180 * t64;
t106 = t111 * pkin(8);
t112 = t181 * t129;
t224 = t186 * t118 - t112;
t63 = t106 + t224;
t56 = pkin(3) * t175 + t63;
t36 = t185 * t56 - t277;
t34 = -pkin(4) * t170 - t36;
t287 = t34 * t80;
t162 = cos(t173);
t288 = g(3) * t162;
t91 = qJDD(2) * pkin(2) - t293 * t305;
t241 = t182 * t252;
t92 = t293 * (-t241 + t250);
t200 = qJD(3) * t214 - t181 * t92 + t186 * t91;
t23 = pkin(3) * t174 - pkin(8) * t65 + t200;
t258 = qJD(3) * t181;
t295 = (qJD(3) * t118 + t92) * t186 - t129 * t258 + t181 * t91;
t28 = -pkin(8) * t66 + t295;
t274 = t185 * t64;
t37 = t180 * t56 + t274;
t297 = qJD(4) * t37 + t180 * t28 - t185 * t23;
t4 = -pkin(4) * t169 + t297;
t309 = t4 + t288;
t283 = t215 * t80;
t24 = t215 ^ 2 - t80 ^ 2;
t20 = -t170 * t80 + t30;
t156 = g(3) * t161;
t296 = (qJD(4) * t56 + t28) * t185 + t180 * t23 - t64 * t256;
t166 = -pkin(2) * t187 - pkin(1);
t137 = t166 * qJD(1);
t93 = -pkin(3) * t110 + t137;
t195 = t162 * t220 - t93 * t80 + t156 - t296;
t303 = pkin(4) * t215;
t163 = pkin(3) * t180 + pkin(9);
t292 = pkin(3) * t111;
t46 = -pkin(9) * t80 - t292 + t303;
t302 = (qJD(5) * t163 + t46) * t307;
t165 = pkin(2) * t186 + pkin(3);
t268 = t181 * t185;
t262 = pkin(2) * t268 + t180 * t165;
t105 = pkin(9) + t262;
t167 = pkin(2) * t260;
t301 = (qJD(5) * t105 + t167 + t46) * t307;
t300 = (t307 * pkin(9) + t303) * t307;
t284 = t307 * t215;
t263 = -t181 * t138 + t186 * t139;
t35 = pkin(9) * t170 + t37;
t40 = -pkin(4) * t80 - pkin(9) * t215 + t93;
t13 = -t179 * t35 + t184 * t40;
t213 = -t13 * t215 + t184 * t311 + t34 * t254;
t14 = t179 * t40 + t184 * t35;
t216 = t14 * t215 + t309 * t179 + t34 * t253;
t192 = -t215 * t93 - t288 - t297 + t311;
t21 = t170 * t215 - t31;
t244 = qJD(2) * t293;
t128 = t182 * t244;
t130 = t187 * t244;
t257 = qJD(3) * t186;
t204 = -t186 * t128 - t181 * t130 - t138 * t257 - t139 * t258;
t47 = -pkin(8) * t90 + t204;
t198 = -t263 * qJD(3) + t128 * t181 - t186 * t130;
t120 = t181 * t182 - t186 * t187;
t89 = t175 * t120;
t48 = pkin(8) * t89 + t198;
t222 = -t186 * t138 - t139 * t181;
t74 = -pkin(8) * t121 + t222;
t75 = -pkin(8) * t120 + t263;
t50 = t180 * t75 - t185 * t74;
t10 = -qJD(4) * t50 + t180 * t48 + t185 * t47;
t238 = pkin(9) * t169 + qJD(5) * t40 + t296;
t87 = t185 * t120 + t121 * t180;
t44 = -qJD(4) * t87 - t180 * t90 - t185 * t89;
t88 = -t120 * t180 + t121 * t185;
t96 = pkin(3) * t120 + t166;
t49 = pkin(4) * t87 - pkin(9) * t88 + t96;
t51 = t180 * t74 + t185 * t75;
t294 = -(qJD(5) * t49 + t10) * t307 - t238 * t87 - t51 * t29 + t34 * t44 + t4 * t88;
t286 = t34 * t88;
t285 = t49 * t29;
t269 = t180 * t181;
t223 = t127 * t181 - t116;
t71 = t223 - t291;
t264 = -t186 * t127 - t112;
t72 = t106 + t264;
t280 = t180 * t71 + t185 * t72 - t165 * t255 - (-t181 * t256 + (t185 * t186 - t269) * qJD(3)) * pkin(2);
t279 = -t180 * t72 + t185 * t71 + t165 * t256 + (t181 * t255 + (t180 * t186 + t268) * qJD(3)) * pkin(2);
t272 = t111 * t110;
t271 = t179 * t183;
t270 = t179 * t188;
t267 = t183 * t184;
t266 = t184 * t188;
t176 = t182 ^ 2;
t261 = -t187 ^ 2 + t176;
t168 = t182 * t278;
t246 = t88 * t254;
t82 = pkin(3) * t90 + t168;
t107 = pkin(2) * t241 + t166 * qJDD(1);
t55 = pkin(3) * t66 + t107;
t6 = pkin(4) * t31 - pkin(9) * t30 + t55;
t237 = qJD(5) * t35 - t6;
t38 = t180 * t63 + t274;
t221 = pkin(3) * t256 - t38;
t219 = g(1) * t183 - g(2) * t188;
t218 = t29 * t88 + t307 * t44;
t212 = -t238 + t156;
t211 = -pkin(2) * t269 + t165 * t185;
t209 = -0.2e1 * pkin(1) * t252 - pkin(6) * qJDD(2);
t206 = -pkin(9) * t29 + t307 * t36 - t287;
t203 = -t105 * t29 + t280 * t307 - t287;
t189 = qJD(2) ^ 2;
t202 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t189 + t219;
t190 = qJD(1) ^ 2;
t201 = pkin(1) * t190 - pkin(6) * qJDD(1) + t220;
t39 = t185 * t63 - t277;
t197 = -t163 * t29 - t287 + (-pkin(3) * t255 + t39) * t307;
t171 = sin(t178);
t172 = cos(t178);
t193 = g(3) * t171 - t137 * t110 + t172 * t220 - t295;
t191 = -g(3) * t172 + t137 * t111 + t171 * t220 + t200;
t164 = -pkin(3) * t185 - pkin(4);
t104 = -pkin(4) - t211;
t102 = t162 * t266 + t271;
t101 = -t162 * t270 + t267;
t100 = -t162 * t267 + t270;
t99 = t162 * t271 + t266;
t94 = t167 - t292;
t70 = -t110 ^ 2 + t111 ^ 2;
t54 = -t217 + (-qJD(1) * t121 - t111) * t175;
t53 = -t110 * t175 + t65;
t45 = qJD(4) * t88 - t180 * t89 + t185 * t90;
t12 = pkin(4) * t45 - pkin(9) * t44 + t82;
t11 = qJD(4) * t51 + t180 * t47 - t185 * t48;
t5 = t184 * t6;
t2 = [qJDD(1), t219, t220, qJDD(1) * t176 + 0.2e1 * t182 * t240, 0.2e1 * t182 * t250 - 0.2e1 * t261 * t252, qJDD(2) * t182 + t187 * t189, qJDD(2) * t187 - t182 * t189, 0, t182 * t209 + t187 * t202, -t182 * t202 + t187 * t209, t111 * t89 + t121 * t65, -t110 * t89 + t111 * t90 - t120 * t65 - t121 * t66, t121 * t174 - t175 * t89, -t120 * t174 - t175 * t90, 0, t107 * t120 - t110 * t168 + t137 * t90 + t166 * t66 + t172 * t219 + t174 * t222 + t175 * t198, t107 * t121 - t111 * t168 - t137 * t89 + t166 * t65 - t219 * t171 - t263 * t174 - t204 * t175, t215 * t44 + t30 * t88, -t215 * t45 - t30 * t87 - t31 * t88 + t44 * t80, t169 * t88 + t170 * t44, -t169 * t87 - t170 * t45, 0, -t11 * t170 + t162 * t219 - t169 * t50 + t31 * t96 + t45 * t93 + t55 * t87 - t80 * t82, -t10 * t170 - t161 * t219 - t169 * t51 + t215 * t82 + t30 * t96 + t44 * t93 + t55 * t88, -t69 * t246 + (t18 * t88 + t44 * t69) * t184, (-t179 * t69 - t184 * t67) * t44 + (-t16 - t184 * t19 + (t179 * t67 - t184 * t69) * qJD(5)) * t88, t18 * t87 + t184 * t218 - t246 * t307 + t45 * t69, -t179 * t218 - t19 * t87 - t45 * t67 - t73 * t88, t29 * t87 + t307 * t45, -g(1) * t100 - g(2) * t102 + t11 * t67 + t13 * t45 + t50 * t19 + t5 * t87 + (t12 * t307 + t285 + (-t307 * t51 - t35 * t87 + t286) * qJD(5)) * t184 + t294 * t179, -g(1) * t99 - g(2) * t101 + t11 * t69 - t14 * t45 + t50 * t18 + (-(-qJD(5) * t51 + t12) * t307 - t285 + t237 * t87 - qJD(5) * t286) * t179 + t294 * t184; 0, 0, 0, -t182 * t190 * t187, t261 * t190, t251, t250, qJDD(2), -g(3) * t187 + t182 * t201, g(3) * t182 + t187 * t201, t272, t70, t53, t54, t174, -t223 * t175 + (t110 * t260 + t174 * t186 - t175 * t258) * pkin(2) + t191, t264 * t175 + (t111 * t260 - t174 * t181 - t175 * t257) * pkin(2) + t193, -t283, t24, t20, t21, t169, t211 * t169 - t279 * t170 + t80 * t94 + t192, -t262 * t169 + t280 * t170 - t215 * t94 + t195, t9, t1, t8, t7, -t284, t104 * t19 + t279 * t67 + (-t309 - t301) * t184 + t203 * t179 + t213, t104 * t18 + t279 * t69 + t203 * t184 + (-t311 + t301) * t179 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t70, t53, t54, t174, -t175 * t214 + t191, t175 * t224 + t193, -t283, t24, t20, t21, t169, t170 * t38 + (-t111 * t80 + t169 * t185 - t170 * t256) * pkin(3) + t192, t170 * t39 + (t111 * t215 - t169 * t180 - t170 * t255) * pkin(3) + t195, t9, t1, t8, t7, -t284, t164 * t19 + t221 * t67 + (-t309 - t302) * t184 + t197 * t179 + t213, t164 * t18 + t221 * t69 + t197 * t184 + (-t311 + t302) * t179 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283, t24, t20, t21, t169, t170 * t37 + t192, t170 * t36 + t195, t9, t1, t8, t7, -t284, -pkin(4) * t19 - t37 * t67 + t206 * t179 + (-t309 - t300) * t184 + t213, -pkin(4) * t18 - t37 * t69 + t206 * t184 + (-t311 + t300) * t179 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67, -t67 ^ 2 + t69 ^ 2, t307 * t67 + t18, t307 * t69 - t19, t29, -g(1) * t101 + g(2) * t99 + t14 * t307 + t179 * t212 - t253 * t35 - t34 * t69 + t5, g(1) * t102 - g(2) * t100 + t13 * t307 + t179 * t237 + t184 * t212 + t34 * t67;];
tau_reg = t2;
