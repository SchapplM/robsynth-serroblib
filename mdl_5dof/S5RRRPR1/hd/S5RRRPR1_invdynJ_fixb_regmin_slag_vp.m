% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:54
% EndTime: 2019-12-05 18:39:01
% DurationCPUTime: 2.61s
% Computational Cost: add. (4168->305), mult. (10255->409), div. (0->0), fcn. (7724->14), ass. (0->185)
t177 = qJD(2) + qJD(3);
t171 = qJD(5) + t177;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t188 = cos(qJ(3));
t189 = cos(qJ(2));
t240 = qJD(1) * t189;
t230 = t188 * t240;
t184 = sin(qJ(3));
t185 = sin(qJ(2));
t241 = qJD(1) * t185;
t231 = t184 * t241;
t115 = -t230 + t231;
t117 = -t184 * t240 - t188 * t241;
t181 = sin(pkin(9));
t182 = cos(pkin(9));
t211 = -t115 * t182 + t117 * t181;
t254 = t187 * t211;
t93 = t115 * t181 + t117 * t182;
t54 = t183 * t93 + t254;
t256 = t171 * t54;
t237 = qJD(5) * t183;
t234 = t189 * qJDD(1);
t236 = qJD(1) * qJD(2);
t228 = t189 * t236;
t235 = t185 * qJDD(1);
t273 = t228 + t235;
t79 = qJD(3) * t230 - t177 * t231 + t184 * t234 + t273 * t188;
t127 = t184 * t189 + t185 * t188;
t102 = t177 * t127;
t215 = t184 * t235 - t188 * t234;
t80 = qJD(1) * t102 + t215;
t39 = -t181 * t79 - t182 * t80;
t40 = -t181 * t80 + t182 * t79;
t9 = qJD(5) * t254 + t183 * t39 + t187 * t40 + t237 * t93;
t4 = t9 - t256;
t272 = t183 * t211 - t187 * t93;
t261 = t272 * t54;
t8 = t272 ^ 2 - t54 ^ 2;
t180 = qJ(2) + qJ(3);
t164 = pkin(9) + qJ(5) + t180;
t155 = sin(t164);
t156 = cos(t164);
t186 = sin(qJ(1));
t190 = cos(qJ(1));
t218 = g(1) * t190 + g(2) * t186;
t267 = pkin(8) * t211;
t268 = pkin(6) + pkin(7);
t144 = t268 * t189;
t134 = qJD(1) * t144;
t122 = t188 * t134;
t143 = t268 * t185;
t132 = qJD(1) * t143;
t258 = qJD(2) * pkin(2);
t124 = -t132 + t258;
t210 = -t124 * t184 - t122;
t249 = qJ(4) * t115;
t78 = -t210 - t249;
t255 = t182 * t78;
t112 = t117 * qJ(4);
t118 = t184 * t134;
t224 = t188 * t124 - t118;
t77 = t112 + t224;
t68 = pkin(3) * t177 + t77;
t33 = t181 * t68 + t255;
t25 = t33 + t267;
t174 = t189 * pkin(2);
t260 = pkin(1) + t174;
t142 = t260 * qJD(1);
t104 = pkin(3) * t115 + qJD(4) - t142;
t63 = -pkin(4) * t211 + t104;
t206 = g(3) * t155 + t218 * t156 + t25 * t237 - t63 * t54;
t10 = qJD(5) * t272 + t183 * t40 - t187 * t39;
t257 = t171 * t272;
t5 = -t10 + t257;
t172 = sin(t180);
t173 = cos(t180);
t275 = -g(3) * t173 + t172 * t218;
t175 = qJDD(2) + qJDD(3);
t103 = qJDD(2) * pkin(2) - t268 * t273;
t229 = t185 * t236;
t105 = t268 * (-t229 + t234);
t196 = qJD(3) * t210 + t188 * t103 - t184 * t105;
t20 = pkin(3) * t175 - qJ(4) * t79 + qJD(4) * t117 + t196;
t239 = qJD(3) * t184;
t269 = (qJD(3) * t124 + t105) * t188 + t184 * t103 - t134 * t239;
t22 = -qJ(4) * t80 - qJD(4) * t115 + t269;
t6 = -t181 * t22 + t182 * t20;
t2 = pkin(4) * t175 - pkin(8) * t40 + t6;
t7 = t181 * t20 + t182 * t22;
t3 = pkin(8) * t39 + t7;
t199 = -g(3) * t156 + t218 * t155 - t183 * t3 + t187 * t2 - t272 * t63;
t90 = t93 * pkin(8);
t246 = t182 * t184;
t259 = pkin(2) * qJD(3);
t223 = t132 * t184 - t122;
t81 = t223 + t249;
t245 = -t188 * t132 - t118;
t82 = t112 + t245;
t252 = t181 * t82 - t182 * t81 + (-t181 * t188 - t246) * t259;
t247 = t181 * t184;
t250 = -t181 * t81 - t182 * t82 + (t182 * t188 - t247) * t259;
t244 = -t184 * t143 + t188 * t144;
t270 = qJD(5) - t171;
t266 = pkin(3) * t117;
t265 = pkin(3) * t181;
t126 = t184 * t185 - t188 * t189;
t232 = qJD(2) * t268;
t133 = t185 * t232;
t135 = t189 * t232;
t238 = qJD(3) * t188;
t203 = -t188 * t133 - t184 * t135 - t143 * t238 - t144 * t239;
t45 = -qJ(4) * t102 - qJD(4) * t126 + t203;
t101 = t177 * t126;
t195 = -qJD(3) * t244 + t133 * t184 - t188 * t135;
t46 = qJ(4) * t101 - qJD(4) * t127 + t195;
t16 = t181 * t46 + t182 * t45;
t69 = t181 * t78;
t38 = t182 * t77 - t69;
t222 = -t188 * t143 - t184 * t144;
t91 = -qJ(4) * t127 + t222;
t92 = -qJ(4) * t126 + t244;
t50 = t181 * t91 + t182 * t92;
t253 = t267 + t252;
t251 = t90 - t250;
t248 = t117 * t115;
t243 = pkin(3) * t173 + t174;
t178 = t185 ^ 2;
t242 = -t189 ^ 2 + t178;
t169 = t185 * t258;
t227 = pkin(3) * t102 + t169;
t15 = -t181 * t45 + t182 * t46;
t32 = t182 * t68 - t69;
t37 = -t181 * t77 - t255;
t49 = -t181 * t92 + t182 * t91;
t219 = pkin(3) * t126 - t260;
t166 = pkin(2) * t188 + pkin(3);
t110 = -pkin(2) * t247 + t182 * t166;
t67 = -pkin(4) * t93 - t266;
t217 = g(1) * t186 - g(2) * t190;
t216 = t211 * t32 - t33 * t93;
t23 = pkin(4) * t177 + t32 + t90;
t214 = -t183 * t23 - t187 * t25;
t99 = -t126 * t181 + t127 * t182;
t30 = -pkin(8) * t99 + t49;
t98 = -t126 * t182 - t127 * t181;
t31 = pkin(8) * t98 + t50;
t213 = -t183 * t31 + t187 * t30;
t212 = t183 * t30 + t187 * t31;
t56 = t183 * t99 - t187 * t98;
t57 = t183 * t98 + t187 * t99;
t161 = pkin(3) * t182 + pkin(4);
t209 = t161 * t183 + t187 * t265;
t208 = t161 * t187 - t183 * t265;
t207 = -0.2e1 * pkin(1) * t236 - pkin(6) * qJDD(2);
t113 = pkin(2) * t229 - qJDD(1) * t260;
t191 = qJD(2) ^ 2;
t201 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t191 + t217;
t192 = qJD(1) ^ 2;
t200 = pkin(1) * t192 - pkin(6) * qJDD(1) + t218;
t198 = t80 * pkin(3) + qJDD(4) + t113;
t194 = g(3) * t172 - t142 * t115 + t218 * t173 - t269;
t193 = -t142 * t117 + t196 + t275;
t176 = -qJ(4) - t268;
t170 = qJDD(5) + t175;
t168 = pkin(2) * t241;
t131 = pkin(1) + t243;
t111 = pkin(2) * t246 + t166 * t181;
t106 = pkin(4) + t110;
t83 = -t115 ^ 2 + t117 ^ 2;
t73 = -pkin(4) * t98 + t219;
t64 = t168 + t67;
t62 = -t101 * t182 - t102 * t181;
t61 = t101 * t181 - t102 * t182;
t60 = -t215 + (-qJD(1) * t127 - t117) * t177;
t59 = t115 * t177 + t79;
t44 = -pkin(4) * t61 + t227;
t27 = t90 + t38;
t26 = t37 - t267;
t17 = -pkin(4) * t39 + t198;
t14 = qJD(5) * t57 + t183 * t62 - t187 * t61;
t13 = -qJD(5) * t56 + t183 * t61 + t187 * t62;
t12 = pkin(8) * t61 + t16;
t11 = -pkin(8) * t62 + t15;
t1 = [qJDD(1), t217, t218, qJDD(1) * t178 + 0.2e1 * t185 * t228, 0.2e1 * t185 * t234 - 0.2e1 * t242 * t236, qJDD(2) * t185 + t189 * t191, qJDD(2) * t189 - t185 * t191, 0, t185 * t207 + t189 * t201, -t185 * t201 + t189 * t207, t101 * t117 + t127 * t79, t101 * t115 + t102 * t117 - t126 * t79 - t127 * t80, -t101 * t177 + t127 * t175, -t102 * t177 - t126 * t175, 0, -t142 * t102 + t113 * t126 + t115 * t169 + t173 * t217 + t175 * t222 + t177 * t195 - t260 * t80, t142 * t101 + t113 * t127 - t117 * t169 - t217 * t172 - t244 * t175 - t203 * t177 - t260 * t79, t15 * t93 + t16 * t211 - t32 * t62 + t33 * t61 + t39 * t50 - t40 * t49 - t6 * t99 + t7 * t98 - t218, t7 * t50 + t33 * t16 + t6 * t49 + t32 * t15 + t198 * t219 + t104 * t227 - g(1) * (-t131 * t186 - t176 * t190) - g(2) * (t131 * t190 - t176 * t186), t13 * t272 + t57 * t9, -t10 * t57 + t13 * t54 - t14 * t272 - t56 * t9, t13 * t171 + t170 * t57, -t14 * t171 - t170 * t56, 0, -t44 * t54 + t73 * t10 + t17 * t56 + t63 * t14 + (-qJD(5) * t212 + t11 * t187 - t12 * t183) * t171 + t213 * t170 + t217 * t156, t44 * t272 + t73 * t9 + t17 * t57 + t63 * t13 - (qJD(5) * t213 + t11 * t183 + t12 * t187) * t171 - t212 * t170 - t217 * t155; 0, 0, 0, -t185 * t192 * t189, t242 * t192, t235, t234, qJDD(2), -g(3) * t189 + t185 * t200, g(3) * t185 + t189 * t200, -t248, t83, t59, t60, t175, -t223 * t177 + (-t115 * t241 + t188 * t175 - t177 * t239) * pkin(2) + t193, t245 * t177 + (t117 * t241 - t184 * t175 - t177 * t238) * pkin(2) + t194, -t110 * t40 + t111 * t39 + t211 * t250 + t252 * t93 + t216, t7 * t111 + t6 * t110 - t104 * (t168 - t266) - g(3) * t243 + t250 * t33 + t252 * t32 - t218 * (-pkin(2) * t185 - pkin(3) * t172), -t261, t8, t4, t5, t170, (t106 * t187 - t111 * t183) * t170 + t64 * t54 + (t251 * t183 + t253 * t187) * t171 + ((-t106 * t183 - t111 * t187) * t171 + t214) * qJD(5) + t199, -t64 * t272 + (-t106 * t170 - t2 + (qJD(5) * t111 - t253) * t171) * t183 + (-qJD(5) * t23 - t111 * t170 - t3 + (-qJD(5) * t106 + t251) * t171) * t187 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248, t83, t59, t60, t175, -t177 * t210 + t193, t224 * t177 + t194, -t37 * t93 - t38 * t211 + (t181 * t39 - t182 * t40) * pkin(3) + t216, -t32 * t37 - t33 * t38 + (t104 * t117 + t181 * t7 + t182 * t6 + t275) * pkin(3), -t261, t8, t4, t5, t170, t208 * t170 + t67 * t54 - (-t183 * t27 + t187 * t26) * t171 + (-t171 * t209 + t214) * qJD(5) + t199, -t209 * t170 - t187 * t3 - t183 * t2 - t67 * t272 + (t183 * t26 + t187 * t27) * t171 + (-t171 * t208 - t187 * t23) * qJD(5) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 ^ 2 - t93 ^ 2, -t211 * t33 - t32 * t93 + t198 - t217, 0, 0, 0, 0, 0, t10 + t257, t9 + t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t8, t4, t5, t170, t214 * t270 + t199, (-t25 * t171 - t2) * t183 + (-t23 * t270 - t3) * t187 + t206;];
tau_reg = t1;
