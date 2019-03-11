% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:37
% EndTime: 2019-03-09 02:39:45
% DurationCPUTime: 2.85s
% Computational Cost: add. (4432->376), mult. (9744->510), div. (0->0), fcn. (7112->18), ass. (0->197)
t172 = sin(pkin(10));
t245 = cos(pkin(10));
t181 = cos(qJ(3));
t162 = t181 * qJDD(2);
t178 = sin(qJ(3));
t173 = sin(pkin(9));
t150 = pkin(1) * t173 + pkin(7);
t138 = t150 * qJDD(1);
t188 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t138;
t237 = qJ(4) + t150;
t218 = t237 * qJD(1);
t198 = t218 * qJD(3);
t50 = qJDD(3) * pkin(3) - t178 * t188 - t181 * t198 + t162;
t56 = (qJDD(2) - t198) * t178 + t188 * t181;
t20 = -t172 * t56 + t245 * t50;
t19 = -qJDD(3) * pkin(4) + qJDD(5) - t20;
t167 = qJ(3) + pkin(10);
t157 = sin(t167);
t160 = cos(t167);
t168 = qJ(1) + pkin(9);
t158 = sin(t168);
t161 = cos(t168);
t215 = g(1) * t161 + g(2) * t158;
t191 = -g(3) * t160 + t157 * t215;
t268 = t19 - t191;
t220 = t245 * t181;
t144 = qJD(1) * t220;
t232 = qJD(1) * t178;
t118 = t172 * t232 - t144;
t115 = qJD(6) + t118;
t132 = t172 * t181 + t245 * t178;
t121 = t132 * qJD(1);
t171 = sin(pkin(11));
t174 = cos(pkin(11));
t107 = qJD(3) * t171 + t121 * t174;
t108 = t174 * qJD(3) - t121 * t171;
t177 = sin(qJ(6));
t180 = cos(qJ(6));
t57 = t107 * t177 - t108 * t180;
t270 = t115 * t57;
t133 = t171 * t180 + t174 * t177;
t125 = t133 * qJD(6);
t246 = t133 * t118 + t125;
t200 = -t107 * t180 - t108 * t177;
t269 = t115 * t200;
t223 = -g(1) * t158 + g(2) * t161;
t267 = t160 * t223;
t175 = cos(pkin(9));
t152 = -pkin(1) * t175 - pkin(2);
t164 = t181 * pkin(3);
t266 = t152 - t164;
t265 = -qJD(6) + t115;
t131 = t171 * t177 - t180 * t174;
t247 = t115 * t131;
t120 = t132 * qJD(3);
t228 = t178 * qJDD(1);
t90 = qJD(1) * t120 - qJDD(1) * t220 + t172 * t228;
t89 = qJDD(6) + t90;
t264 = t247 * t115 - t133 * t89;
t229 = qJD(1) * qJD(3);
t224 = t178 * t229;
t91 = qJD(3) * t144 + qJDD(1) * t132 - t172 * t224;
t78 = -t174 * qJDD(3) + t171 * t91;
t79 = qJDD(3) * t171 + t174 * t91;
t12 = -qJD(6) * t200 + t177 * t79 + t180 * t78;
t116 = t118 ^ 2;
t230 = qJD(6) * t180;
t231 = qJD(6) * t177;
t11 = -t107 * t231 + t108 * t230 - t177 * t78 + t180 * t79;
t195 = -t172 * t178 + t220;
t263 = -t11 * t195 - t120 * t200;
t262 = pkin(8) * t174;
t258 = g(3) * t157;
t256 = g(3) * t181;
t146 = pkin(3) * t172 + qJ(5);
t255 = pkin(8) + t146;
t21 = t172 * t50 + t245 * t56;
t18 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t21;
t187 = pkin(3) * t224 + qJDD(1) * t266 + qJDD(4);
t30 = pkin(4) * t90 - qJ(5) * t91 - qJD(5) * t121 + t187;
t7 = t171 * t30 + t174 * t18;
t123 = t195 * qJD(3);
t240 = t132 * t174;
t241 = t132 * t171;
t37 = t123 * t133 + t230 * t240 - t231 * t241;
t85 = t133 * t132;
t254 = -t37 * t115 - t85 * t89;
t110 = qJD(2) * t178 + t181 * t218;
t221 = t245 * t110;
t109 = t181 * qJD(2) - t178 * t218;
t253 = qJD(3) * pkin(3);
t99 = t109 + t253;
t55 = t172 * t99 + t221;
t47 = qJD(3) * qJ(5) + t55;
t117 = qJD(1) * t266 + qJD(4);
t69 = pkin(4) * t118 - qJ(5) * t121 + t117;
t24 = t171 * t69 + t174 * t47;
t96 = t172 * t110;
t62 = t245 * t109 - t96;
t82 = pkin(3) * t232 + pkin(4) * t121 + qJ(5) * t118;
t32 = t171 * t82 + t174 * t62;
t226 = t178 * t253;
t60 = pkin(4) * t120 - qJ(5) * t123 - qJD(5) * t132 + t226;
t219 = qJD(3) * t237;
t111 = qJD(4) * t181 - t178 * t219;
t112 = -qJD(4) * t178 - t181 * t219;
t68 = t245 * t111 + t172 * t112;
t29 = t171 * t60 + t174 * t68;
t84 = -pkin(4) * t195 - qJ(5) * t132 + t266;
t128 = t237 * t178;
t129 = t237 * t181;
t88 = -t172 * t128 + t245 * t129;
t39 = t171 * t84 + t174 * t88;
t252 = t121 * t57;
t251 = t121 * t200;
t249 = t171 * t90;
t248 = t174 * t90;
t244 = qJ(5) * t157;
t243 = t118 * t171;
t242 = t123 * t171;
t239 = t158 * t160;
t238 = t160 * t161;
t54 = t245 * t99 - t96;
t46 = -qJD(3) * pkin(4) + qJD(5) - t54;
t236 = -qJD(5) + t46;
t235 = qJDD(2) - g(3);
t153 = t164 + pkin(2);
t182 = cos(qJ(1));
t234 = t182 * pkin(1) + t161 * t153;
t169 = t178 ^ 2;
t233 = -t181 ^ 2 + t169;
t141 = qJD(1) * t152;
t227 = t181 * qJDD(1);
t6 = -t171 * t18 + t174 * t30;
t2 = pkin(5) * t90 - pkin(8) * t79 + t6;
t5 = -pkin(8) * t78 + t7;
t225 = -t177 * t5 + t180 * t2;
t23 = -t171 * t47 + t174 * t69;
t28 = -t171 * t68 + t174 * t60;
t31 = -t171 * t62 + t174 * t82;
t38 = -t171 * t88 + t174 * t84;
t61 = t109 * t172 + t221;
t67 = t111 * t172 - t245 * t112;
t87 = t245 * t128 + t129 * t172;
t216 = -t115 * t246 - t131 * t89;
t151 = -t245 * pkin(3) - pkin(4);
t179 = sin(qJ(1));
t213 = g(1) * t179 - g(2) * t182;
t212 = -t171 * t7 - t174 * t6;
t211 = -t171 * t6 + t174 * t7;
t210 = t177 * t2 + t180 * t5;
t36 = -t123 * t131 - t132 * t125;
t86 = t131 * t132;
t209 = -t115 * t36 + t86 * t89;
t176 = -qJ(4) - pkin(7);
t208 = -pkin(1) * t179 - t161 * t176;
t207 = -pkin(4) * t160 - t244;
t10 = pkin(5) * t118 - pkin(8) * t107 + t23;
t14 = pkin(8) * t108 + t24;
t3 = t10 * t180 - t14 * t177;
t4 = t10 * t177 + t14 * t180;
t206 = t12 * t195 - t120 * t57;
t205 = -t171 * t23 + t174 * t24;
t27 = -pkin(5) * t195 - pkin(8) * t240 + t38;
t34 = -pkin(8) * t241 + t39;
t204 = -t177 * t34 + t180 * t27;
t203 = t177 * t27 + t180 * t34;
t202 = -t118 * t123 - t132 * t90;
t201 = t107 * t171 + t108 * t174;
t126 = t255 * t171;
t197 = pkin(8) * t243 - qJD(5) * t174 + qJD(6) * t126 + t32;
t127 = t255 * t174;
t196 = pkin(5) * t121 + qJD(5) * t171 + qJD(6) * t127 + t118 * t262 + t31;
t193 = -qJD(1) * t141 - t138 + t215;
t192 = 0.2e1 * qJD(3) * t141 - qJDD(3) * t150;
t189 = t123 * t46 + t132 * t19 - t215;
t183 = qJD(3) ^ 2;
t186 = -0.2e1 * qJDD(1) * t152 - t150 * t183 - t223;
t184 = qJD(1) ^ 2;
t166 = pkin(11) + qJ(6);
t159 = cos(t166);
t156 = sin(t166);
t137 = qJDD(3) * t181 - t178 * t183;
t136 = qJDD(3) * t178 + t181 * t183;
t135 = -t174 * pkin(5) + t151;
t103 = t156 * t158 + t159 * t238;
t102 = -t156 * t238 + t158 * t159;
t101 = t156 * t161 - t159 * t239;
t100 = t156 * t239 + t159 * t161;
t66 = pkin(5) * t241 + t87;
t43 = pkin(5) * t242 + t67;
t40 = -pkin(5) * t243 + t61;
t35 = -pkin(5) * t108 + t46;
t17 = -pkin(8) * t242 + t29;
t13 = pkin(5) * t120 - t123 * t262 + t28;
t8 = t78 * pkin(5) + t19;
t1 = [qJDD(1), t213, g(1) * t182 + g(2) * t179 (t213 + (t173 ^ 2 + t175 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t169 + 0.2e1 * t181 * t224, 0.2e1 * t178 * t227 - 0.2e1 * t229 * t233, t136, t137, 0, t178 * t192 + t181 * t186, -t178 * t186 + t181 * t192, -t118 * t68 - t120 * t55 + t121 * t67 - t123 * t54 - t132 * t20 + t195 * t21 + t87 * t91 - t88 * t90 - t215, t21 * t88 + t55 * t68 - t20 * t87 - t54 * t67 + t187 * t266 + t117 * t226 - g(1) * (-t153 * t158 + t208) - g(2) * (-t158 * t176 + t234) -t108 * t67 + t28 * t118 + t23 * t120 + t171 * t189 - t174 * t267 - t195 * t6 + t38 * t90 + t87 * t78, t67 * t107 - t29 * t118 - t24 * t120 + t171 * t267 + t174 * t189 + t195 * t7 - t39 * t90 + t87 * t79, t108 * t29 - t107 * t28 - t38 * t79 - t39 * t78 - t223 * t157 + t212 * t132 + (-t171 * t24 - t174 * t23) * t123, t7 * t39 + t24 * t29 + t6 * t38 + t23 * t28 + t19 * t87 + t46 * t67 - g(1) * t208 - g(2) * (pkin(4) * t238 + t161 * t244 + t234) + (-g(1) * (-t153 + t207) + g(2) * t176) * t158, -t11 * t86 - t200 * t36, -t11 * t85 + t12 * t86 + t200 * t37 - t36 * t57, -t209 + t263, t206 + t254, t115 * t120 - t195 * t89 (t13 * t180 - t17 * t177) * t115 + t204 * t89 - t225 * t195 + t3 * t120 + t43 * t57 + t66 * t12 + t8 * t85 + t35 * t37 - g(1) * t101 - g(2) * t103 + (-t115 * t203 + t195 * t4) * qJD(6) -(t13 * t177 + t17 * t180) * t115 - t203 * t89 + t210 * t195 - t4 * t120 - t43 * t200 + t66 * t11 - t8 * t86 + t35 * t36 - g(1) * t100 - g(2) * t102 + (-t115 * t204 + t195 * t3) * qJD(6); 0, 0, 0, t235, 0, 0, 0, 0, 0, t137, -t136, t120 * t121 - t195 * t91 + t202, -t120 * t54 + t123 * t55 + t132 * t21 + t195 * t20 - g(3), -t108 * t120 + t171 * t202 - t195 * t78, t107 * t120 + t174 * t202 - t195 * t79 (t171 * t79 - t174 * t78) * t132 + t201 * t123, t120 * t46 + t205 * t123 + t211 * t132 - t19 * t195 - g(3), 0, 0, 0, 0, 0, -t206 + t254, t209 + t263; 0, 0, 0, 0, -t178 * t184 * t181, t233 * t184, t228, t227, qJDD(3), t178 * t193 + t162 - t256, -t178 * t235 + t181 * t193 (t55 - t61) * t121 + (-t54 + t62) * t118 + (-t172 * t90 - t245 * t91) * pkin(3), t54 * t61 - t55 * t62 + (t245 * t20 - t256 + t172 * t21 + (-qJD(1) * t117 + t215) * t178) * pkin(3), -t146 * t249 + t108 * t61 - t121 * t23 + t151 * t78 + (t171 * t236 - t31) * t118 - t268 * t174, -t146 * t248 - t107 * t61 + t121 * t24 + t151 * t79 + (t174 * t236 + t32) * t118 + t268 * t171, -t258 - t108 * t32 + t107 * t31 - t215 * t160 + (qJD(5) * t108 - t118 * t23 - t146 * t78 + t7) * t174 + (qJD(5) * t107 - t118 * t24 + t146 * t79 - t6) * t171, t19 * t151 - t24 * t32 - t23 * t31 - t46 * t61 - g(3) * (t164 - t207) + t211 * t146 + t205 * qJD(5) + t215 * (pkin(3) * t178 + pkin(4) * t157 - qJ(5) * t160) t11 * t133 + t200 * t247, -t11 * t131 - t12 * t133 + t200 * t246 + t247 * t57, t251 - t264, t216 + t252, -t115 * t121 (-t126 * t180 - t127 * t177) * t89 + t135 * t12 + t8 * t131 - t3 * t121 - t40 * t57 + t246 * t35 + (t177 * t197 - t180 * t196) * t115 + t191 * t159 -(-t126 * t177 + t127 * t180) * t89 + t135 * t11 + t8 * t133 + t4 * t121 + t40 * t200 - t247 * t35 + (t177 * t196 + t180 * t197) * t115 - t191 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 ^ 2 - t116, t118 * t55 + t121 * t54 + t187 + t223, t108 * t121 - t116 * t171 + t248, -t107 * t121 - t116 * t174 - t249, t118 * t201 - t171 * t78 - t174 * t79, t118 * t205 - t121 * t46 - t212 + t223, 0, 0, 0, 0, 0, t216 - t252, t251 + t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t118 + t78, t108 * t118 + t79, -t107 ^ 2 - t108 ^ 2, t107 * t23 - t24 * t108 + t268, 0, 0, 0, 0, 0, t12 - t269, t11 - t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200 * t57, t200 ^ 2 - t57 ^ 2, t11 + t270, -t12 - t269, t89, -g(1) * t102 + g(2) * t100 + t156 * t258 + t200 * t35 + t265 * t4 + t225, g(1) * t103 - g(2) * t101 + t159 * t258 + t265 * t3 + t35 * t57 - t210;];
tau_reg  = t1;
