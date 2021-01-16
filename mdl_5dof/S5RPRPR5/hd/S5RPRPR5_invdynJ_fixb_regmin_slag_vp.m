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
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 11:55:12
% EndTime: 2021-01-15 11:55:25
% DurationCPUTime: 2.73s
% Computational Cost: add. (2649->337), mult. (6599->470), div. (0->0), fcn. (4900->14), ass. (0->193)
t155 = sin(pkin(9));
t156 = sin(pkin(8));
t157 = cos(pkin(9));
t164 = cos(qJ(3));
t213 = qJDD(1) * t164;
t200 = t156 * t213;
t161 = sin(qJ(3));
t214 = qJDD(1) * t161;
t110 = t155 * t164 + t157 * t161;
t262 = t110 * qJD(3);
t50 = (qJD(1) * t262 + t155 * t214) * t156 - t157 * t200;
t150 = t156 ^ 2;
t257 = 0.2e1 * t150;
t158 = cos(pkin(8));
t224 = qJD(1) * t158;
t133 = -qJD(3) + t224;
t127 = -qJD(5) + t133;
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t176 = qJD(1) * t110;
t78 = t156 * t176;
t246 = t163 * t78;
t225 = qJD(1) * t156;
t206 = t161 * t225;
t191 = t155 * t206;
t223 = qJD(1) * t164;
t205 = t156 * t223;
t82 = t157 * t205 - t191;
t37 = t160 * t82 + t246;
t248 = t127 * t37;
t219 = qJD(5) * t160;
t201 = t156 * t214;
t215 = qJD(1) * qJD(3);
t203 = t164 * t215;
t270 = t156 * t203 + t201;
t49 = qJD(3) * t191 - t155 * t200 - t270 * t157;
t6 = -qJD(5) * t246 + t160 * t49 - t163 * t50 - t82 * t219;
t273 = t6 - t248;
t184 = -t160 * t78 + t163 * t82;
t272 = t184 * t37;
t241 = qJDD(1) * pkin(1);
t162 = sin(qJ(1));
t165 = cos(qJ(1));
t264 = -g(2) * t165 - g(3) * t162;
t181 = -qJDD(2) + t241 + t264;
t249 = t127 * t184;
t7 = t184 * qJD(5) - t160 * t50 - t163 * t49;
t271 = -t7 - t249;
t216 = qJD(1) * qJD(2);
t217 = qJ(2) * qJDD(1);
t180 = t216 + t217;
t269 = t184 ^ 2 - t37 ^ 2;
t256 = pkin(7) * t78;
t243 = qJ(2) * t164;
t132 = t158 * t243;
t116 = -pkin(2) * t158 - pkin(6) * t156 - pkin(1);
t97 = qJD(1) * t116 + qJD(2);
t59 = -qJ(4) * t206 + qJD(1) * t132 + t161 * t97;
t247 = t157 * t59;
t238 = t156 * t164;
t209 = qJ(4) * t238;
t244 = qJ(2) * t161;
t211 = t158 * t244;
t174 = -t209 - t211;
t89 = t164 * t97;
t58 = qJD(1) * t174 + t89;
t48 = -pkin(3) * t133 + t58;
t26 = t155 * t48 + t247;
t12 = t26 - t256;
t11 = t12 * t219;
t152 = qJ(3) + pkin(9);
t147 = qJ(5) + t152;
t139 = cos(t147);
t253 = g(1) * t156;
t98 = pkin(3) * t206 + qJ(2) * t225 + qJD(4);
t57 = pkin(4) * t78 + t98;
t138 = sin(t147);
t231 = t165 * t138;
t237 = t158 * t162;
t74 = t139 * t237 - t231;
t236 = t158 * t165;
t76 = t138 * t162 + t139 * t236;
t268 = g(2) * t74 - g(3) * t76 + t139 * t253 + t57 * t37 + t11;
t265 = t180 * t158;
t232 = t164 * t165;
t235 = t161 * t162;
t101 = -t158 * t235 - t232;
t233 = t162 * t164;
t234 = t161 * t165;
t103 = t158 * t234 - t233;
t263 = -g(2) * t101 - g(3) * t103;
t183 = t155 * t161 - t157 * t164;
t261 = t183 * qJD(3);
t260 = qJD(5) + t127;
t212 = t158 * qJDD(1);
t131 = -qJDD(3) + t212;
t220 = qJD(4) * t156;
t222 = qJD(2) * t158;
t173 = -t161 * t222 - t164 * t220;
t239 = t156 * t161;
t210 = qJ(4) * t239;
t242 = qJD(3) * t97;
t96 = qJDD(1) * t116 + qJDD(2);
t88 = t164 * t96;
t15 = -t161 * t242 - pkin(3) * t131 + t88 + t174 * qJDD(1) + ((-t132 + t210) * qJD(3) + t173) * qJD(1);
t168 = qJD(3) * t174 - t161 * t220;
t204 = t164 * t216;
t221 = qJD(3) * t164;
t196 = qJDD(1) * t132 + t158 * t204 + t161 * t96 + t97 * t221;
t21 = -qJ(4) * t201 + qJD(1) * t168 + t196;
t4 = t157 * t15 - t155 * t21;
t2 = -pkin(4) * t131 + pkin(7) * t50 + t4;
t5 = t155 * t15 + t157 * t21;
t3 = pkin(7) * t49 + t5;
t207 = -t160 * t3 + t163 * t2;
t73 = -t138 * t237 - t139 * t165;
t75 = -t139 * t162 + t158 * t231;
t259 = -g(2) * t73 - g(3) * t75 + t138 * t253 - t57 * t184 + t207;
t258 = (pkin(3) * t164 + pkin(2)) * t158 + t156 * (qJ(4) + pkin(6)) + pkin(1);
t151 = t158 ^ 2;
t255 = pkin(7) * t82;
t254 = pkin(3) * t155;
t230 = t116 * t221 + t164 * t222;
t46 = t168 + t230;
t47 = (-t132 + (qJ(4) * t156 - t116) * t161) * qJD(3) + t173;
t19 = t155 * t47 + t157 * t46;
t53 = t155 * t59;
t29 = t157 * t58 - t53;
t108 = t164 * t116;
t63 = -t209 + t108 + (-pkin(3) - t244) * t158;
t229 = t161 * t116 + t132;
t68 = -t210 + t229;
t32 = t155 * t63 + t157 * t68;
t250 = -t158 * t176 + t262;
t245 = t183 * t224 - t261;
t166 = qJD(1) ^ 2;
t240 = t150 * t166;
t106 = (pkin(3) * t221 + qJD(2)) * t156;
t111 = pkin(3) * t239 + t156 * qJ(2);
t227 = t150 + t151;
t154 = t164 ^ 2;
t226 = t161 ^ 2 - t154;
t218 = qJD(3) + t133;
t18 = -t155 * t46 + t157 * t47;
t25 = t157 * t48 - t53;
t28 = -t155 * t58 - t247;
t31 = -t155 * t68 + t157 * t63;
t197 = t227 * t166;
t195 = qJD(1) * t218;
t194 = t131 + t212;
t193 = 0.2e1 * t227;
t192 = qJD(3) * t211;
t189 = qJD(5) * t110 + t250;
t187 = g(2) * t162 - g(3) * t165;
t186 = -qJD(5) * t183 + t245;
t10 = -pkin(4) * t133 + t25 - t255;
t185 = -t10 * t160 - t12 * t163;
t94 = t110 * t156;
t95 = t183 * t156;
t51 = -t160 * t95 + t163 * t94;
t52 = -t160 * t94 - t163 * t95;
t64 = t270 * pkin(3) + t180 * t156 + qJDD(4);
t182 = qJD(3) * (t133 + t224);
t141 = pkin(3) * t157 + pkin(4);
t179 = t141 * t160 + t163 * t254;
t178 = t141 * t163 - t160 * t254;
t171 = -t133 ^ 2 - t240;
t169 = t193 * t216 - t187;
t145 = cos(t152);
t144 = sin(t152);
t140 = -pkin(3) * t161 - qJ(2);
t122 = -qJDD(5) + t131;
t104 = t158 * t232 + t235;
t102 = t158 * t233 - t234;
t93 = t144 * t162 + t145 * t236;
t92 = t144 * t236 - t145 * t162;
t91 = -t144 * t165 + t145 * t237;
t90 = -t144 * t237 - t145 * t165;
t85 = t156 * t262;
t81 = t156 * t261;
t67 = pkin(3) * t205 + pkin(4) * t82;
t65 = pkin(4) * t94 + t111;
t60 = -pkin(4) * t81 + t106;
t30 = -pkin(4) * t49 + t64;
t27 = -pkin(7) * t94 + t32;
t24 = -pkin(4) * t158 + pkin(7) * t95 + t31;
t23 = qJD(5) * t52 - t160 * t85 - t163 * t81;
t22 = -t51 * qJD(5) + t160 * t81 - t163 * t85;
t17 = t29 - t255;
t16 = t28 + t256;
t9 = pkin(7) * t81 + t19;
t8 = pkin(7) * t85 + t18;
t1 = [qJDD(1), t264, t187, (t181 + t241) * t158, t193 * t217 + t169, t181 * pkin(1) + (t227 * t217 + t169) * qJ(2), (qJDD(1) * t154 - 0.2e1 * t161 * t203) * t150, (-t161 * t213 + t226 * t215) * t257, (t161 * t182 - t194 * t164) * t156, (t194 * t161 + t164 * t182) * t156, t131 * t158, -g(2) * t104 - g(3) * t102 - t108 * t131 - t88 * t158 + (t133 * t158 + (t257 + t151) * qJD(1)) * qJ(2) * t221 + (qJD(3) * t116 * t133 + t180 * t257 + (qJ(2) * t131 + qJD(2) * t133 + t242 + t265) * t158) * t161, (-t192 + t230) * t133 + t229 * t131 + (-qJD(1) * t192 + t196) * t158 + g(2) * t103 - g(3) * t101 + (t204 + (-t161 * t215 + t213) * qJ(2)) * t257, -g(2) * t93 - g(3) * t91 + t106 * t78 - t111 * t49 - t131 * t31 - t133 * t18 - t158 * t4 + t64 * t94 - t81 * t98, g(2) * t92 - g(3) * t90 + t106 * t82 - t111 * t50 + t131 * t32 + t133 * t19 + t158 * t5 - t64 * t95 - t85 * t98, t156 * t264 - t18 * t82 - t19 * t78 + t25 * t85 + t26 * t81 + t31 * t50 + t32 * t49 + t4 * t95 - t5 * t94, t5 * t32 + t26 * t19 + t4 * t31 + t25 * t18 + t64 * t111 + t98 * t106 - g(2) * (-t140 * t162 + t258 * t165) - g(3) * (t140 * t165 + t258 * t162), t184 * t22 + t52 * t6, -t184 * t23 - t22 * t37 - t51 * t6 - t52 * t7, -t122 * t52 - t127 * t22 - t158 * t6, t122 * t51 + t127 * t23 + t158 * t7, t122 * t158, -(-t160 * t9 + t163 * t8) * t127 - (-t160 * t27 + t163 * t24) * t122 - t207 * t158 + t60 * t37 + t65 * t7 + t30 * t51 + t57 * t23 - g(2) * t76 - g(3) * t74 + (-(-t160 * t24 - t163 * t27) * t127 - t185 * t158) * qJD(5), g(2) * t75 - g(3) * t73 - t11 * t158 + t57 * t22 + t30 * t52 + t60 * t184 + t65 * t6 + ((-qJD(5) * t27 + t8) * t127 + t24 * t122 + t2 * t158) * t160 + ((qJD(5) * t24 + t9) * t127 + t27 * t122 + (qJD(5) * t10 + t3) * t158) * t163; 0, 0, 0, -t212, -t197, -qJ(2) * t197 - t181, 0, 0, 0, 0, 0, -t131 * t164 + t161 * t171, t131 * t161 + t164 * t171, t131 * t183 + t250 * t133 - t78 * t225, t110 * t131 + t245 * t133 - t82 * t225, t110 * t49 - t183 * t50 - t245 * t78 + t250 * t82, t110 * t5 - t183 * t4 - t98 * t225 + t245 * t26 - t250 * t25 - t264, 0, 0, 0, 0, 0, -(-t110 * t160 - t163 * t183) * t122 - t37 * t225 + (t160 * t186 + t163 * t189) * t127, (t110 * t163 - t160 * t183) * t122 - t184 * t225 + (-t160 * t189 + t163 * t186) * t127; 0, 0, 0, 0, 0, 0, t164 * t161 * t240, -t226 * t240, (-t161 * t195 + t213) * t156, (-t164 * t195 - t214) * t156, -t131, t88 + (-t158 * t195 - t240) * t243 + (-t218 * t97 + t253 - t265) * t161 + t263, g(1) * t238 + g(2) * t102 - g(3) * t104 - t89 * t133 + (t218 * t224 + t240) * t244 - t196, t144 * t253 - g(2) * t90 - g(3) * t92 + t133 * t28 - t82 * t98 + (-t131 * t157 - t78 * t205) * pkin(3) + t4, t145 * t253 + g(2) * t91 - g(3) * t93 - t133 * t29 + t78 * t98 + (t131 * t155 - t82 * t205) * pkin(3) - t5, (t26 + t28) * t82 + (-t25 + t29) * t78 + (t155 * t49 + t157 * t50) * pkin(3), -t25 * t28 - t26 * t29 + (t5 * t155 + t4 * t157 + (g(1) * t161 - t98 * t223) * t156 + t263) * pkin(3), t272, t269, t273, t271, -t122, -t178 * t122 + (t16 * t163 - t160 * t17) * t127 - t67 * t37 + (t127 * t179 + t185) * qJD(5) + t259, t179 * t122 - t163 * t3 - t160 * t2 - (t16 * t160 + t163 * t17) * t127 - t67 * t184 + (-t163 * t10 + t127 * t178) * qJD(5) + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t82 - t49, t133 * t78 - t50, -t78 ^ 2 - t82 ^ 2, g(1) * t158 - t156 * t187 + t25 * t82 + t26 * t78 + t64, 0, 0, 0, 0, 0, t7 - t249, t6 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t269, t273, t271, -t122, t185 * t260 + t259, (t12 * t127 - t2) * t160 + (-t10 * t260 - t3) * t163 + t268;];
tau_reg = t1;
