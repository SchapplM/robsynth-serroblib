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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:06
% DurationCPUTime: 2.59s
% Computational Cost: add. (2649->337), mult. (6599->470), div. (0->0), fcn. (4900->14), ass. (0->192)
t154 = sin(pkin(9));
t155 = sin(pkin(8));
t156 = cos(pkin(9));
t162 = cos(qJ(3));
t210 = qJDD(1) * t162;
t198 = t155 * t210;
t159 = sin(qJ(3));
t211 = qJDD(1) * t159;
t111 = t154 * t162 + t156 * t159;
t257 = t111 * qJD(3);
t50 = t155 * (qJD(1) * t257 + t154 * t211) - t156 * t198;
t149 = t155 ^ 2;
t253 = 0.2e1 * t149;
t157 = cos(pkin(8));
t221 = qJD(1) * t157;
t133 = -qJD(3) + t221;
t127 = -qJD(5) + t133;
t158 = sin(qJ(5));
t161 = cos(qJ(5));
t174 = qJD(1) * t111;
t78 = t155 * t174;
t242 = t161 * t78;
t222 = qJD(1) * t155;
t204 = t159 * t222;
t189 = t154 * t204;
t220 = qJD(1) * t162;
t203 = t155 * t220;
t82 = t156 * t203 - t189;
t37 = t158 * t82 + t242;
t244 = t127 * t37;
t216 = qJD(5) * t158;
t199 = t155 * t211;
t212 = qJD(1) * qJD(3);
t201 = t162 * t212;
t265 = t155 * t201 + t199;
t49 = qJD(3) * t189 - t154 * t198 - t156 * t265;
t6 = -qJD(5) * t242 + t158 * t49 - t161 * t50 - t82 * t216;
t268 = t6 - t244;
t182 = -t158 * t78 + t161 * t82;
t267 = t182 * t37;
t236 = qJDD(1) * pkin(1);
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t259 = g(1) * t160 - g(2) * t163;
t179 = -qJDD(2) + t236 + t259;
t245 = t127 * t182;
t7 = qJD(5) * t182 - t158 * t50 - t161 * t49;
t266 = -t7 - t245;
t213 = qJD(1) * qJD(2);
t214 = qJ(2) * qJDD(1);
t178 = t213 + t214;
t264 = t182 ^ 2 - t37 ^ 2;
t252 = pkin(7) * t78;
t238 = qJ(2) * t162;
t132 = t157 * t238;
t117 = -pkin(2) * t157 - pkin(6) * t155 - pkin(1);
t98 = qJD(1) * t117 + qJD(2);
t59 = -qJ(4) * t204 + qJD(1) * t132 + t159 * t98;
t243 = t156 * t59;
t233 = t155 * t162;
t206 = qJ(4) * t233;
t239 = qJ(2) * t159;
t208 = t157 * t239;
t172 = -t206 - t208;
t89 = t162 * t98;
t58 = qJD(1) * t172 + t89;
t48 = -pkin(3) * t133 + t58;
t26 = t154 * t48 + t243;
t12 = t26 - t252;
t11 = t12 * t216;
t151 = qJ(3) + pkin(9);
t147 = qJ(5) + t151;
t139 = cos(t147);
t246 = g(3) * t155;
t99 = pkin(3) * t204 + qJ(2) * t222 + qJD(4);
t57 = pkin(4) * t78 + t99;
t138 = sin(t147);
t232 = t157 * t160;
t74 = t138 * t163 - t139 * t232;
t231 = t157 * t163;
t76 = t138 * t160 + t139 * t231;
t263 = g(1) * t76 - g(2) * t74 + t139 * t246 + t57 * t37 + t11;
t260 = t178 * t157;
t227 = t162 * t163;
t230 = t159 * t160;
t102 = t157 * t230 + t227;
t228 = t160 * t162;
t229 = t159 * t163;
t104 = -t157 * t229 + t228;
t258 = -g(1) * t104 + g(2) * t102;
t181 = t154 * t159 - t156 * t162;
t256 = t181 * qJD(3);
t255 = qJD(5) + t127;
t209 = t157 * qJDD(1);
t131 = -qJDD(3) + t209;
t217 = qJD(4) * t155;
t219 = qJD(2) * t157;
t171 = -t159 * t219 - t162 * t217;
t234 = t155 * t159;
t207 = qJ(4) * t234;
t237 = qJD(3) * t98;
t97 = qJDD(1) * t117 + qJDD(2);
t88 = t162 * t97;
t15 = -t159 * t237 - pkin(3) * t131 + t88 + t172 * qJDD(1) + ((-t132 + t207) * qJD(3) + t171) * qJD(1);
t166 = qJD(3) * t172 - t159 * t217;
t202 = t162 * t213;
t218 = qJD(3) * t162;
t194 = qJDD(1) * t132 + t157 * t202 + t159 * t97 + t98 * t218;
t21 = -qJ(4) * t199 + qJD(1) * t166 + t194;
t4 = t156 * t15 - t154 * t21;
t2 = -pkin(4) * t131 + pkin(7) * t50 + t4;
t5 = t154 * t15 + t156 * t21;
t3 = pkin(7) * t49 + t5;
t205 = -t158 * t3 + t161 * t2;
t73 = t138 * t232 + t139 * t163;
t75 = -t138 * t231 + t139 * t160;
t254 = -g(1) * t75 + g(2) * t73 + t138 * t246 - t57 * t182 + t205;
t150 = t157 ^ 2;
t251 = pkin(7) * t82;
t250 = pkin(3) * t154;
t226 = t117 * t218 + t162 * t219;
t46 = t166 + t226;
t47 = (-t132 + (qJ(4) * t155 - t117) * t159) * qJD(3) + t171;
t19 = t154 * t47 + t156 * t46;
t53 = t154 * t59;
t29 = t156 * t58 - t53;
t109 = t162 * t117;
t63 = -t206 + t109 + (-pkin(3) - t239) * t157;
t225 = t159 * t117 + t132;
t68 = -t207 + t225;
t32 = t154 * t63 + t156 * t68;
t241 = t157 * t174 - t257;
t240 = t181 * t221 - t256;
t164 = qJD(1) ^ 2;
t235 = t149 * t164;
t107 = (pkin(3) * t218 + qJD(2)) * t155;
t112 = pkin(3) * t234 + t155 * qJ(2);
t224 = t149 + t150;
t153 = t162 ^ 2;
t223 = t159 ^ 2 - t153;
t215 = qJD(3) + t133;
t18 = -t154 * t46 + t156 * t47;
t25 = t156 * t48 - t53;
t28 = -t154 * t58 - t243;
t31 = -t154 * t68 + t156 * t63;
t195 = t224 * t164;
t193 = qJD(1) * t215;
t192 = t131 + t209;
t191 = 0.2e1 * t224;
t190 = qJD(3) * t208;
t187 = g(1) * t163 + g(2) * t160;
t185 = qJD(5) * t111 - t241;
t184 = -qJD(5) * t181 + t240;
t10 = -pkin(4) * t133 + t25 - t251;
t183 = -t10 * t158 - t161 * t12;
t95 = t111 * t155;
t96 = t181 * t155;
t51 = -t158 * t96 + t161 * t95;
t52 = -t158 * t95 - t161 * t96;
t64 = pkin(3) * t265 + t155 * t178 + qJDD(4);
t180 = qJD(3) * (t133 + t221);
t141 = pkin(3) * t156 + pkin(4);
t177 = t141 * t158 + t161 * t250;
t176 = t141 * t161 - t158 * t250;
t169 = -t133 ^ 2 - t235;
t167 = t191 * t213 - t187;
t145 = cos(t151);
t144 = sin(t151);
t140 = pkin(3) * t159 + qJ(2);
t122 = -qJDD(5) + t131;
t105 = t157 * t227 + t230;
t103 = -t157 * t228 + t229;
t94 = (qJ(4) + pkin(6)) * t155 + pkin(1) + (pkin(3) * t162 + pkin(2)) * t157;
t93 = t144 * t160 + t145 * t231;
t92 = -t144 * t231 + t145 * t160;
t91 = t144 * t163 - t145 * t232;
t90 = t144 * t232 + t145 * t163;
t85 = t155 * t257;
t81 = t155 * t256;
t67 = pkin(3) * t203 + pkin(4) * t82;
t65 = pkin(4) * t95 + t112;
t60 = -pkin(4) * t81 + t107;
t30 = -pkin(4) * t49 + t64;
t27 = -pkin(7) * t95 + t32;
t24 = -pkin(4) * t157 + pkin(7) * t96 + t31;
t23 = qJD(5) * t52 - t158 * t85 - t161 * t81;
t22 = -qJD(5) * t51 + t158 * t81 - t161 * t85;
t17 = t29 - t251;
t16 = t28 + t252;
t9 = pkin(7) * t81 + t19;
t8 = pkin(7) * t85 + t18;
t1 = [qJDD(1), t259, t187, (t179 + t236) * t157, t191 * t214 + t167, t179 * pkin(1) + (t224 * t214 + t167) * qJ(2), (qJDD(1) * t153 - 0.2e1 * t159 * t201) * t149, (-t159 * t210 + t223 * t212) * t253, (t159 * t180 - t192 * t162) * t155, (t192 * t159 + t162 * t180) * t155, t131 * t157, -g(1) * t103 - g(2) * t105 - t109 * t131 - t88 * t157 + (t133 * t157 + (t253 + t150) * qJD(1)) * qJ(2) * t218 + (qJD(3) * t117 * t133 + t178 * t253 + (qJ(2) * t131 + qJD(2) * t133 + t237 + t260) * t157) * t159, (-t190 + t226) * t133 + t225 * t131 + (-qJD(1) * t190 + t194) * t157 - g(1) * t102 - g(2) * t104 + (t202 + (-t159 * t212 + t210) * qJ(2)) * t253, -g(1) * t91 - g(2) * t93 + t107 * t78 - t112 * t49 - t131 * t31 - t133 * t18 - t157 * t4 + t64 * t95 - t81 * t99, -g(1) * t90 - g(2) * t92 + t107 * t82 - t112 * t50 + t131 * t32 + t133 * t19 + t157 * t5 - t64 * t96 - t85 * t99, t155 * t259 - t18 * t82 - t19 * t78 + t25 * t85 + t26 * t81 + t31 * t50 + t32 * t49 + t4 * t96 - t5 * t95, t5 * t32 + t26 * t19 + t4 * t31 + t25 * t18 + t64 * t112 + t99 * t107 - g(1) * (t140 * t163 - t160 * t94) - g(2) * (t140 * t160 + t163 * t94), t182 * t22 + t52 * t6, -t182 * t23 - t22 * t37 - t51 * t6 - t52 * t7, -t122 * t52 - t127 * t22 - t157 * t6, t122 * t51 + t127 * t23 + t157 * t7, t122 * t157, -(-t158 * t9 + t161 * t8) * t127 - (-t158 * t27 + t161 * t24) * t122 - t205 * t157 + t60 * t37 + t65 * t7 + t30 * t51 + t57 * t23 - g(1) * t74 - g(2) * t76 + (-(-t158 * t24 - t161 * t27) * t127 - t183 * t157) * qJD(5), -g(1) * t73 - g(2) * t75 - t11 * t157 + t57 * t22 + t30 * t52 + t60 * t182 + t65 * t6 + ((-qJD(5) * t27 + t8) * t127 + t24 * t122 + t2 * t157) * t158 + ((qJD(5) * t24 + t9) * t127 + t27 * t122 + (qJD(5) * t10 + t3) * t157) * t161; 0, 0, 0, -t209, -t195, -qJ(2) * t195 - t179, 0, 0, 0, 0, 0, -t131 * t162 + t159 * t169, t131 * t159 + t162 * t169, t131 * t181 - t133 * t241 - t222 * t78, t111 * t131 + t133 * t240 - t222 * t82, t111 * t49 - t181 * t50 - t240 * t78 - t241 * t82, t111 * t5 - t181 * t4 - t222 * t99 + t240 * t26 + t241 * t25 - t259, 0, 0, 0, 0, 0, -(-t111 * t158 - t161 * t181) * t122 - t37 * t222 + (t158 * t184 + t161 * t185) * t127, (t111 * t161 - t158 * t181) * t122 - t182 * t222 + (-t158 * t185 + t161 * t184) * t127; 0, 0, 0, 0, 0, 0, t162 * t159 * t235, -t223 * t235, (-t159 * t193 + t210) * t155, (-t162 * t193 - t211) * t155, -t131, t88 + (-t157 * t193 - t235) * t238 + (-t215 * t98 + t246 - t260) * t159 + t258, g(3) * t233 + g(1) * t105 - g(2) * t103 - t89 * t133 + (t215 * t221 + t235) * t239 - t194, t144 * t246 - g(1) * t92 + g(2) * t90 + t133 * t28 - t82 * t99 + (-t131 * t156 - t203 * t78) * pkin(3) + t4, t145 * t246 + g(1) * t93 - g(2) * t91 - t133 * t29 + t78 * t99 + (t131 * t154 - t203 * t82) * pkin(3) - t5, (t26 + t28) * t82 + (-t25 + t29) * t78 + (t154 * t49 + t156 * t50) * pkin(3), -t25 * t28 - t26 * t29 + (t5 * t154 + t4 * t156 + (g(3) * t159 - t220 * t99) * t155 + t258) * pkin(3), t267, t264, t268, t266, -t122, -t176 * t122 + (-t158 * t17 + t16 * t161) * t127 - t67 * t37 + (t127 * t177 + t183) * qJD(5) + t254, t177 * t122 - t161 * t3 - t158 * t2 - (t158 * t16 + t161 * t17) * t127 - t67 * t182 + (-t161 * t10 + t127 * t176) * qJD(5) + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t82 - t49, t133 * t78 - t50, -t78 ^ 2 - t82 ^ 2, g(3) * t157 - t155 * t187 + t25 * t82 + t26 * t78 + t64, 0, 0, 0, 0, 0, t7 - t245, t6 + t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t264, t268, t266, -t122, t255 * t183 + t254, (t12 * t127 - t2) * t158 + (-t255 * t10 - t3) * t161 + t263;];
tau_reg = t1;
