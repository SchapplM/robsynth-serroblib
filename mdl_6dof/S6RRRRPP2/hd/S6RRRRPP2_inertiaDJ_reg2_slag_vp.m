% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:40
% EndTime: 2019-03-09 20:51:51
% DurationCPUTime: 4.49s
% Computational Cost: add. (4687->375), mult. (10612->539), div. (0->0), fcn. (9388->6), ass. (0->185)
t151 = sin(qJ(4));
t149 = t151 ^ 2;
t154 = cos(qJ(4));
t150 = t154 ^ 2;
t260 = t149 + t150;
t155 = cos(qJ(2));
t238 = cos(qJ(3));
t187 = t238 * t155;
t152 = sin(qJ(3));
t153 = sin(qJ(2));
t217 = t152 * t153;
t113 = -t187 + t217;
t259 = 0.2e1 * t113;
t145 = t151 * qJ(5);
t242 = pkin(4) + pkin(5);
t169 = t242 * t154 + t145;
t205 = t154 * qJD(5);
t258 = -qJD(4) * t169 + t205;
t216 = t152 * t155;
t114 = t153 * t238 + t216;
t219 = t114 * t151;
t144 = qJD(4) * t154;
t186 = t114 * t144;
t176 = qJD(2) * t187;
t183 = qJD(3) * t238;
t248 = qJD(2) + qJD(3);
t66 = -t155 * t183 + t248 * t217 - t176;
t226 = t151 * t66;
t44 = t186 - t226;
t67 = t248 * t114;
t14 = t113 * t44 + t67 * t219;
t257 = -0.2e1 * t14;
t218 = t114 * t154;
t143 = qJD(4) * t151;
t224 = t154 * t66;
t41 = t114 * t143 + t224;
t13 = -0.2e1 * t67 * t218 + t259 * t41;
t251 = t67 * qJ(5) + t113 * qJD(5);
t256 = qJ(6) * t144 + t151 * qJD(6);
t178 = pkin(2) * t183;
t255 = t260 * t178;
t249 = (t149 - t150) * qJD(4);
t223 = qJ(5) * t154;
t245 = t242 * t151 - t223;
t253 = 0.4e1 * t114;
t252 = 0.2e1 * t251;
t209 = qJD(3) * t152;
t197 = pkin(2) * t209;
t142 = t151 * qJD(5);
t93 = pkin(4) * t143 - qJ(5) * t144 - t142;
t76 = -pkin(5) * t143 - t93;
t70 = t76 - t197;
t60 = t70 * t154;
t138 = -pkin(2) * t238 - pkin(3);
t175 = t154 * pkin(4) + t145;
t101 = t138 - t175;
t146 = t154 * pkin(5);
t85 = t146 - t101;
t250 = -t85 * t143 + t60;
t207 = qJD(6) * t154;
t206 = t153 * qJD(2);
t200 = pkin(2) * t206;
t165 = pkin(3) * t67 + pkin(9) * t66 + t200;
t241 = -pkin(8) - pkin(7);
t190 = qJD(2) * t241;
t116 = t153 * t190;
t123 = t241 * t153;
t124 = t241 * t155;
t37 = -t116 * t238 - t123 * t183 - t124 * t209 - t190 * t216;
t139 = -pkin(2) * t155 - pkin(1);
t52 = pkin(3) * t113 - pkin(9) * t114 + t139;
t74 = t152 * t123 - t124 * t238;
t11 = -t52 * t143 - t74 * t144 + t151 * t37 + t154 * t165;
t58 = t67 * pkin(4);
t8 = -t58 - t11;
t247 = -t41 * qJ(6) + t114 * t207 - t8;
t38 = t74 * qJD(3) + t152 * t116 - t241 * t176;
t244 = qJD(4) * t175 - t205;
t110 = t114 ^ 2;
t192 = t151 * t224;
t19 = 0.2e1 * t110 * t249 + t192 * t253;
t243 = -0.2e1 * t249;
t157 = 0.2e1 * qJD(5);
t240 = pkin(9) * t67;
t73 = -t238 * t123 - t152 * t124;
t24 = -t245 * t114 - t73;
t9 = t114 * t258 + t245 * t66 - t38;
t239 = t24 * t144 + t9 * t151;
t237 = pkin(9) * t113;
t236 = t38 * t73;
t7 = t9 * t154;
t235 = pkin(9) - qJ(6);
t234 = t73 * t144 + t38 * t151;
t233 = t85 * t144 + t70 * t151;
t31 = t151 * t52 + t154 * t74;
t201 = pkin(3) + t175;
t102 = t146 + t201;
t83 = t102 * t144;
t232 = t76 * t151 + t83;
t77 = t93 + t197;
t231 = -t77 - t93;
t230 = pkin(2) * qJD(3);
t137 = pkin(2) * t152 + pkin(9);
t229 = t137 * t67;
t228 = t149 * t66;
t227 = t150 * t66;
t225 = t152 * t73;
t222 = qJD(4) * t24;
t174 = pkin(4) * t151 - t223;
t36 = t114 * t174 + t73;
t221 = qJD(4) * t36;
t220 = t113 * t137;
t215 = -qJ(6) + t137;
t214 = t255 * t137;
t213 = t255 * pkin(9);
t212 = t138 * t144 + t151 * t197;
t204 = t155 * qJD(2);
t203 = -0.2e1 * pkin(1) * qJD(2);
t202 = t67 * t259;
t21 = t113 * qJ(5) + t31;
t199 = pkin(3) * t143;
t198 = pkin(3) * t144;
t196 = pkin(9) * t143;
t195 = pkin(9) * t144;
t32 = t36 * t143;
t61 = t73 * t143;
t189 = t151 * t238;
t188 = t154 * t238;
t185 = t151 * t144;
t184 = t153 * t204;
t68 = t151 * t74;
t30 = t154 * t52 - t68;
t104 = t215 * t151;
t179 = qJ(6) * t143 - t207;
t80 = t137 * t143 - t154 * t178;
t64 = t179 - t80;
t181 = qJD(4) * t104 + t64;
t121 = t235 * t151;
t92 = t179 - t196;
t180 = qJD(4) * t121 + t92;
t177 = t110 * t185;
t22 = -t113 * pkin(4) - t30;
t173 = t151 * t21 - t154 * t22;
t172 = t151 * t31 + t154 * t30;
t168 = -t114 * t138 + t220;
t167 = t138 * t143 - t154 * t197;
t42 = t113 * t144 + t151 * t67;
t39 = t113 * t143 - t154 * t67;
t10 = t143 * t74 - t52 * t144 - t151 * t165 + t154 * t37;
t166 = t114 * t93 + t201 * t66 - t240;
t164 = (-t113 * t238 + t114 * t152) * qJD(3);
t12 = t244 * t114 - t174 * t66 + t38;
t163 = -t12 + (-t114 * t201 - t237) * qJD(4);
t161 = -t12 + (t101 * t114 - t220) * qJD(4);
t160 = -qJ(6) * t226 + t256 * t114 - t10;
t81 = t137 * t144 + t151 * t178;
t5 = -t10 + t251;
t1 = -qJD(4) * t173 + t8 * t151 + t5 * t154;
t2 = -qJD(4) * t172 - t10 * t154 - t11 * t151;
t159 = -t101 * t66 - t113 * t178 + t114 * t77 - t229;
t158 = pkin(2) * t164 - t138 * t66 - t229;
t148 = qJ(5) * t157;
t130 = -0.2e1 * t185;
t129 = 0.2e1 * t185;
t122 = t235 * t154;
t106 = t122 * t143;
t105 = t215 * t154;
t99 = t201 * t143;
t94 = t195 - t256;
t90 = t260 * t238 * t230;
t84 = t105 * t143;
t82 = t101 * t143;
t79 = 0.2e1 * t90;
t72 = t76 * t154;
t65 = t81 - t256;
t29 = -0.2e1 * t114 * t227 - 0.2e1 * t177;
t28 = -0.2e1 * t114 * t228 + 0.2e1 * t177;
t27 = t114 * t249 + t192;
t18 = t185 * t253 + t227 - t228;
t17 = qJ(6) * t219 + t21;
t16 = t17 * t143;
t15 = t68 + (-qJ(6) * t114 - t52) * t154 - t242 * t113;
t4 = t160 + t251;
t3 = -pkin(5) * t67 - t247;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t184, 0.2e1 * (-t153 ^ 2 + t155 ^ 2) * qJD(2), 0, -0.2e1 * t184, 0, 0, t153 * t203, t155 * t203, 0, 0, -0.2e1 * t114 * t66, 0.2e1 * t113 * t66 - 0.2e1 * t114 * t67, 0, t202, 0, 0, 0.2e1 * t113 * t200 + 0.2e1 * t139 * t67, 0.2e1 * t114 * t200 - 0.2e1 * t139 * t66, 0.2e1 * t113 * t37 + 0.2e1 * t114 * t38 - 0.2e1 * t66 * t73 - 0.2e1 * t67 * t74, 0.2e1 * t139 * t200 - 0.2e1 * t37 * t74 + 0.2e1 * t236, t29, t19, -t13, t28, t257, t202, 0.2e1 * t11 * t113 + 0.2e1 * t114 * t234 - 0.2e1 * t226 * t73 + 0.2e1 * t30 * t67, -0.2e1 * t73 * t224 + 0.2e1 * t10 * t113 - 0.2e1 * t31 * t67 + 0.2e1 * (t38 * t154 - t61) * t114, 0.2e1 * t172 * t66 + 0.2e1 * (t10 * t151 - t11 * t154 + (t151 * t30 - t154 * t31) * qJD(4)) * t114, -0.2e1 * t10 * t31 + 0.2e1 * t11 * t30 + 0.2e1 * t236, t29, -t13, -t19, t202, 0.2e1 * t14, t28, -0.2e1 * t36 * t226 - 0.2e1 * t113 * t8 - 0.2e1 * t22 * t67 + 0.2e1 * (t12 * t151 + t144 * t36) * t114, 0.2e1 * t173 * t66 + 0.2e1 * (-t151 * t5 + t154 * t8 + (-t151 * t22 - t154 * t21) * qJD(4)) * t114, 0.2e1 * t36 * t224 + 0.2e1 * t113 * t5 + 0.2e1 * t21 * t67 + 0.2e1 * (-t12 * t154 + t32) * t114, 0.2e1 * t12 * t36 + 0.2e1 * t21 * t5 + 0.2e1 * t22 * t8, t29, -t19, t13, t28, t257, t202, -0.2e1 * t113 * t3 - 0.2e1 * t114 * t239 - 0.2e1 * t15 * t67 + 0.2e1 * t226 * t24, -0.2e1 * t24 * t224 + 0.2e1 * t113 * t4 + 0.2e1 * t17 * t67 + 0.2e1 * (-t143 * t24 + t7) * t114, 0.2e1 * (t15 * t154 - t151 * t17) * t66 + 0.2e1 * (t151 * t4 - t154 * t3 + (t15 * t151 + t154 * t17) * qJD(4)) * t114, 0.2e1 * t15 * t3 + 0.2e1 * t17 * t4 + 0.2e1 * t24 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, -t206, 0, -pkin(7) * t204, pkin(7) * t206, 0, 0, 0, 0, -t66, 0, -t67, 0, -t38, t37 (-t152 * t67 + t238 * t66 + t164) * pkin(2) (-t238 * t38 - t152 * t37 + (t238 * t74 + t225) * qJD(3)) * pkin(2), -t27, -t18, t42, t27, -t39, 0, t61 + (-qJD(4) * t168 - t38) * t154 + t158 * t151, t143 * t168 + t154 * t158 + t234, t2, t38 * t138 + (t188 * t31 - t189 * t30 + t225) * t230 + t2 * t137, -t27, t42, t18, 0, t39, t27, t151 * t159 + t154 * t161 + t32, t1, t161 * t151 + (-t159 - t221) * t154, t12 * t101 + t36 * t77 + (t188 * t21 + t189 * t22) * t230 + t1 * t137, -t27, t18, -t42, t27, -t39, 0, -t85 * t186 - t104 * t67 - t113 * t65 + t7 + (-t114 * t70 + t66 * t85 - t222) * t151, t105 * t67 + t113 * t64 + t250 * t114 - t85 * t224 + t239, t16 + (-t105 * t66 + t114 * t181 - t3) * t151 + (t104 * t66 - t114 * t65 - t4 + (t105 * t114 - t15) * qJD(4)) * t154, t104 * t3 + t105 * t4 + t15 * t65 + t17 * t64 + t24 * t70 + t85 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t197, -0.2e1 * t178, 0, 0, t129, t243, 0, t130, 0, 0, 0.2e1 * t167, 0.2e1 * t212, t79, 0.2e1 * t138 * t197 + 0.2e1 * t214, t129, 0, -t243, 0, 0, t130, -0.2e1 * t154 * t77 + 0.2e1 * t82, t79, -0.2e1 * t101 * t144 - 0.2e1 * t77 * t151, 0.2e1 * t101 * t77 + 0.2e1 * t214, t129, -t243, 0, t130, 0, 0, 0.2e1 * t250, 0.2e1 * t233, -0.2e1 * t65 * t151 - 0.2e1 * t154 * t181 + 0.2e1 * t84, 0.2e1 * t104 * t65 + 0.2e1 * t105 * t64 + 0.2e1 * t70 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, -t67, 0, -t38, t37, 0, 0, -t27, -t18, t42, t27, -t39, 0, t61 + (pkin(3) * t66 - t240) * t151 + (-t38 + (-pkin(3) * t114 - t237) * qJD(4)) * t154, pkin(3) * t41 + pkin(9) * t39 + t234, t2, -pkin(3) * t38 + pkin(9) * t2, -t27, t42, t18, 0, t39, t27, t151 * t166 + t154 * t163 + t32, t1, t163 * t151 + (-t166 - t221) * t154, pkin(9) * t1 - t12 * t201 + t36 * t93, -t27, t18, -t42, t27, -t39, 0, -t114 * t83 - t113 * t94 - t121 * t67 + t7 + (t102 * t66 - t114 * t76 - t222) * t151, -t102 * t41 + t113 * t92 + t122 * t67 + t218 * t76 + t239, t16 + (t114 * t180 - t122 * t66 - t3) * t151 + (-t114 * t94 + t121 * t66 - t4 + (t114 * t122 - t15) * qJD(4)) * t154, t102 * t9 + t121 * t3 + t122 * t4 + t15 * t94 + t17 * t92 + t24 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, -t178, 0, 0, t129, t243, 0, t130, 0, 0, t167 - t199, -t198 + t212, t90, -pkin(3) * t197 + t213, t129, 0, -t243, 0, 0, t130, t154 * t231 + t82 - t99, t90, t231 * t151 + (-t101 + t201) * t144, t101 * t93 - t201 * t77 + t213, t129, -t243, 0, t130, 0, 0, t60 + t72 + (-t102 - t85) * t143, t232 + t233, t106 + t84 + (-t65 - t94) * t151 + (-t64 - t92 + (-t104 - t121) * qJD(4)) * t154, t102 * t70 + t104 * t94 + t105 * t92 + t121 * t65 + t122 * t64 + t76 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t243, 0, t130, 0, 0, -0.2e1 * t199, -0.2e1 * t198, 0, 0, t129, 0, -t243, 0, 0, t130, -0.2e1 * t154 * t93 - 0.2e1 * t99, 0, 0.2e1 * t144 * t201 - 0.2e1 * t93 * t151, -0.2e1 * t201 * t93, t129, -t243, 0, t130, 0, 0, -0.2e1 * t102 * t143 + 0.2e1 * t72, 0.2e1 * t232, -0.2e1 * t94 * t151 - 0.2e1 * t154 * t180 + 0.2e1 * t106, 0.2e1 * t102 * t76 + 0.2e1 * t121 * t94 + 0.2e1 * t122 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, -t44, t67, t11, t10, 0, 0, 0, -t41, 0, t67, t44, 0, -t8 + t58, t175 * t66 + (qJD(4) * t174 - t142) * t114, -t10 + t252, -pkin(4) * t8 + qJ(5) * t5 + qJD(5) * t21, 0, 0, t41, 0, -t44, t67 (pkin(5) + t242) * t67 + t247, t160 + t252, -t169 * t66 + (-t245 * qJD(4) + t142) * t114, qJ(5) * t4 + qJD(5) * t17 - t242 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, -t143, 0, -t81, t80, 0, 0, 0, t144, 0, 0, t143, 0, -t81, -t244, -t80 (-pkin(4) * t189 + qJ(5) * t188) * t230 - t244 * t137, 0, 0, -t144, 0, -t143, 0, -t65, t64, -t258, qJ(5) * t64 + qJD(5) * t105 - t242 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, -t143, 0, -t195, t196, 0, 0, 0, t144, 0, 0, t143, 0, -t195, -t244, -t196, -t244 * pkin(9), 0, 0, -t144, 0, -t143, 0, -t94, t92, -t258, qJ(5) * t92 + qJD(5) * t122 - t242 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t148, 0, 0, 0, 0, 0, 0, 0, t157, 0, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t41, 0, t8, 0, 0, 0, 0, 0, 0, -t67, 0, t41, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, t81, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, t195, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t41, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t144, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t144, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
