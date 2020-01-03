% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:34
% DurationCPUTime: 3.29s
% Computational Cost: add. (2512->385), mult. (6285->471), div. (0->0), fcn. (3721->4), ass. (0->195)
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t123 = cos(qJ(2));
t204 = qJD(1) * t123;
t103 = -qJD(3) + t204;
t121 = sin(qJ(2));
t205 = qJD(1) * t121;
t186 = t122 * t205;
t203 = qJD(2) * t120;
t78 = t186 + t203;
t228 = t103 * t78;
t187 = t120 * t205;
t198 = t122 * qJD(2);
t76 = t187 - t198;
t229 = t103 * t76;
t197 = qJD(1) * qJD(2);
t179 = t123 * t197;
t200 = qJD(3) * t120;
t183 = t121 * t200;
t196 = qJD(2) * qJD(3);
t50 = qJD(1) * t183 + (-t179 - t196) * t122;
t199 = qJD(3) * t122;
t182 = t121 * t199;
t201 = qJD(2) * t123;
t264 = t120 * t201 + t182;
t51 = qJD(1) * t264 + t120 * t196;
t270 = (t50 - t229) * t122 + (t51 - t228) * t120;
t108 = t121 * t197;
t241 = pkin(3) + qJ(5);
t267 = t241 * t108;
t85 = -pkin(2) * t123 - pkin(7) * t121 - pkin(1);
t139 = t85 * qJD(1);
t266 = -pkin(6) * t108 + qJD(3) * t139;
t160 = pkin(2) * t121 - pkin(7) * t123;
t83 = t160 * qJD(2);
t68 = qJD(1) * t83;
t111 = pkin(6) * t204;
t89 = qJD(2) * pkin(7) + t111;
t173 = t266 * t120 - t122 * t68 + t89 * t199;
t42 = t120 * t139 + t122 * t89;
t142 = -t103 * t42 - t173;
t174 = -t120 * t68 - t266 * t122 + t89 * t200;
t41 = t120 * t89 - t122 * t139;
t265 = -t103 * t41 - t174;
t148 = pkin(4) * t78 + t41;
t209 = qJD(4) + t148;
t221 = t51 * t122;
t222 = t50 * t120;
t126 = (qJD(3) * (t120 * t76 - t122 * t78) - t221 + t222) * t121 - (t120 * t78 + t122 * t76) * t201;
t170 = -t78 + t203;
t181 = t103 * t199;
t214 = t122 * t123;
t24 = (t103 * t214 + t121 * t170) * qJD(1) - t181;
t169 = t76 + t198;
t216 = t120 * t123;
t263 = (-t103 * t216 + t121 * t169) * qJD(1) + t103 * t200;
t117 = t121 ^ 2;
t149 = qJD(1) * t117 - t103 * t123;
t5 = (t120 * t149 + t121 * t76) * qJD(2) - t121 * t181 - t123 * t51;
t4 = (t121 * t78 + t122 * t149) * qJD(2) + t103 * t183 + t123 * t50;
t262 = -0.2e1 * t197;
t261 = pkin(6) * (t121 * t198 + t123 * t200);
t73 = t76 ^ 2;
t74 = t78 ^ 2;
t192 = -t74 + t73;
t101 = qJ(4) * t108;
t92 = qJD(4) * t103;
t260 = t101 - t92;
t167 = pkin(3) * t108;
t13 = -t167 + t173;
t35 = t103 * qJ(4) - t42;
t259 = -t103 * t35 + t13;
t31 = t51 + t228;
t100 = t103 ^ 2;
t258 = -t74 - t100;
t256 = t120 * qJD(4) + t111 + (t120 * t204 - t200) * pkin(3);
t250 = pkin(4) + pkin(7);
t249 = t76 * pkin(4);
t248 = pkin(6) * t120;
t247 = pkin(7) * t103;
t141 = pkin(6) * t179 + qJ(4) * t50 - qJD(4) * t78;
t8 = pkin(3) * t51 + t141;
t246 = t120 * t8;
t245 = t122 * t8;
t110 = pkin(6) * t205;
t234 = qJD(2) * pkin(2);
t88 = t110 - t234;
t143 = -t78 * qJ(4) + t88;
t21 = t241 * t76 + t143;
t244 = t21 * t78;
t37 = pkin(3) * t76 + t143;
t243 = t37 * t78;
t242 = t76 * t78;
t188 = -pkin(3) - t248;
t130 = pkin(4) * t214 + (-qJ(5) + t188) * t121;
t80 = t160 * qJD(1);
t224 = t122 * t80;
t91 = t250 * t122;
t240 = -qJD(1) * t130 + qJD(3) * t91 + t224;
t176 = pkin(6) * t122 - qJ(4);
t132 = -pkin(4) * t216 - t121 * t176;
t65 = t120 * t80;
t239 = qJD(1) * t132 + t250 * t200 + t65;
t218 = qJ(4) * t122;
t152 = qJ(5) * t120 - t218;
t140 = t152 * t123;
t238 = qJD(1) * t140 - qJD(3) * t152 + t122 * qJD(5) + t256;
t237 = qJ(4) * t199 - t204 * t218 + t256;
t236 = t120 * t83 + t85 * t199;
t235 = qJ(4) * t51;
t16 = t241 * t103 + t209;
t233 = t103 * t16;
t227 = t120 * t88;
t223 = t122 * t88;
t220 = t76 * qJ(4);
t107 = pkin(6) * t214;
t57 = t120 * t85 + t107;
t217 = t120 * t121;
t215 = t121 * t122;
t125 = qJD(1) ^ 2;
t213 = t123 * t125;
t124 = qJD(2) ^ 2;
t212 = t124 * t121;
t211 = t124 * t123;
t210 = -qJD(4) - t41;
t29 = t42 - t249;
t208 = -qJD(5) - t29;
t207 = pkin(3) * t217 + t121 * pkin(6);
t206 = -t123 ^ 2 + t117;
t202 = qJD(2) * t121;
t195 = t120 * t247;
t194 = t122 * t247;
t106 = pkin(6) * t216;
t193 = pkin(7) * t198;
t112 = pkin(6) * t201;
t189 = t121 * t213;
t180 = t103 * t205;
t177 = -t120 * qJ(4) - pkin(2);
t175 = -t108 + t242;
t56 = t122 * t85 - t106;
t172 = pkin(1) * t262;
t171 = t264 * pkin(3) + qJ(4) * t183 + t112;
t165 = qJD(3) * t107 - t122 * t83 + t85 * t200;
t164 = t123 * t108;
t163 = t92 + t174;
t48 = qJ(4) * t123 - t57;
t162 = -qJD(4) * t123 + t236;
t161 = t188 * t121;
t159 = (qJD(3) * t76 - t50) * pkin(7);
t158 = (qJD(3) * t78 - t51) * pkin(7);
t34 = pkin(3) * t103 - t210;
t157 = t120 * t35 + t122 * t34;
t156 = -t120 * t42 + t122 * t41;
t147 = -pkin(4) * t50 + t173;
t146 = -pkin(4) * t51 - t174;
t1 = qJD(5) * t76 + t241 * t51 + t141;
t145 = t1 * t120 + t199 * t21;
t144 = -t1 * t122 + t200 * t21;
t134 = -t21 * t76 + t146;
t131 = t50 + t229;
t14 = -t122 * t228 - t222;
t15 = -t120 * t229 - t221;
t129 = t147 - t267;
t12 = t76 * t182 + (t121 * t51 + t201 * t76) * t120;
t11 = t78 * t123 * t198 + (-t122 * t50 - t200 * t78) * t121;
t116 = t123 * pkin(3);
t99 = 0.2e1 * t101;
t90 = t250 * t120;
t84 = -pkin(3) * t122 + t177;
t67 = -t241 * t122 + t177;
t60 = -qJ(4) * t215 + t207;
t55 = (-t103 - t204) * t202;
t53 = -pkin(6) * t186 + t65;
t52 = pkin(6) * t187 + t224;
t49 = t116 - t56;
t47 = t121 * t152 + t207;
t45 = pkin(3) * t78 + t220;
t44 = qJD(1) * t161 - t224;
t43 = t176 * t205 - t65;
t39 = -pkin(4) * t217 - t48;
t38 = t123 * qJ(5) + t106 + t116 + (pkin(4) * t121 - t85) * t122;
t32 = t241 * t78 + t220;
t27 = t202 * t248 - t165;
t26 = t236 - t261;
t25 = (-qJ(4) * t201 - qJD(4) * t121) * t122 + t171;
t22 = qJD(2) * t161 + t165;
t20 = qJD(5) - t35 - t249;
t19 = -qJ(4) * t202 - t162 + t261;
t10 = qJD(2) * t140 + (qJD(5) * t120 + (qJ(5) * qJD(3) - qJD(4)) * t122) * t121 + t171;
t9 = (-pkin(4) * t215 - t106) * qJD(3) + t132 * qJD(2) + t162;
t7 = -t101 + t163;
t6 = -pkin(4) * t183 + qJD(2) * t130 + qJD(5) * t123 + t165;
t3 = t146 + t260;
t2 = qJD(5) * t103 + t129;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t164, t206 * t262, t211, -0.2e1 * t164, -t212, 0, -pkin(6) * t211 + t121 * t172, pkin(6) * t212 + t123 * t172, 0, 0, t11, t126, t4, t12, -t5, t55, -t103 * t27 + t123 * t173 + (pkin(6) * t51 + t199 * t88) * t121 + ((pkin(6) * t76 + t227) * t123 + (-t41 + (t56 + t106) * qJD(1)) * t121) * qJD(2), t103 * t26 - t123 * t174 + (-pkin(6) * t50 - t200 * t88) * t121 + ((pkin(6) * t78 + t223) * t123 + (-t42 + (-t57 + t107) * qJD(1)) * t121) * qJD(2), -t26 * t76 - t27 * t78 + t50 * t56 - t51 * t57 + t156 * t201 + (t120 * t174 + t122 * t173 + (-t120 * t41 - t122 * t42) * qJD(3)) * t121, -t174 * t57 - t173 * t56 + t26 * t42 - t27 * t41 + (t88 + t110) * t112, t55, -t4, t5, t11, t126, t12, t19 * t76 + t22 * t78 + t48 * t51 - t49 * t50 + t157 * t201 + (t120 * t7 + t122 * t13 + (-t120 * t34 + t122 * t35) * qJD(3)) * t121, -t103 * t22 - t25 * t76 - t51 * t60 + (-t203 * t37 - t13) * t123 + (-t37 * t199 - t246 + (qJD(1) * t49 + t34) * qJD(2)) * t121, t103 * t19 - t25 * t78 + t50 * t60 + (-t198 * t37 + t7) * t123 + (t37 * t200 - t245 + (-qJD(1) * t48 - t35) * qJD(2)) * t121, t13 * t49 + t19 * t35 + t22 * t34 + t25 * t37 + t48 * t7 + t60 * t8, t55, t5, t4, t12, -t126, t11, -t38 * t50 - t39 * t51 + t6 * t78 - t76 * t9 + (-t120 * t20 + t122 * t16) * t201 + (-t120 * t3 + t122 * t2 + (-t120 * t16 - t122 * t20) * qJD(3)) * t121, -t10 * t78 - t103 * t9 + t47 * t50 + (-t198 * t21 - t3) * t123 + ((qJD(1) * t39 + t20) * qJD(2) + t144) * t121, t10 * t76 + t103 * t6 + t47 * t51 + (t203 * t21 + t2) * t123 + ((-qJD(1) * t38 - t16) * qJD(2) + t145) * t121, t1 * t47 + t10 * t21 + t16 * t6 + t2 * t38 + t20 * t9 + t3 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t206 * t125, 0, t189, 0, 0, t125 * pkin(1) * t121, pkin(1) * t213, 0, 0, t14, -t270, t24, t15, t263, t180, -pkin(2) * t51 + t103 * t52 + (t194 + t227) * qJD(3) + ((-pkin(7) * t203 + t41) * t121 + (-pkin(6) * t169 - t227) * t123) * qJD(1), pkin(2) * t50 - t103 * t53 + (-t195 + t223) * qJD(3) + ((t42 - t193) * t121 + (pkin(6) * t170 - t223) * t123) * qJD(1), t52 * t78 + t53 * t76 + (t158 + t265) * t122 + (t159 - t142) * t120, t41 * t52 - t42 * t53 + (-t88 - t234) * t111 + (qJD(3) * t156 + t120 * t173 - t122 * t174) * pkin(7), t180, -t24, -t263, t14, -t270, t15, -t43 * t76 - t44 * t78 + (-t103 * t34 + t158 - t7) * t122 + (t159 + t259) * t120, t103 * t44 + t245 - t51 * t84 + t237 * t76 + (-t120 * t37 - t194) * qJD(3) + (-t121 * t34 + (pkin(7) * t202 + t123 * t37) * t120) * qJD(1), -t103 * t43 - t246 + t50 * t84 + t237 * t78 + (-t122 * t37 + t195) * qJD(3) + (t37 * t214 + (t35 + t193) * t121) * qJD(1), -t34 * t44 - t35 * t43 + t8 * t84 - t237 * t37 + (qJD(3) * t157 + t120 * t13 - t122 * t7) * pkin(7), t180, -t263, t24, t15, t270, t14, -t50 * t90 - t51 * t91 + t240 * t78 + t239 * t76 + (t3 - t233) * t122 + (t103 * t20 + t2) * t120, t50 * t67 + t238 * t78 + t239 * t103 + (t21 * t214 + (qJD(2) * t91 - t20) * t121) * qJD(1) - t145, t51 * t67 - t238 * t76 + t240 * t103 + (-t21 * t216 + (-qJD(2) * t90 + t16) * t121) * qJD(1) + t144, t1 * t67 + t16 * t240 + t2 * t90 - t20 * t239 - t21 * t238 + t3 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, -t192, -t131, -t242, -t31, t108, -t78 * t88 + t142, t76 * t88 - t265, 0, 0, t108, t131, t31, t242, -t192, -t242, pkin(3) * t50 - t235 + (-t35 - t42) * t78 + (t34 + t210) * t76, t45 * t76 - t142 - 0.2e1 * t167 + t243, t103 * t210 - t37 * t76 + t45 * t78 - t163 + t99, -pkin(3) * t13 - qJ(4) * t7 + t210 * t35 - t34 * t42 - t37 * t45, t108, t31, -t131, -t242, t192, t242, -t235 + t241 * t50 + (t20 + t208) * t78 + (t16 - t209) * t76, -t103 * t148 + t32 * t78 + t134 - 0.2e1 * t92 + t99, -t244 - t32 * t76 + (-0.2e1 * qJD(5) - t29) * t103 + 0.2e1 * t267 - t147, qJ(4) * t3 + t16 * t208 - t2 * t241 + t20 * t209 - t21 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t175, t258, t243 + t259, 0, 0, 0, 0, 0, 0, -t131, t258, t175, t244 + (qJD(5) + t20) * t103 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t108 + t242, -t73 - t100, t134 - t233 + t260;];
tauc_reg = t17;
