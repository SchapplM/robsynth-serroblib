% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP2
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
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:01
% EndTime: 2019-12-31 20:52:06
% DurationCPUTime: 1.62s
% Computational Cost: add. (1459->245), mult. (2646->270), div. (0->0), fcn. (1902->4), ass. (0->189)
t143 = sin(qJ(3));
t145 = cos(qJ(3));
t137 = t145 * pkin(4);
t135 = t143 * qJ(4);
t197 = -t145 * pkin(3) - t135;
t188 = pkin(2) - t197;
t86 = t137 + t188;
t225 = t86 * t143;
t146 = cos(qJ(2));
t237 = t146 * pkin(1);
t124 = -pkin(2) - t237;
t85 = t124 + t197;
t71 = t137 - t85;
t229 = t71 * t143;
t236 = t229 / 0.2e1 + t225 / 0.2e1;
t147 = -pkin(3) - pkin(4);
t210 = t147 * t143;
t211 = t145 * qJ(4);
t94 = t211 + t210;
t76 = t94 * t145;
t247 = t236 - t76;
t239 = t143 * pkin(3);
t105 = -t211 + t239;
t215 = t105 * t145;
t72 = t85 * t143;
t95 = t188 * t143;
t234 = t72 / 0.2e1 - t95 / 0.2e1;
t246 = t234 - t215;
t192 = t145 * qJD(4);
t245 = qJD(3) * t197 + t192;
t190 = qJD(1) + qJD(2);
t189 = -pkin(3) / 0.2e1 - pkin(4) / 0.2e1;
t54 = t211 + (t147 / 0.2e1 + t189) * t143;
t244 = t190 * t54;
t180 = pkin(2) / 0.2e1 - t124 / 0.2e1;
t243 = t143 * t180;
t141 = t143 ^ 2;
t142 = t145 ^ 2;
t118 = t142 + t141;
t242 = t190 * t118;
t119 = t142 - t141;
t241 = t190 * t119;
t240 = pkin(2) * t145;
t144 = sin(qJ(2));
t238 = t144 * pkin(1);
t3 = t71 * t94;
t6 = t94 * t86;
t232 = pkin(1) * qJD(1);
t231 = pkin(1) * qJD(2);
t230 = t3 * qJD(1);
t228 = t71 * t145;
t227 = t85 * t105;
t226 = t85 * t145;
t224 = t86 * t145;
t123 = pkin(7) + t238;
t134 = t143 * qJ(5);
t88 = t143 * t123 - t134;
t223 = t88 * t143;
t136 = t145 * qJ(5);
t89 = t145 * t123 - t136;
t222 = t89 * t145;
t108 = t118 * qJD(5);
t167 = t118 * t146;
t87 = pkin(1) * t167;
t202 = t87 * qJD(2);
t221 = t108 - t202;
t132 = t143 * qJD(4);
t220 = t54 * qJD(3) + t132;
t219 = t188 * t145;
t104 = t143 * pkin(7) - t134;
t218 = t104 * t143;
t217 = t105 * t188;
t216 = t105 * t143;
t106 = t145 * pkin(7) - t136;
t214 = t106 * t145;
t28 = t222 + t223;
t12 = (-t144 * t71 + t146 * t28) * pkin(1);
t213 = t12 * qJD(1);
t212 = t124 * t145;
t22 = (t123 * t167 + t144 * t85) * pkin(1);
t209 = t22 * qJD(1);
t75 = t94 * t143;
t25 = t75 + t228;
t208 = t25 * qJD(1);
t26 = t76 - t229;
t207 = t26 * qJD(1);
t206 = t28 * qJD(1);
t32 = t216 + t226;
t205 = t32 * qJD(1);
t33 = -t72 + t215;
t204 = t33 * qJD(1);
t203 = t87 * qJD(1);
t201 = t89 * qJD(3);
t185 = t144 * t232;
t115 = t143 * t185;
t130 = t141 * qJD(4);
t200 = t115 + t130;
t116 = t145 * t185;
t120 = t143 * t192;
t199 = t116 + t120;
t187 = t144 * t231;
t117 = t143 * t187;
t198 = t130 - t117;
t196 = t106 * qJD(3);
t195 = t143 * qJD(1);
t194 = t143 * qJD(2);
t193 = t143 * qJD(3);
t133 = t145 * qJD(3);
t191 = t145 * qJD(5);
t186 = pkin(7) * t193;
t184 = pkin(7) * t133;
t183 = -t238 / 0.2e1;
t182 = -t237 / 0.2e1;
t181 = t237 / 0.2e1;
t179 = -t86 / 0.2e1 - t71 / 0.2e1;
t178 = t71 * t195;
t177 = qJD(1) * t227;
t176 = t85 * t195;
t175 = -t188 / 0.2e1 + t85 / 0.2e1;
t174 = -t115 + t198;
t173 = t124 * t195;
t172 = qJD(1) * t212;
t171 = t123 * t193;
t170 = t123 * t133;
t169 = t211 / 0.2e1;
t168 = pkin(1) * t190;
t166 = (-t147 * t145 + t135) * qJD(3) - t192;
t103 = t190 * t145;
t165 = t145 * t187;
t149 = (t169 + t210 / 0.2e1) * t237;
t1 = t179 * t94 + t149;
t164 = -t1 * qJD(1) + t6 * qJD(2);
t29 = t75 + t224;
t113 = t145 * t181;
t8 = t145 * t179 + t113 - t75;
t163 = -t8 * qJD(1) + t29 * qJD(2);
t30 = t76 - t225;
t112 = t143 * t182;
t7 = t112 + t247;
t162 = t7 * qJD(1) - t30 * qJD(2);
t38 = t214 + t218;
t4 = t183 + (t106 / 0.2e1 + t89 / 0.2e1) * t145 + (t104 / 0.2e1 + t88 / 0.2e1) * t143;
t161 = -t4 * qJD(1) - t38 * qJD(2);
t16 = t112 - t246;
t36 = t95 + t215;
t160 = t16 * qJD(1) + t36 * qJD(2);
t17 = t145 * t175 + t113 + t216;
t35 = t216 - t219;
t159 = t17 * qJD(1) + t35 * qJD(2);
t158 = qJD(3) * t105 - t132;
t157 = t120 - t165;
t49 = t112 + t243;
t156 = pkin(2) * t194 + t49 * qJD(1);
t114 = t145 * t182;
t50 = t145 * t180 + t114;
t155 = t50 * qJD(1) + qJD(2) * t240;
t111 = t143 * t181;
t20 = t111 - t236;
t154 = t20 * qJD(1) - t194 * t86;
t150 = (t169 - t239 / 0.2e1) * t237;
t13 = -t105 * t175 + t150;
t153 = t13 * qJD(1) + qJD(2) * t217;
t23 = t111 + t234;
t152 = t23 * qJD(1) - t188 * t194;
t151 = -t116 + t157;
t140 = qJ(4) * qJD(4);
t139 = qJD(3) * qJ(4);
t131 = t143 * qJD(5);
t121 = t143 * t133;
t109 = t119 * qJD(3);
t102 = t190 * t143;
t101 = t190 * t141;
t81 = t143 * t103;
t73 = -t210 / 0.2e1 + t189 * t143;
t70 = t73 * qJD(3);
t69 = t73 * qJD(5);
t52 = -t240 / 0.2e1 + t212 / 0.2e1 + t114;
t51 = t112 - t243;
t45 = t54 * qJD(5);
t24 = t111 - t234;
t21 = t111 + t236;
t19 = t112 + t246;
t18 = -t216 - t226 / 0.2e1 + t113 + t219 / 0.2e1;
t15 = -t217 / 0.2e1 + t227 / 0.2e1 + t150;
t10 = t112 - t247;
t9 = t75 + t224 / 0.2e1 + t228 / 0.2e1 + t113;
t5 = -t214 / 0.2e1 - t222 / 0.2e1 - t218 / 0.2e1 - t223 / 0.2e1 + t183;
t2 = t6 / 0.2e1 + t3 / 0.2e1 + t149;
t11 = [0, 0, 0, 0, -t187, -t146 * t231, t121, t109, 0, 0, 0, t124 * t193 - t165, t124 * t133 + t117, -t33 * qJD(3) + t157, t202, -t32 * qJD(3) + t198, t22 * qJD(2) + t158 * t85, t26 * qJD(3) + t157, t25 * qJD(3) + t198, t221, t12 * qJD(2) + t3 * qJD(3) - t28 * qJD(5) + t132 * t71; 0, 0, 0, 0, -t144 * t168, -t146 * t168, t121, t109, 0, 0, 0, t51 * qJD(3) - t116 - t165, t52 * qJD(3) + t115 + t117, t19 * qJD(3) + t151, t203 + t202, t18 * qJD(3) + t174, t209 + t15 * qJD(3) + t24 * qJD(4) + (pkin(7) * t167 - t144 * t188) * t231, t10 * qJD(3) + t151, t9 * qJD(3) + t174, -t203 + t221, t213 + t2 * qJD(3) + t21 * qJD(4) + t5 * qJD(5) + (-t144 * t86 + t146 * t38) * t231; 0, 0, 0, 0, 0, 0, t81, t241, t133, -t193, 0, t51 * qJD(2) - t170 + t173, t52 * qJD(2) + t171 + t172, t19 * qJD(2) - t170 - t204, t245, t18 * qJD(2) - t171 - t205, t15 * qJD(2) + t123 * t245 + t177, t10 * qJD(2) - t201 + t207, t9 * qJD(2) - t88 * qJD(3) + t208, t166, t230 + t2 * qJD(2) + (-t88 * qJ(4) + t89 * t147) * qJD(3) + t89 * qJD(4) + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t133, t101, t24 * qJD(2) + t170 - t176, t81, t101, -t133, t21 * qJD(2) + t178 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t5 * qJD(2) - t206 + t70; 0, 0, 0, 0, t185, t146 * t232, t121, t109, 0, 0, 0, -t49 * qJD(3) + t116, -t50 * qJD(3) - t115, -t16 * qJD(3) + t199, -t203, -t17 * qJD(3) + t200, -t13 * qJD(3) - t23 * qJD(4) - t209, -t7 * qJD(3) + t199, -t8 * qJD(3) + t200, t203 + t108, -t1 * qJD(3) - t20 * qJD(4) - t4 * qJD(5) - t213; 0, 0, 0, 0, 0, 0, t121, t109, 0, 0, 0, -pkin(2) * t193, -pkin(2) * t133, -t36 * qJD(3) + t120, 0, -t35 * qJD(3) + t130, -t158 * t188, t30 * qJD(3) + t120, t29 * qJD(3) + t130, t108, t6 * qJD(3) - t38 * qJD(5) + t132 * t86; 0, 0, 0, 0, 0, 0, t81, t241, t133, -t193, 0, -t156 - t184, -t155 + t186, -t160 - t184, t245, -t159 - t186, pkin(7) * t245 - t153, -t162 - t196, -t104 * qJD(3) + t163, t166, (-t104 * qJ(4) + t106 * t147) * qJD(3) + t106 * qJD(4) + t69 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t133, t101, -t152 + t184, t81, t101, -t133, -t154 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t161 + t70; 0, 0, 0, 0, 0, 0, -t81, -t241, 0, 0, 0, t49 * qJD(2) - t173, t50 * qJD(2) - t172, t16 * qJD(2) + t204, 0, t17 * qJD(2) + t205, t13 * qJD(2) - t177, t7 * qJD(2) + t131 - t207, t8 * qJD(2) - t191 - t208, 0, t1 * qJD(2) - t230 - t45; 0, 0, 0, 0, 0, 0, -t81, -t241, 0, 0, 0, t156, t155, t160, 0, t159, t153, t131 + t162, -t163 - t191, 0, -t164 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t140, 0, qJD(4), 0, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t139, 0, qJD(3), 0, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t103, 0, -t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, -t101, t23 * qJD(2) + t176, -t81, -t101, 0, t20 * qJD(2) - t131 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, 0, -t101, t152, -t81, -t101, 0, -t131 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t139, 0, -qJD(3), 0, -t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, t133, -t242, t4 * qJD(2) + t206 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, t133, -t242, -t161 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t103, 0, t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
