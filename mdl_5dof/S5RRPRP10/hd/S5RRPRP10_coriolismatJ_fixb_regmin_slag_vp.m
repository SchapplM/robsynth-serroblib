% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:08
% EndTime: 2019-12-31 20:11:13
% DurationCPUTime: 1.88s
% Computational Cost: add. (1949->226), mult. (3512->339), div. (0->0), fcn. (2938->4), ass. (0->192)
t139 = sin(qJ(4));
t247 = t139 / 0.2e1;
t249 = pkin(3) + pkin(6);
t140 = sin(qJ(2));
t244 = t140 * pkin(4);
t142 = cos(qJ(2));
t223 = t139 * t142;
t143 = -pkin(2) - pkin(7);
t222 = t140 * qJ(3);
t159 = -t142 * t143 + t222;
t78 = -pkin(1) - t159;
t101 = t249 * t140;
t141 = cos(qJ(4));
t90 = t141 * t101;
t54 = t139 * t78 - t90;
t48 = qJ(5) * t223 - t54;
t45 = t48 + t244;
t256 = -t45 + t48;
t235 = t45 * t139;
t221 = t141 * t142;
t225 = t139 * t101;
t55 = t141 * t78 + t225;
t49 = -qJ(5) * t221 + t55;
t17 = (-t49 * t141 + t235) * t140;
t255 = -qJ(5) + t143;
t253 = t48 * t247 - t235 / 0.2e1;
t135 = t139 ^ 2;
t137 = t141 ^ 2;
t110 = t135 - t137;
t191 = t142 * qJD(1);
t175 = t141 * t191;
t147 = t110 * qJD(2) + 0.2e1 * t139 * t175;
t190 = t142 * qJD(3);
t252 = qJD(2) * t159 - t190;
t227 = qJ(5) * t140;
t242 = t142 * pkin(4);
t220 = t142 * qJ(3);
t100 = t140 * pkin(2) - t220;
t82 = t140 * pkin(7) + t100;
t133 = t142 * pkin(6);
t134 = t142 * pkin(3);
t102 = t133 + t134;
t92 = t102 * t141;
t47 = t242 + t92 + (-t82 - t227) * t139;
t251 = t47 / 0.2e1;
t93 = t255 * t139;
t250 = -t93 / 0.2e1;
t248 = -t135 / 0.2e1;
t246 = -t140 / 0.2e1;
t245 = -t141 / 0.2e1;
t243 = t141 * pkin(4);
t71 = t141 * t82;
t91 = t102 * t139;
t241 = t71 + t91;
t240 = pkin(4) * qJD(4);
t239 = t139 * t82;
t238 = t140 * t93;
t94 = t255 * t141;
t237 = t140 * t94;
t183 = -t244 / 0.2e1;
t106 = t139 * t183;
t3 = t106 + t253;
t236 = t3 * qJD(1);
t50 = t141 * t227 + t241;
t73 = (-t243 - t249) * t140;
t74 = pkin(4) * t221 + t102;
t5 = t45 * t47 + t49 * t50 + t74 * t73;
t232 = t5 * qJD(1);
t6 = t50 * t221 - t47 * t223 + t17;
t231 = t6 * qJD(1);
t188 = pkin(4) * t223;
t7 = -t74 * t188 + t256 * t49;
t230 = t7 * qJD(1);
t8 = t256 * t221;
t229 = t8 * qJD(1);
t224 = t139 * t140;
t105 = pkin(4) * t224 / 0.2e1;
t9 = t105 - t253;
t228 = t9 * qJD(1);
t226 = qJD(1) * t17;
t15 = (-t54 - t90) * t142 - t140 * t239;
t219 = t15 * qJD(1);
t16 = t241 * t140 - t102 * t224 + (t55 - t225) * t142;
t218 = t16 * qJD(1);
t18 = t49 * t221 - t45 * t223;
t217 = t18 * qJD(1);
t24 = t102 * t221 - t54 * t140;
t216 = t24 * qJD(1);
t25 = -t102 * t223 - t55 * t140;
t215 = t25 * qJD(1);
t162 = -t142 * pkin(2) - t222;
t95 = -pkin(1) + t162;
t59 = t100 * t140 + t95 * t142;
t214 = t59 * qJD(1);
t60 = t100 * t142 - t95 * t140;
t213 = t60 * qJD(1);
t109 = t135 + t137;
t70 = t109 * t142 * t140;
t212 = t70 * qJD(1);
t136 = t140 ^ 2;
t138 = t142 ^ 2;
t111 = t138 - t136;
t86 = t111 * t139;
t211 = t86 * qJD(1);
t87 = t109 * t138;
t210 = t87 * qJD(1);
t89 = t111 * t141;
t209 = t89 * qJD(1);
t168 = t248 - t137 / 0.2e1;
t97 = -0.1e1 / 0.2e1 + t168;
t208 = t97 * qJD(2);
t207 = qJD(1) * t140;
t206 = qJD(1) * t141;
t205 = qJD(2) * qJ(3);
t204 = qJD(3) * t140;
t203 = qJD(3) * t141;
t202 = qJD(4) * t139;
t201 = qJD(4) * t140;
t200 = qJD(4) * t141;
t199 = qJD(4) * t143;
t198 = t109 * qJD(2);
t197 = t111 * qJD(1);
t196 = t136 * qJD(1);
t195 = t136 * qJD(3);
t194 = t139 * qJD(2);
t193 = t140 * qJD(2);
t192 = t141 * qJD(2);
t130 = t142 * qJD(2);
t189 = t142 * qJD(4);
t187 = pkin(1) * t207;
t186 = pkin(1) * t191;
t185 = pkin(4) * t202;
t184 = pkin(6) * t193;
t182 = t243 / 0.2e1;
t181 = t242 / 0.2e1;
t121 = t139 * pkin(4) + qJ(3);
t180 = t121 * t223;
t179 = t95 * t100 * qJD(1);
t178 = t95 * t207;
t174 = t139 * t189;
t173 = t141 * t189;
t172 = t139 * t130;
t116 = t140 * t130;
t115 = t140 * t191;
t171 = t139 * t200;
t170 = t139 * t192;
t169 = t133 / 0.2e1 + t134 / 0.2e1;
t165 = qJD(4) + t207;
t163 = t142 * t170;
t1 = (-t48 / 0.2e1 + t45 / 0.2e1) * t93 + (t251 + t180 / 0.2e1 + t74 * t245) * pkin(4);
t19 = t121 * t243;
t160 = -t1 * qJD(1) + t19 * qJD(2);
t144 = (t142 * t250 - t45 / 0.2e1) * t141 + (-t49 / 0.2e1 + t94 * t142 / 0.2e1) * t139;
t11 = (pkin(6) / 0.2e1 + pkin(3) / 0.2e1 + t182) * t140 + t144;
t57 = t93 * t139 + t94 * t141;
t158 = t11 * qJD(1) - t57 * qJD(2);
t14 = (t237 / 0.2e1 - t50 / 0.2e1) * t139 + (-t238 / 0.2e1 + t181 - t47 / 0.2e1) * t141 + t169;
t157 = t14 * qJD(1) + t121 * qJD(2);
t156 = -t196 - t201;
t155 = t143 * t246 - t220 / 0.2e1;
t146 = t155 * t141;
t52 = t71 / 0.2e1 + t146;
t154 = qJ(3) * t194 - t52 * qJD(1);
t51 = (t82 / 0.2e1 + t155) * t139;
t153 = -qJ(3) * t192 - t51 * qJD(1);
t81 = (t137 / 0.2e1 + t248) * t142;
t152 = t81 * qJD(1) + t170;
t151 = t165 * t223;
t150 = t139 * t191 - t192;
t149 = t139 * t138 * t206 - t81 * qJD(2);
t88 = t110 * t138;
t148 = -t88 * qJD(1) + 0.2e1 * t163;
t145 = t162 * qJD(2) + t190;
t127 = pkin(6) * t130;
t122 = t130 / 0.2e1;
t118 = t141 * t130;
t117 = t140 * t206;
t114 = t139 * t207;
t96 = 0.1e1 / 0.2e1 + t168;
t85 = -t117 - t200;
t84 = -t114 - t202;
t83 = t115 + t189 / 0.2e1;
t79 = t150 * pkin(4);
t72 = t81 * qJD(4);
t29 = -t91 - t71 / 0.2e1 + t146;
t28 = t92 - t239 / 0.2e1 + t155 * t139;
t13 = t50 * t247 + (t93 * t245 + t94 * t247) * t140 + t169 + (t181 + t251) * t141;
t12 = t141 * t183 + t249 * t246 + t144;
t10 = t105 + t253;
t4 = t106 - t253;
t2 = t48 * t93 / 0.2e1 + t45 * t250 + t74 * t182 + (-t180 / 0.2e1 + t251) * pkin(4);
t20 = [0, 0, 0, t116, t111 * qJD(2), 0, 0, 0, -pkin(1) * t193, -pkin(1) * t130, 0, t60 * qJD(2) - t140 * t190, -t59 * qJD(2) + t195, (qJD(2) * t100 - t204) * t95, -t116 * t135 + t138 * t171, -t88 * qJD(4) - 0.2e1 * t140 * t163, -t86 * qJD(2) - t140 * t173, -t89 * qJD(2) + t140 * t174, t116, t15 * qJD(2) + t25 * qJD(4) + t139 * t195, -t16 * qJD(2) - t24 * qJD(4) + t141 * t195, -t6 * qJD(2) + t70 * qJD(3) - t8 * qJD(4) + t87 * qJD(5), qJD(2) * t5 + qJD(3) * t17 + qJD(4) * t7 - qJD(5) * t18; 0, 0, 0, t115, t197, t130, -t193, 0, -t127 - t187, t184 - t186, t145, t127 + t213, -t184 - t214, t145 * pkin(6) + t179, -t72 + (-t135 * t191 + t170) * t140, -t147 * t140 + 0.2e1 * t142 * t171, t118 - t211, -t172 - t209, t83, t28 * qJD(4) - t101 * t194 - t252 * t141 + t219, t29 * qJD(4) - t101 * t192 + t252 * t139 - t218, -t231 + ((-t47 + t238) * t141 + (-t50 - t237) * t139) * qJD(2) + t4 * qJD(4), t232 + (t73 * t121 + t47 * t94 + t50 * t93) * qJD(2) + t13 * qJD(3) + t2 * qJD(4) + t12 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t115, t196, t127 - t178, 0, 0, 0, 0, 0, t139 * t196 + t118, t141 * t196 - t172, t212, qJD(2) * t13 + qJD(4) * t10 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t148, -t165 * t221, t151, t122, t28 * qJD(2) - t55 * qJD(4) + t215, t29 * qJD(2) + t54 * qJD(4) - t216, pkin(4) * t173 + t4 * qJD(2) - t229, t2 * qJD(2) + t10 * qJD(3) - t49 * t240 + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, qJD(2) * t12 - t217; 0, 0, 0, -t115, -t197, 0, 0, 0, t187, t186, 0, -t213, t214, -t179, t115 * t135 - t72, 0.2e1 * t141 * t151, -t139 * t201 + t211, -t140 * t200 + t209, -t83, qJD(4) * t51 - t219, qJD(4) * t52 + t218, -qJD(4) * t3 + t231, qJD(3) * t14 - qJD(4) * t1 + qJD(5) * t11 - t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t171, t110 * qJD(4), 0, 0, 0, qJ(3) * t200 + qJD(3) * t139, -qJ(3) * t202 + t203, t109 * qJD(5), t121 * qJD(3) + t19 * qJD(4) - t57 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t205, 0, 0, 0, 0, 0, t194, t192, 0, t96 * qJD(5) + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t147, t84, t85, -t191 / 0.2e1, -t139 * t199 - t153, -t141 * t199 - t154, t185 - t236, -t93 * t240 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t96 * qJD(3) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t196, t178, 0, 0, 0, 0, 0, t156 * t139, t156 * t141, -t212, -qJD(2) * t14 - qJD(4) * t9 - t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t205, 0, 0, 0, 0, 0, -t194, -t192, 0, t97 * qJD(5) - t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t85, 0, -t185 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t148, (t175 + t194) * t140, -t150 * t140, t122, -t51 * qJD(2) + t139 * t204 - t215, -t52 * qJD(2) + t140 * t203 + t216, qJD(2) * t3 + t229, t1 * qJD(2) + t9 * qJD(3) + qJD(5) * t188 - t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t147, t114, t117, t191 / 0.2e1, t153, t154, t236, -qJD(5) * t243 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t117, 0, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, -pkin(4) * t174 - t11 * qJD(2) + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, pkin(4) * t200 - t97 * qJD(3) - t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t20;
