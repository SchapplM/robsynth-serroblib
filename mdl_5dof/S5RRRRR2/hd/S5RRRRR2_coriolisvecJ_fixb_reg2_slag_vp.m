% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:31
% EndTime: 2019-03-29 15:26:38
% DurationCPUTime: 2.84s
% Computational Cost: add. (3105->281), mult. (7466->423), div. (0->0), fcn. (5292->8), ass. (0->194)
t241 = cos(qJ(4));
t171 = t241 * qJD(4);
t251 = t241 * qJD(3) + t171;
t113 = qJD(1) + qJD(2);
t192 = qJD(3) + qJD(4);
t117 = sin(qJ(4));
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t90 = t117 * t121 + t241 * t118;
t249 = t192 * t90;
t51 = t249 * t113;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t194 = qJD(5) * t120;
t122 = cos(qJ(2));
t198 = qJD(2) * t122;
t170 = qJD(1) * t198;
t161 = pkin(1) * t170;
t142 = t121 * t161;
t119 = sin(qJ(2));
t205 = t118 * t119;
t189 = pkin(1) * t205;
t163 = qJD(3) * t189;
t143 = qJD(1) * t163;
t175 = t118 * t198;
t204 = t119 * t121;
t225 = pkin(1) * qJD(1);
t79 = (-qJD(3) * t204 - t175) * t225;
t226 = -t117 * t143 - t241 * t79;
t131 = t117 * t142 + t226;
t186 = t119 * t225;
t165 = t118 * t186;
t135 = -qJD(3) * pkin(2) + t165;
t164 = t121 * t186;
t138 = t241 * t164;
t71 = -t117 * t135 + t138;
t32 = t71 * qJD(4) + t131;
t70 = t117 * t164 + t241 * t135;
t250 = t32 * t116 + t70 * t194;
t195 = qJD(5) * t116;
t248 = t32 * t120 - t70 * t195;
t247 = t70 * qJD(4);
t114 = t118 ^ 2;
t115 = t121 ^ 2;
t200 = t114 + t115;
t246 = t251 * t121;
t199 = qJD(1) * t122;
t187 = pkin(1) * t199;
t207 = t113 * t121;
t88 = -pkin(2) * t207 - t187;
t144 = t116 * t71 - t120 * t88;
t245 = qJD(5) * t144;
t244 = t113 * t200;
t243 = t119 * pkin(1) ^ 2 * (-0.1e1 + t200);
t206 = t117 * t118;
t150 = t192 * t206;
t210 = t246 * t113;
t50 = t113 * t150 - t210;
t85 = t90 * t113;
t65 = t116 * t192 + t120 * t85;
t27 = qJD(5) * t65 - t116 * t50;
t31 = t117 * t79 + (t142 - t143) * t241 - t247;
t162 = qJD(2) * t186;
t197 = qJD(3) * t118;
t174 = t113 * t197;
t86 = pkin(2) * t174 + t162;
t15 = t116 * t86 + t120 * t31 - t245;
t242 = -0.2e1 * qJD(3);
t176 = t241 * t121;
t136 = t176 - t206;
t240 = t51 * t136;
t239 = t51 * t90;
t166 = t120 * t192;
t63 = t116 * t85 - t166;
t238 = t63 * t70;
t181 = t113 * t206;
t83 = -t113 * t176 + t181;
t77 = qJD(5) + t83;
t237 = t63 * t77;
t236 = t65 * t63;
t235 = t65 * t70;
t234 = t65 * t77;
t66 = t150 - t246;
t43 = t70 * t66;
t73 = t136 * t186;
t233 = t70 * t73;
t134 = pkin(1) * t90;
t74 = t134 * t199;
t232 = t70 * t74;
t231 = t77 * t85;
t230 = t85 * t83;
t229 = t88 * t85;
t228 = -t136 * t86 + t249 * t88;
t227 = -t88 * t66 + t86 * t90;
t224 = pkin(1) * qJD(2);
t26 = -qJD(5) * t166 + t120 * t50 + t195 * t85;
t223 = t116 * t26;
t221 = t116 * t51;
t220 = t116 * t63;
t219 = t116 * t65;
t218 = t116 * t83;
t217 = t120 * t27;
t216 = t120 * t51;
t215 = t120 * t63;
t214 = t120 * t65;
t213 = t120 * t83;
t212 = t121 * t51;
t14 = t15 * t120;
t46 = t116 * t88 + t120 * t71;
t209 = qJD(5) * t46;
t208 = t113 * t118;
t123 = qJD(3) ^ 2;
t203 = t123 * t118;
t111 = t123 * t121;
t202 = t200 * t161;
t201 = t114 - t115;
t196 = qJD(3) * t122;
t193 = -qJD(1) - t113;
t190 = t144 * t213 - t46 * t218 + t14;
t188 = pkin(2) * t208;
t185 = t241 * t32;
t184 = t90 * t195;
t183 = t90 * t194;
t182 = qJD(5) * t121 * t77;
t112 = t113 ^ 2;
t180 = t118 * t112 * t121;
t179 = t116 * t241;
t178 = t120 * t241;
t177 = (-t77 + t83) * t70;
t173 = t118 * t196;
t169 = t136 * t31 - t249 * t71 + t32 * t90 - t43;
t16 = -t116 * t31 + t120 * t86 - t209;
t168 = -t16 - t209;
t167 = t120 * t77;
t160 = t46 * t85 + t250;
t159 = t77 * t171;
t158 = t121 * t174;
t155 = -t116 * t43 - t16 * t136 - t144 * t249 + t250 * t90;
t154 = (-qJD(2) + t113) * t225;
t153 = t193 * t224;
t38 = -t117 * t163 - qJD(4) * t117 * t189 + (t241 * t175 + (t117 * t198 + t251 * t119) * t121) * pkin(1);
t81 = t119 * t134;
t152 = t32 * t81 + t38 * t70;
t151 = -t66 * t77 + t239;
t147 = t116 * t46 - t120 * t144;
t146 = -t116 * t144 - t120 * t46;
t145 = t214 + t220;
t105 = -pkin(1) * t122 - pkin(2) * t121;
t82 = t136 * t119 * pkin(1);
t55 = t105 * t120 - t116 * t82;
t56 = t105 * t116 + t120 * t82;
t141 = t119 * t154;
t140 = t122 * t154;
t139 = t144 * t85 - t248;
t137 = t197 * t77 - t212;
t133 = t136 * t122;
t132 = -t120 * t43 + t15 * t136 + t248 * t90 - t249 * t46;
t129 = -qJD(5) * t147 - t16 * t116;
t127 = -t144 * t184 + (-t144 * t66 + t168 * t90) * t120 + (-t15 * t90 + t46 * t66) * t116;
t126 = t88 * t83 - t31;
t103 = t118 * t162;
t96 = -0.2e1 * t158;
t95 = 0.2e1 * t158;
t94 = pkin(2) * t197 + t119 * t224;
t80 = t201 * t113 * t242;
t76 = t133 * t225;
t75 = t90 * t186;
t60 = t249 * t192;
t59 = t66 * t192;
t58 = t116 * t186 + t120 * t76;
t57 = -t116 * t76 + t120 * t186;
t54 = t116 * t188 - t120 * t75;
t53 = t116 * t75 + t120 * t188;
t40 = -t83 ^ 2 + t85 ^ 2;
t37 = (qJD(2) * t133 - t119 * t249) * pkin(1);
t36 = t192 * t85 - t51;
t35 = t210 + (-t181 + t83) * t192;
t21 = t249 * t83 - t240;
t20 = -t50 * t90 - t66 * t85;
t19 = t249 * t77 - t240;
t18 = -qJD(5) * t56 - t116 * t37 + t120 * t94;
t17 = qJD(5) * t55 + t116 * t94 + t120 * t37;
t11 = t167 * t77 - t65 * t85 + t221;
t10 = -t116 * t77 ^ 2 + t63 * t85 + t216;
t9 = t220 * t77 - t217;
t8 = t167 * t65 - t223;
t7 = t63 * t183 + (t27 * t90 - t63 * t66) * t116;
t6 = -t65 * t184 + (-t26 * t90 - t65 * t66) * t120;
t5 = -t136 * t50 - t249 * t85 + t66 * t83 - t239;
t4 = -t116 * t151 + t136 * t27 - t183 * t77 - t249 * t63;
t3 = t120 * t151 + t136 * t26 - t184 * t77 + t249 * t65;
t2 = (-t26 - t237) * t120 + (-t27 - t234) * t116;
t1 = (t215 + t219) * t66 + (t223 - t217 + (-t214 + t220) * qJD(5)) * t90;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t153, t122 * t153, 0, 0, t95, t80, t111, t96, -t203, 0 (-t119 * t111 + t193 * t173 + (t193 * t204 - t173) * qJD(2)) * pkin(1), t103 + ((qJD(2) * t113 + t123) * t205 - 0.2e1 * t196 * t207) * pkin(1), pkin(1) * t198 * t244 + t202, 0.2e1 * t170 * t243, t20, t5, -t59, t21, -t60, 0, t105 * t51 - t192 * t38 + t94 * t83 + t228, -t105 * t50 - t192 * t37 + t94 * t85 + t227, -t37 * t83 + t38 * t85 - t50 * t81 - t51 * t82 + t169, t105 * t86 + t31 * t82 + t37 * t71 + t88 * t94 + t152, t6, t1, t3, t7, t4, t19, t18 * t77 + t27 * t81 + t38 * t63 + t51 * t55 + t155, -t17 * t77 - t26 * t81 + t38 * t65 - t51 * t56 + t132, -t17 * t63 - t18 * t65 + t26 * t55 - t27 * t56 + t127, -t144 * t18 + t15 * t56 + t16 * t55 + t17 * t46 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t140, 0, 0, t95, t80, t111, t96, -t203, 0, t121 * t141, -t113 * t165 + t103, -t187 * t244 + t202, -qJD(1) ^ 2 * t122 * t243, t20, t5, -t59, t21, -t60, 0, -t83 * t186 + t74 * t192 + (t197 * t83 - t212) * pkin(2) + t228, -t85 * t186 + t76 * t192 + (t121 * t50 + t197 * t85) * pkin(2) + t227, -t74 * t85 + t76 * t83 + t169, -t88 * t186 - t232 - t71 * t76 + (-t121 * t86 + t197 * t88) * pkin(2), t6, t1, t3, t7, t4, t19, -t57 * t77 - t63 * t74 + (t116 * t182 + t120 * t137) * pkin(2) + t155, t58 * t77 - t65 * t74 + (-t116 * t137 + t120 * t182) * pkin(2) + t132, t57 * t65 + t58 * t63 + (-t145 * t197 + (t116 * t27 - t120 * t26 + (t215 - t219) * qJD(5)) * t121) * pkin(2) + t127, t144 * t57 - t46 * t58 - t232 + (t147 * t197 + (qJD(5) * t146 - t116 * t15 - t120 * t16) * t121) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t201 * t112, 0, t180, 0, 0, t118 * t140, t121 * t140, 0, 0, t230, t40, t35, -t230, t36, 0, -qJD(4) * t138 - t83 * t188 - t229 + t73 * t192 + (-t142 + (t165 + (t242 - qJD(4)) * pkin(2)) * qJD(4)) * t117 - t226, -t188 * t85 + t126 + (-pkin(2) * t171 - t75) * t192 (t71 - t73) * t85 + (t70 - t75) * t83 + (t241 * t50 - t117 * t51 + (t117 * t85 - t241 * t83) * qJD(4)) * pkin(2), -t233 + t71 * t75 + (-t88 * t208 - t185 + t117 * t31 + (t117 * t70 + t241 * t71) * qJD(4)) * pkin(2), t8, t2, t11, t9, t10, -t231, t70 * t218 - t53 * t77 - t73 * t63 + (-t116 * t159 - t241 * t27 + (qJD(4) * t63 - t194 * t77 - t221) * t117) * pkin(2) + t139, t70 * t213 + t54 * t77 - t73 * t65 + (-t120 * t159 + t241 * t26 + (qJD(4) * t65 + t195 * t77 - t216) * t117) * pkin(2) + t160, t53 * t65 + t54 * t63 + ((-t178 * t63 + t179 * t65) * qJD(4) + (qJD(5) * t145 - t217 - t223) * t117) * pkin(2) + t129 + t190, t144 * t53 - t46 * t54 - t233 + (-t185 + (t144 * t179 + t178 * t46) * qJD(4) + (t129 + t14 + t247) * t117) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t40, t35, -t230, t36, 0, t71 * qJD(3) - t131 - t229, -t192 * t70 + t126, 0, 0, t8, t2, t11, t9, t10, -t231, t116 * t177 - t63 * t71 + t139, t120 * t177 - t65 * t71 + t160 (-t238 + t245) * t120 + (t168 + t235) * t116 + t190 (-t146 - t71) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t63 ^ 2 + t65 ^ 2, -t26 + t237, -t236, t234 - t27, t51, t46 * t77 + t16 - t235, -t144 * t77 - t15 + t238, 0, 0;];
tauc_reg  = t12;
