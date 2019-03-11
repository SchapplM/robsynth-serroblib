% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:19
% EndTime: 2019-03-09 03:35:26
% DurationCPUTime: 2.32s
% Computational Cost: add. (4074->262), mult. (9963->362), div. (0->0), fcn. (7632->10), ass. (0->165)
t150 = cos(qJ(6));
t192 = qJD(6) * t150;
t143 = sin(pkin(11));
t149 = sin(qJ(3));
t145 = cos(pkin(11));
t152 = cos(qJ(3));
t200 = t145 * t152;
t126 = -t143 * t149 + t200;
t119 = t126 * qJD(1);
t151 = cos(qJ(5));
t107 = t151 * t119;
t127 = t143 * t152 + t145 * t149;
t121 = t127 * qJD(1);
t148 = sin(qJ(5));
t71 = -t148 * t121 + t107;
t241 = t150 * t71;
t248 = t192 - t241;
t140 = qJD(3) + qJD(5);
t204 = t71 * t140;
t120 = t127 * qJD(3);
t110 = qJD(1) * t120;
t189 = qJD(1) * qJD(3);
t184 = t149 * t189;
t111 = -t143 * t184 + t189 * t200;
t194 = qJD(5) * t148;
t37 = qJD(5) * t107 - t148 * t110 + t151 * t111 - t121 * t194;
t247 = t37 - t204;
t196 = -qJD(6) + t71;
t246 = qJD(6) + t196;
t147 = sin(qJ(6));
t169 = t148 * t119 + t151 * t121;
t193 = qJD(6) * t147;
t23 = t140 * t192 + t150 * t37 - t169 * t193;
t65 = t147 * t140 + t150 * t169;
t24 = qJD(6) * t65 + t147 * t37;
t63 = -t150 * t140 + t147 * t169;
t245 = -t147 * t24 + t23 * t150 - t248 * t63;
t21 = t23 * t147;
t244 = t248 * t65 + t21;
t38 = qJD(5) * t169 + t151 * t110 + t148 * t111;
t34 = t147 * t38;
t66 = t196 * t192;
t214 = t34 - t66;
t219 = t65 * t169;
t243 = t196 * t241 + t214 - t219;
t222 = t121 * pkin(8);
t134 = sin(pkin(10)) * pkin(1) + pkin(7);
t197 = qJ(4) + t134;
t174 = t197 * qJD(1);
t96 = t149 * qJD(2) + t152 * t174;
t86 = t143 * t96;
t213 = qJD(3) * pkin(3);
t95 = t152 * qJD(2) - t174 * t149;
t90 = t95 + t213;
t53 = t145 * t90 - t86;
t42 = qJD(3) * pkin(4) - t222 + t53;
t223 = t119 * pkin(8);
t212 = t145 * t96;
t54 = t143 * t90 + t212;
t44 = t54 + t223;
t16 = -t148 * t44 + t151 * t42;
t14 = -t140 * pkin(5) - t16;
t242 = t14 * t71;
t240 = t169 * t71;
t239 = t147 * t196;
t205 = t169 * t140;
t238 = -t38 + t205;
t236 = t169 ^ 2 - t71 ^ 2;
t40 = pkin(5) * t169 - t71 * pkin(9);
t188 = qJD(1) * qJD(4);
t78 = t95 * qJD(3) + t152 * t188;
t79 = -t96 * qJD(3) - t149 * t188;
t45 = -t143 * t78 + t145 * t79;
t31 = -t111 * pkin(8) + t45;
t46 = t143 * t79 + t145 * t78;
t32 = -t110 * pkin(8) + t46;
t2 = (qJD(5) * t42 + t32) * t151 + t148 * t31 - t44 * t194;
t136 = -cos(pkin(10)) * pkin(1) - pkin(2);
t166 = -t152 * pkin(3) + t136;
t160 = t166 * qJD(1);
t117 = qJD(4) + t160;
t80 = -t119 * pkin(4) + t117;
t235 = -t80 * t71 - t2;
t218 = t169 * t63;
t232 = t196 * t169;
t36 = t150 * t38;
t231 = -t193 * t196 - t36;
t17 = t148 * t42 + t151 * t44;
t15 = t140 * pkin(9) + t17;
t26 = -pkin(5) * t71 - pkin(9) * t169 + t80;
t170 = t147 * t15 - t150 * t26;
t230 = t14 * t193 + t169 * t170;
t3 = qJD(5) * t17 + t148 * t32 - t151 * t31;
t5 = t147 * t26 + t150 * t15;
t229 = t14 * t192 + t3 * t147 + t5 * t169;
t228 = -t169 * t80 - t3;
t82 = t148 * t126 + t151 * t127;
t217 = t82 * t38;
t123 = t126 * qJD(3);
t168 = t151 * t126 - t148 * t127;
t47 = qJD(5) * t168 - t148 * t120 + t151 * t123;
t171 = -t196 * t47 + t217;
t185 = t82 * t193;
t227 = -t150 * t171 - t185 * t196;
t124 = t197 * t149;
t125 = t197 * t152;
t76 = -t145 * t124 - t143 * t125;
t61 = -t127 * pkin(8) + t76;
t77 = -t143 * t124 + t145 * t125;
t62 = t126 * pkin(8) + t77;
t28 = t148 * t61 + t151 * t62;
t94 = -t126 * pkin(4) + t166;
t33 = -pkin(5) * t168 - t82 * pkin(9) + t94;
t27 = t148 * t62 - t151 * t61;
t175 = qJD(3) * t197;
t100 = -t149 * qJD(4) - t152 * t175;
t99 = t152 * qJD(4) - t149 * t175;
t59 = t145 * t100 - t143 * t99;
t51 = -t123 * pkin(8) + t59;
t60 = t143 * t100 + t145 * t99;
t52 = -t120 * pkin(8) + t60;
t6 = -qJD(5) * t27 + t148 * t51 + t151 * t52;
t226 = t14 * t47 + (qJD(6) * t33 + t6) * t196 + (qJD(6) * t26 + t2) * t168 - t28 * t38 + t3 * t82;
t224 = pkin(3) * t143;
t221 = t14 * t82;
t220 = t33 * t38;
t48 = qJD(5) * t82 + t151 * t120 + t148 * t123;
t216 = -t168 * t23 + t65 * t48;
t56 = t145 * t95 - t86;
t210 = t147 * t65;
t206 = t47 * t140;
t135 = t145 * pkin(3) + pkin(4);
t162 = t151 * t135 - t148 * t224;
t55 = -t143 * t95 - t212;
t49 = t55 - t223;
t50 = t56 - t222;
t203 = -t162 * qJD(5) + t148 * t49 + t151 * t50;
t163 = t148 * t135 + t151 * t224;
t202 = t163 * qJD(5) - t148 * t50 + t151 * t49;
t153 = qJD(3) ^ 2;
t199 = t153 * t149;
t198 = t153 * t152;
t195 = t149 ^ 2 - t152 ^ 2;
t129 = qJD(1) * t136;
t191 = t149 * qJD(1);
t138 = t149 * t213;
t132 = pkin(3) * t184;
t85 = t110 * pkin(4) + t132;
t98 = t120 * pkin(4) + t138;
t97 = pkin(3) * t191 + t121 * pkin(4);
t116 = pkin(9) + t163;
t177 = qJD(6) * t116 + t40 + t97;
t172 = t168 * t24 - t48 * t63;
t165 = -t239 * t71 - t231;
t164 = 0.2e1 * qJD(3) * t129;
t159 = -t116 * t38 - t196 * t203 - t242;
t155 = -t171 * t147 + t82 * t66;
t154 = qJD(1) ^ 2;
t115 = -pkin(5) - t162;
t43 = t48 * t140;
t10 = t48 * pkin(5) - t47 * pkin(9) + t98;
t9 = t38 * pkin(5) - t37 * pkin(9) + t85;
t8 = t150 * t9;
t7 = qJD(5) * t28 + t148 * t52 - t151 * t51;
t1 = [0, 0, 0, 0, 0.2e1 * t152 * t184, -0.2e1 * t195 * t189, t198, -t199, 0, -t134 * t198 + t149 * t164, t134 * t199 + t152 * t164, -t77 * t110 - t76 * t111 + t60 * t119 - t54 * t120 - t59 * t121 - t53 * t123 + t46 * t126 - t45 * t127, t45 * t76 + t46 * t77 + t53 * t59 + t54 * t60 + (t117 + t160) * t138, t169 * t47 + t37 * t82, t168 * t37 - t169 * t48 + t47 * t71 - t217, t206, -t43, 0, -t7 * t140 - t168 * t85 + t94 * t38 + t80 * t48 - t71 * t98, -t6 * t140 + t169 * t98 + t94 * t37 + t80 * t47 + t85 * t82, -t65 * t185 + (t23 * t82 + t47 * t65) * t150 (-t150 * t63 - t210) * t47 + (-t21 - t150 * t24 + (t147 * t63 - t150 * t65) * qJD(6)) * t82, t216 - t227, t155 + t172, -t168 * t38 - t196 * t48, t27 * t24 - t170 * t48 + t7 * t63 - t8 * t168 + (-t10 * t196 + t220 + (t15 * t168 + t196 * t28 + t221) * qJD(6)) * t150 + t226 * t147, t27 * t23 - t5 * t48 + t7 * t65 + ((-qJD(6) * t28 + t10) * t196 - t220 + (-qJD(6) * t15 + t9) * t168 - qJD(6) * t221) * t147 + t226 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, -t198, -t127 * t110 - t126 * t111 + t123 * t119 + t120 * t121, -t53 * t120 + t54 * t123 + t45 * t126 + t46 * t127, 0, 0, 0, 0, 0, -t43, -t206, 0, 0, 0, 0, 0, t155 - t172, t216 + t227; 0, 0, 0, 0, -t149 * t154 * t152, t195 * t154, 0, 0, 0, -t129 * t191, -t129 * t152 * qJD(1) (t54 + t55) * t121 + (t53 - t56) * t119 + (-t110 * t143 - t111 * t145) * pkin(3), -t53 * t55 - t54 * t56 + (-t117 * t191 + t143 * t46 + t145 * t45) * pkin(3), -t240, t236, t247, t238, 0, -t140 * t202 + t71 * t97 + t228, t140 * t203 - t169 * t97 + t235, t244, t196 * t210 + t245, t243, t165 + t218, t232, t115 * t24 + t202 * t63 + (t177 * t196 - t3) * t150 + t159 * t147 + t230, t115 * t23 + t150 * t159 - t177 * t239 + t202 * t65 + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 ^ 2 - t121 ^ 2, -t54 * t119 + t53 * t121 + t132, 0, 0, 0, 0, 0, t38 + t205, t37 + t204, 0, 0, 0, 0, 0, t165 - t218, -t150 * t196 ^ 2 - t219 - t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t236, t247, t238, 0, t17 * t140 + t228, t16 * t140 + t235, t244, t239 * t65 + t245, t243, -t196 * t239 + t218 + t36, t232, -pkin(5) * t24 - t3 * t150 + (-t147 * t16 + t150 * t40) * t196 - t17 * t63 - t147 * t242 - t214 * pkin(9) + t230, -pkin(5) * t23 - (t147 * t40 + t150 * t16) * t196 - t17 * t65 - t14 * t241 + t231 * pkin(9) + t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t63 ^ 2 + t65 ^ 2, -t196 * t63 + t23, -t196 * t65 - t24, t38, -t14 * t65 - t147 * t2 - t246 * t5 + t8, t14 * t63 - t147 * t9 - t150 * t2 + t170 * t246;];
tauc_reg  = t1;
