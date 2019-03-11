% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:01
% EndTime: 2019-03-09 05:21:07
% DurationCPUTime: 2.33s
% Computational Cost: add. (4853->289), mult. (10700->398), div. (0->0), fcn. (7651->8), ass. (0->181)
t142 = sin(qJ(4));
t145 = cos(qJ(4));
t146 = cos(qJ(3));
t202 = qJD(1) * t146;
t143 = sin(qJ(3));
t203 = qJD(1) * t143;
t103 = t142 * t203 - t145 * t202;
t109 = t142 * t146 + t145 * t143;
t104 = t109 * qJD(1);
t140 = sin(pkin(10));
t216 = cos(pkin(10));
t179 = t140 * t103 - t216 * t104;
t238 = qJD(6) - t179;
t251 = qJD(6) - t238;
t147 = -pkin(1) - pkin(7);
t118 = t147 * qJD(1) + qJD(2);
t183 = pkin(8) * qJD(1) - t118;
t201 = qJD(3) * t143;
t85 = t183 * t201;
t200 = qJD(3) * t146;
t86 = t183 * t200;
t168 = t142 * t86 + t145 * t85;
t93 = -pkin(8) * t203 + t143 * t118;
t82 = t145 * t93;
t94 = -pkin(8) * t202 + t146 * t118;
t84 = qJD(3) * pkin(3) + t94;
t169 = -t142 * t84 - t82;
t154 = t169 * qJD(4) + t168;
t135 = qJD(3) + qJD(4);
t76 = t135 * t109;
t68 = t76 * qJD(1);
t151 = t68 * qJ(5) + t103 * qJD(5) + t154;
t210 = t145 * t146;
t166 = t142 * t143 - t210;
t239 = qJD(1) * t166;
t152 = t135 * t239;
t199 = qJD(4) * t142;
t186 = t142 * t85 - t93 * t199;
t241 = (qJD(4) * t84 - t86) * t145;
t16 = qJ(5) * t152 - t104 * qJD(5) + t186 + t241;
t3 = t140 * t16 - t216 * t151;
t72 = t140 * t109 + t166 * t216;
t250 = t3 * t72;
t141 = sin(qJ(6));
t144 = cos(qJ(6));
t161 = -t216 * t103 - t140 * t104;
t56 = -t144 * t135 + t141 * t161;
t249 = t238 * t56;
t178 = t144 * t238;
t40 = -t140 * t68 - t216 * t152;
t221 = t141 * t40;
t248 = -t238 * t178 - t221;
t160 = t216 * t109 - t140 * t166;
t77 = t135 * t210 - t142 * t201 - t143 * t199;
t162 = -t140 * t76 + t216 * t77;
t214 = t104 * qJ(5);
t52 = -t169 - t214;
t222 = t140 * t52;
t81 = t142 * t93;
t187 = t145 * t84 - t81;
t97 = t103 * qJ(5);
t51 = t187 + t97;
t47 = t135 * pkin(4) + t51;
t21 = t216 * t47 - t222;
t48 = t216 * t52;
t22 = t140 * t47 + t48;
t4 = t140 * t151 + t216 * t16;
t44 = t140 * t77 + t216 * t76;
t247 = t160 * t4 + t162 * t22 - t21 * t44;
t197 = qJD(6) * t144;
t198 = qJD(6) * t141;
t41 = t140 * t152 - t216 * t68;
t23 = t135 * t197 + t144 * t41 - t161 * t198;
t58 = t141 * t135 + t144 * t161;
t246 = -t72 * t23 - t44 * t58;
t245 = t160 * t40;
t132 = pkin(3) * t202;
t234 = t103 * pkin(4);
t33 = pkin(5) * t161 - pkin(9) * t179 - t234;
t129 = t145 * pkin(3) + pkin(4);
t184 = t216 * t142;
t96 = pkin(3) * t184 + t140 * t129;
t90 = pkin(9) + t96;
t244 = (qJD(6) * t90 + t132 + t33) * t238;
t37 = t144 * t40;
t242 = t198 * t238 - t37;
t227 = pkin(8) - t147;
t107 = t227 * t201;
t114 = t227 * t146;
t108 = qJD(3) * t114;
t113 = t227 * t143;
t167 = t145 * t113 + t142 * t114;
t155 = t167 * qJD(4) + t145 * t107 + t142 * t108;
t153 = t76 * qJ(5) + qJD(5) * t166 + t155;
t211 = t145 * t114;
t159 = -qJD(4) * t211 + t142 * t107 - t145 * t108 + t113 * t199;
t30 = -t77 * qJ(5) - t109 * qJD(5) + t159;
t11 = t140 * t153 + t216 * t30;
t157 = qJ(5) * t166 + t142 * t113 - t211;
t61 = -t109 * qJ(5) - t167;
t35 = t140 * t157 + t216 * t61;
t173 = -t35 * t40 - t250;
t19 = -t135 * pkin(5) - t21;
t115 = pkin(3) * t203 + qJD(1) * qJ(2);
t78 = t104 * pkin(4) + qJD(5) + t115;
t31 = -pkin(5) * t179 - pkin(9) * t161 + t78;
t128 = t143 * pkin(3) + qJ(2);
t175 = t109 * pkin(4) + t128;
t36 = pkin(5) * t160 + pkin(9) * t72 + t175;
t236 = -t160 * (qJD(6) * t31 + t4) - (qJD(6) * t36 + t11) * t238 - t19 * t44 + t173;
t136 = qJD(1) * qJD(2);
t235 = 0.2e1 * t136;
t233 = t19 * t179;
t232 = t19 * t72;
t231 = t36 * t40;
t230 = t58 * t161;
t229 = t238 * t161;
t228 = t161 * t56;
t226 = t145 * t94 - t81;
t185 = -t142 * t94 - t82;
t164 = t185 + t214;
t223 = pkin(3) * qJD(4);
t53 = t97 + t226;
t225 = -t140 * t53 + t216 * t164 + (t140 * t145 + t184) * t223;
t212 = t140 * t142;
t224 = -t140 * t164 - t216 * t53 + (t216 * t145 - t212) * t223;
t220 = t141 * t41;
t219 = t141 * t179;
t218 = t23 * t141;
t217 = t77 * t135;
t215 = t103 * t104;
t213 = t115 * t103;
t148 = qJD(3) ^ 2;
t209 = t148 * t143;
t208 = t148 * t146;
t149 = qJD(1) ^ 2;
t207 = t149 * qJ(2);
t196 = qJD(1) * qJD(3);
t189 = t146 * t196;
t112 = pkin(3) * t189 + t136;
t205 = t143 ^ 2 - t146 ^ 2;
t204 = -t148 - t149;
t119 = pkin(3) * t200 + qJD(2);
t195 = 0.2e1 * qJD(1);
t193 = t72 * t198;
t192 = t238 * t197;
t20 = t135 * pkin(9) + t22;
t170 = t141 * t20 - t144 * t31;
t191 = t161 * t170 + t19 * t198;
t190 = -pkin(3) * t135 - t84;
t177 = qJD(6) * t160 + qJD(1);
t9 = t141 * t31 + t144 * t20;
t176 = t3 * t141 + t161 * t9 + t19 * t197;
t174 = t77 * pkin(4) + t119;
t172 = t161 * t22 + t179 * t21;
t171 = -t238 * t44 - t40 * t72;
t165 = t219 * t238 - t242;
t163 = t115 * t104 - t186;
t95 = -pkin(3) * t212 + t216 * t129;
t156 = -t224 * t238 - t90 * t40 - t233;
t150 = -pkin(4) * t152 + t112;
t127 = -t216 * pkin(4) - pkin(5);
t126 = t140 * pkin(4) + pkin(9);
t89 = -pkin(5) - t95;
t69 = t76 * t135;
t59 = t103 ^ 2 - t104 ^ 2;
t55 = -t103 * t135 + t152;
t54 = t104 * t135 - t68;
t34 = t140 * t61 - t216 * t157;
t26 = t216 * t51 - t222;
t25 = t140 * t51 + t48;
t24 = t58 * qJD(6) + t220;
t14 = pkin(5) * t162 + pkin(9) * t44 + t174;
t13 = t40 * pkin(5) - t41 * pkin(9) + t150;
t12 = t144 * t13;
t10 = t140 * t30 - t216 * t153;
t7 = t58 * t178 + t218;
t6 = -t230 - t248;
t5 = t165 + t228;
t1 = (t23 - t249) * t144 + (-t238 * t58 - t24) * t141;
t2 = [0, 0, 0, 0, t235, qJ(2) * t235, -0.2e1 * t143 * t189, 0.2e1 * t205 * t196, -t209, -t208, 0, -t147 * t209 + (qJ(2) * t200 + qJD(2) * t143) * t195, -t147 * t208 + (-qJ(2) * t201 + qJD(2) * t146) * t195, t103 * t76 + t166 * t68, t103 * t77 + t76 * t104 + t68 * t109 - t152 * t166, -t69, -t217, 0, t119 * t104 + t112 * t109 + t115 * t77 + (-t128 * t239 + t155) * t135, -t119 * t103 - t112 * t166 - t115 * t76 - t128 * t68 - t159 * t135, t10 * t161 + t11 * t179 + t34 * t41 + t173 - t247, -t21 * t10 + t22 * t11 + t150 * t175 + t78 * t174 + t3 * t34 + t4 * t35, t246 * t144 + t58 * t193 -(-t141 * t58 - t144 * t56) * t44 - (-t218 - t144 * t24 + (t141 * t56 - t144 * t58) * qJD(6)) * t72, t144 * t171 + t160 * t23 + t162 * t58 + t193 * t238, -t141 * t171 - t160 * t24 - t162 * t56 + t192 * t72, t162 * t238 + t245, t10 * t56 + t12 * t160 + t34 * t24 - t170 * t162 + (t14 * t238 + t231 + (-t160 * t20 - t238 * t35 - t232) * qJD(6)) * t144 + t236 * t141, t10 * t58 + t34 * t23 - t9 * t162 + (-(-qJD(6) * t35 + t14) * t238 - t231 - (-qJD(6) * t20 + t13) * t160 + qJD(6) * t232) * t141 + t236 * t144; 0, 0, 0, 0, -t149, -t207, 0, 0, 0, 0, 0, t204 * t143, t204 * t146, 0, 0, 0, 0, 0, -qJD(1) * t104 - t69, qJD(1) * t103 - t217, t161 * t44 + t162 * t179 + t72 * t41 - t245, -t78 * qJD(1) + t247 + t250, 0, 0, 0, 0, 0, -t160 * t221 + t72 * t24 + t44 * t56 + (-t141 * t162 - t144 * t177) * t238, -t160 * t37 + (t141 * t177 - t144 * t162) * t238 - t246; 0, 0, 0, 0, 0, 0, t146 * t149 * t143, -t205 * t149, 0, 0, 0, -t146 * t207, t143 * t207, -t215, t59, t54, t55, 0, -t104 * t132 + t213 - t185 * t135 + (t190 * t142 - t82) * qJD(4) + t168, t103 * t132 + t226 * t135 + (t190 * qJD(4) + t86) * t145 + t163, t161 * t225 + t179 * t224 - t96 * t40 - t95 * t41 + t172, t4 * t96 - t3 * t95 - t78 * (t132 - t234) + t224 * t22 - t225 * t21, t7, t1, t6, t5, -t229, t89 * t24 + t225 * t56 + (-t3 - t244) * t144 + t156 * t141 + t191, t141 * t244 + t156 * t144 + t225 * t58 + t89 * t23 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t59, t54, t55, 0, -t169 * t135 + t154 + t213, t187 * t135 + t163 - t241, -t25 * t161 - t26 * t179 + (-t140 * t40 - t216 * t41) * pkin(4) + t172, t21 * t25 - t22 * t26 + (t103 * t78 + t140 * t4 - t216 * t3) * pkin(4), t7, t1, t6, t5, -t229, t127 * t24 - t3 * t144 - (-t141 * t26 + t144 * t33) * t238 - t25 * t56 - t19 * t219 + (-t192 - t221) * t126 + t191, t127 * t23 + (t141 * t33 + t144 * t26) * t238 - t25 * t58 - t144 * t233 + t242 * t126 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161 ^ 2 - t179 ^ 2, t161 * t21 - t179 * t22 + t150, 0, 0, 0, 0, 0, t165 - t228, -t230 + t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t56 ^ 2 + t58 ^ 2, t23 + t249, -t251 * t58 - t220, t40, -t141 * t4 - t19 * t58 - t251 * t9 + t12, -t141 * t13 - t144 * t4 + t251 * t170 + t19 * t56;];
tauc_reg  = t2;
