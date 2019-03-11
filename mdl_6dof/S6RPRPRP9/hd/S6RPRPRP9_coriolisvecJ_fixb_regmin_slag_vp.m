% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:37
% EndTime: 2019-03-09 03:29:42
% DurationCPUTime: 2.28s
% Computational Cost: add. (3816->339), mult. (8270->466), div. (0->0), fcn. (5427->6), ass. (0->170)
t140 = sin(pkin(9));
t145 = cos(qJ(3));
t202 = qJD(1) * t145;
t182 = t140 * t202;
t141 = cos(pkin(9));
t195 = t141 * qJD(3);
t109 = -t182 + t195;
t181 = t141 * t202;
t196 = t140 * qJD(3);
t110 = t181 + t196;
t142 = sin(qJ(5));
t144 = cos(qJ(5));
t58 = -t144 * t109 + t142 * t110;
t240 = t58 ^ 2;
t143 = sin(qJ(3));
t194 = t143 * qJD(1);
t189 = t140 * t194;
t211 = t144 * t141;
t197 = qJD(5) * t144;
t198 = qJD(5) * t142;
t235 = -t140 * t198 + t141 * t197;
t220 = -t142 * t189 + t194 * t211 + t235;
t135 = qJD(5) + t194;
t239 = t58 * t135;
t113 = t144 * t140 + t142 * t141;
t100 = t113 * qJD(5);
t98 = t113 * qJD(1);
t216 = t143 * t98 + t100;
t161 = t142 * t109 + t144 * t110;
t231 = t161 ^ 2;
t192 = 0.2e1 * qJD(1);
t229 = t141 * pkin(8);
t191 = t143 * t229;
t153 = (pkin(4) * t145 + t191) * qJD(1);
t169 = pkin(3) * t145 + qJ(4) * t143;
t115 = t169 * qJD(1);
t146 = -pkin(1) - pkin(7);
t234 = qJD(1) * t146;
t130 = qJD(2) + t234;
t210 = t145 * t130;
t68 = t141 * t115 - t140 * t210;
t41 = t153 + t68;
t69 = t140 * t115 + t141 * t210;
t50 = pkin(8) * t189 + t69;
t227 = pkin(8) + qJ(4);
t126 = t227 * t140;
t127 = t227 * t141;
t73 = -t142 * t126 + t144 * t127;
t238 = -qJD(4) * t113 - qJD(5) * t73 + t142 * t50 - t144 * t41;
t112 = t142 * t140 - t211;
t160 = -t144 * t126 - t142 * t127;
t237 = qJD(4) * t112 - qJD(5) * t160 + t142 * t41 + t144 * t50;
t236 = t135 * t161;
t193 = qJD(1) * qJD(3);
t180 = t143 * t193;
t171 = t144 * t180;
t172 = t142 * t180;
t27 = qJD(5) * t161 - t140 * t171 - t141 * t172;
t120 = t143 * pkin(3) - t145 * qJ(4) + qJ(2);
t107 = t141 * t120;
t178 = -t140 * t146 + pkin(4);
t56 = t143 * t178 - t145 * t229 + t107;
t213 = t140 * t145;
t212 = t143 * t146;
t76 = t140 * t120 + t141 * t212;
t67 = -pkin(8) * t213 + t76;
t162 = t142 * t56 + t144 * t67;
t91 = qJD(3) * t169 - t145 * qJD(4) + qJD(2);
t79 = t141 * t91;
t38 = t79 + (t145 * t178 + t191) * qJD(3);
t188 = t143 * t196;
t199 = qJD(3) * t146;
t185 = t145 * t199;
t65 = t140 * t91 + t141 * t185;
t48 = pkin(8) * t188 + t65;
t233 = -qJD(5) * t162 - t142 * t48 + t144 * t38;
t90 = t112 * t145;
t222 = -qJD(3) * t90 - t143 * t100 - t98;
t26 = -t109 * t197 + t110 * t198 - t140 * t172 + t141 * t171;
t89 = t112 * t143;
t232 = qJD(3) * (-t143 * t161 - t202 * t89) + t135 * t222 - t145 * t26;
t219 = qJD(3) * pkin(3);
t176 = -qJD(4) + t219;
t95 = -t176 - t210;
t66 = -t109 * pkin(4) + t95;
t11 = t58 * pkin(5) - qJ(6) * t161 + t66;
t230 = t11 * t161;
t228 = t161 * t58;
t226 = qJ(6) * t202 + t237;
t225 = -pkin(5) * t202 + t238;
t119 = t143 * t130;
t82 = -pkin(4) * t189 + t119;
t224 = -t216 * pkin(5) + t220 * qJ(6) + t113 * qJD(6) + t82;
t88 = t113 * t145;
t221 = -t112 * qJD(1) + qJD(3) * t88 - qJD(5) * t89;
t74 = t91 * qJD(1);
t93 = (qJD(4) + t210) * qJD(3);
t36 = t140 * t74 + t141 * t93;
t103 = t120 * qJD(1);
t104 = qJD(3) * qJ(4) + t119;
t52 = t140 * t103 + t141 * t104;
t51 = t141 * t103 - t140 * t104;
t28 = pkin(4) * t194 - t110 * pkin(8) + t51;
t32 = t109 * pkin(8) + t52;
t9 = -t142 * t32 + t144 * t28;
t217 = qJD(6) - t9;
t215 = qJD(3) * t160;
t214 = qJD(3) * t73;
t147 = qJD(3) ^ 2;
t209 = t147 * t143;
t208 = t147 * t145;
t148 = qJD(1) ^ 2;
t207 = t148 * qJ(2);
t206 = t148 * t145;
t204 = t143 ^ 2 - t145 ^ 2;
t203 = -t147 - t148;
t201 = qJD(3) * t143;
t200 = qJD(3) * t145;
t190 = qJD(2) * t192;
t137 = -t141 * pkin(4) - pkin(3);
t187 = t142 * t201;
t186 = t144 * t201;
t134 = t143 * t199;
t179 = t145 * t193;
t35 = -t140 * t93 + t141 * t74;
t21 = qJD(3) * t153 + t35;
t173 = t140 * t180;
t29 = pkin(8) * t173 + t36;
t177 = t142 * t29 - t144 * t21 + t32 * t197 + t28 * t198;
t108 = pkin(4) * t213 - t145 * t146;
t175 = -t110 + t196;
t174 = pkin(5) * t179;
t170 = qJ(6) * t179;
t168 = -t35 * t140 + t36 * t141;
t167 = -t140 * t52 - t141 * t51;
t166 = -t140 * t51 + t141 * t52;
t10 = t142 * t28 + t144 * t32;
t163 = -t142 * t67 + t144 * t56;
t159 = (t109 - t195) * t143;
t96 = -pkin(4) * t188 + t134;
t158 = t10 * t135 - t177;
t156 = -t142 * t21 - t144 * t29 - t28 * t197 + t198 * t32;
t155 = t142 * t38 + t144 * t48 + t56 * t197 - t198 * t67;
t114 = t130 * t201;
t77 = -pkin(4) * t173 + t114;
t154 = -t95 + (t130 + t234) * t145;
t2 = -t174 + t177;
t152 = -qJ(4) * t200 + (t176 + t95) * t143;
t87 = t113 * t143;
t151 = t58 * t201 + (-t193 * t87 - t27) * t145 - t221 * t135;
t150 = t27 + t236;
t3 = t27 * pkin(5) + t26 * qJ(6) - qJD(6) * t161 + t77;
t149 = t26 + t239;
t75 = -t140 * t212 + t107;
t64 = -t140 * t185 + t79;
t55 = t112 * pkin(5) - t113 * qJ(6) + t137;
t47 = -t140 * t186 - t141 * t187 + t235 * t145;
t45 = t100 * t145 - t140 * t187 + t141 * t186;
t31 = t88 * pkin(5) + t90 * qJ(6) + t108;
t17 = pkin(5) * t161 + t58 * qJ(6);
t16 = -t143 * pkin(5) - t163;
t15 = t143 * qJ(6) + t162;
t12 = -t26 + t239;
t8 = t47 * pkin(5) + t45 * qJ(6) + t90 * qJD(6) + t96;
t7 = t135 * qJ(6) + t10;
t6 = -t135 * pkin(5) + t217;
t5 = -pkin(5) * t200 - t233;
t4 = qJ(6) * t200 + t143 * qJD(6) + t155;
t1 = t135 * qJD(6) - t156 + t170;
t13 = [0, 0, 0, 0, t190, qJ(2) * t190, -0.2e1 * t143 * t179, 0.2e1 * t204 * t193, -t209, -t208, 0, -t146 * t209 + (qJ(2) * t200 + qJD(2) * t143) * t192, -t146 * t208 + (-qJ(2) * t201 + qJD(2) * t145) * t192 (qJD(1) * t64 + t35) * t143 + ((qJD(1) * t75 + t51) * t145 + (-t109 * t146 + t140 * t154) * t143) * qJD(3) (-qJD(1) * t65 - t36) * t143 + ((-qJD(1) * t76 - t52) * t145 + (t110 * t146 + t141 * t154) * t143) * qJD(3), t65 * t109 - t64 * t110 + (-t140 * t36 - t141 * t35) * t145 + ((t140 * t76 + t141 * t75) * qJD(1) - t167) * t201, t35 * t75 + t36 * t76 + t51 * t64 + t52 * t65 + (t95 - t210) * t134, -t161 * t45 + t26 * t90, -t161 * t47 + t26 * t88 + t90 * t27 + t45 * t58, -t45 * t135 - t26 * t143 + (-qJD(1) * t90 + t161) * t200, -t47 * t135 - t27 * t143 + (-qJD(1) * t88 - t58) * t200 (t135 + t194) * t200, t233 * t135 - t177 * t143 + t96 * t58 + t108 * t27 + t77 * t88 + t66 * t47 + (qJD(1) * t163 + t9) * t200, -t155 * t135 + t156 * t143 + t96 * t161 - t108 * t26 - t77 * t90 - t66 * t45 + (-qJD(1) * t162 - t10) * t200, t11 * t47 - t5 * t135 - t2 * t143 + t31 * t27 + t3 * t88 + t8 * t58 + (-qJD(1) * t16 - t6) * t200, -t1 * t88 - t15 * t27 - t16 * t26 + t161 * t5 - t2 * t90 - t4 * t58 - t6 * t45 - t7 * t47, t1 * t143 + t11 * t45 + t4 * t135 + t31 * t26 + t3 * t90 - t8 * t161 + (qJD(1) * t15 + t7) * t200, t1 * t15 + t11 * t8 + t2 * t16 + t3 * t31 + t7 * t4 + t6 * t5; 0, 0, 0, 0, -t148, -t207, 0, 0, 0, 0, 0, t203 * t143, t203 * t145 (-t141 * t148 + (-t109 - t182) * qJD(3)) * t143 (t140 * t148 + (t110 - t181) * qJD(3)) * t143 (t109 * t141 + t110 * t140) * t200 + (-t109 * t140 + t110 * t141) * qJD(1), t168 * t143 + t167 * qJD(1) + (t143 * t95 + (t166 - t119) * t145) * qJD(3), 0, 0, 0, 0, 0, t151, -t232, t151, t161 * t221 - t222 * t58 - t87 * t26 + t89 * t27, t232, -t1 * t89 + t11 * t201 - t3 * t145 + t2 * t87 + t221 * t6 + t222 * t7; 0, 0, 0, 0, 0, 0, t143 * t206, -t204 * t148, 0, 0, 0, -qJ(2) * t206, t143 * t207, t130 * t159 + (t140 * t152 - t143 * t68 - t145 * t51) * qJD(1), t175 * t119 + (t141 * t152 + t143 * t69 + t145 * t52) * qJD(1), -t69 * t109 + t68 * t110 + (qJD(4) * t109 - t194 * t51 + t36) * t141 + (qJD(4) * t110 - t194 * t52 - t35) * t140, -t51 * t68 - t52 * t69 + (-t95 - t219) * t119 + t166 * qJD(4) + t168 * qJ(4), -t26 * t113 + t161 * t220, t26 * t112 - t113 * t27 - t161 * t216 - t220 * t58, t220 * t135 + (qJD(3) * t113 - t161) * t202, -t216 * t135 + (-qJD(3) * t112 + t58) * t202, -t135 * t202, t77 * t112 + t137 * t27 - t82 * t58 + t216 * t66 + t238 * t135 + (-t9 + t215) * t202, t77 * t113 - t137 * t26 - t82 * t161 + t220 * t66 + t237 * t135 + (t10 - t214) * t202, t3 * t112 + t55 * t27 - t224 * t58 + t225 * t135 + t216 * t11 + (t6 + t215) * t202, -t1 * t112 + t2 * t113 + t160 * t26 - t161 * t225 - t216 * t7 + t220 * t6 + t226 * t58 - t73 * t27, -t3 * t113 + t55 * t26 + t224 * t161 - t226 * t135 - t220 * t11 + (-t7 + t214) * t202, t1 * t73 - t224 * t11 - t160 * t2 - t225 * t6 - t226 * t7 + t3 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175 * t194, qJD(1) * t159, -t109 ^ 2 - t110 ^ 2, -t52 * t109 + t51 * t110 + t114, 0, 0, 0, 0, 0, t150, -t149, t150, -t231 - t240, t149, -t161 * t6 + t58 * t7 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t231 - t240, t12, -t27 + t236, t179, -t161 * t66 + t158, t9 * t135 + t66 * t58 + t156, -t17 * t58 + t158 + 0.2e1 * t174 - t230, pkin(5) * t26 - t27 * qJ(6) + (-t10 + t7) * t161 + (t6 - t217) * t58, 0.2e1 * t170 - t11 * t58 + t17 * t161 + (0.2e1 * qJD(6) - t9) * t135 - t156, -t2 * pkin(5) + t1 * qJ(6) - t6 * t10 - t11 * t17 + t217 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179 + t228, t12, -t135 ^ 2 - t231, -t7 * t135 + t2 + t230;];
tauc_reg  = t13;
