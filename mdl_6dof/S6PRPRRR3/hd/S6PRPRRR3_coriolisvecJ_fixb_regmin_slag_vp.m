% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:17
% EndTime: 2019-03-08 20:34:23
% DurationCPUTime: 2.78s
% Computational Cost: add. (3951->289), mult. (10466->415), div. (0->0), fcn. (8853->12), ass. (0->176)
t132 = cos(qJ(6));
t178 = qJD(6) * t132;
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t134 = cos(qJ(4));
t126 = cos(pkin(12));
t183 = qJD(2) * t126;
t168 = t134 * t183;
t124 = sin(pkin(12));
t130 = sin(qJ(4));
t192 = t130 * t124;
t169 = qJD(2) * t192;
t96 = t168 - t169;
t104 = t134 * t124 + t130 * t126;
t97 = qJD(2) * t104;
t67 = t129 * t97 - t133 * t96;
t238 = t132 * t67;
t250 = t178 + t238;
t210 = pkin(8) + qJ(3);
t108 = t210 * t126;
t190 = t134 * t126;
t103 = -t190 + t192;
t125 = sin(pkin(6));
t135 = cos(qJ(2));
t194 = t125 * t135;
t141 = t103 * t194;
t107 = t210 * t124;
t191 = t134 * t107;
t249 = -qJD(1) * t141 + (qJD(3) * t124 + qJD(4) * t108) * t130 - qJD(3) * t190 + qJD(4) * t191;
t123 = qJD(4) + qJD(5);
t197 = t67 * t123;
t180 = qJD(5) * t133;
t181 = qJD(5) * t129;
t110 = qJD(4) * t168;
t90 = -qJD(4) * t169 + t110;
t99 = t104 * qJD(4);
t91 = qJD(2) * t99;
t41 = -t129 * t91 + t133 * t90 + t96 * t180 - t97 * t181;
t248 = t41 + t197;
t142 = t104 * t194;
t149 = t130 * t107 - t134 * t108;
t247 = qJD(1) * t142 - t104 * qJD(3) + t149 * qJD(4);
t235 = -qJD(6) - t67;
t246 = qJD(6) + t235;
t128 = sin(qJ(6));
t151 = t129 * t96 + t133 * t97;
t179 = qJD(6) * t128;
t24 = t123 * t178 + t132 * t41 - t151 * t179;
t59 = t128 * t123 + t132 * t151;
t25 = qJD(6) * t59 + t128 * t41;
t57 = -t132 * t123 + t128 * t151;
t245 = -t128 * t25 + t24 * t132 - t250 * t57;
t22 = t24 * t128;
t244 = t250 * t59 + t22;
t42 = t151 * qJD(5) + t129 * t90 + t133 * t91;
t34 = t128 * t42;
t60 = t235 * t178;
t208 = t34 - t60;
t213 = t59 * t151;
t243 = -t235 * t238 + t208 - t213;
t131 = sin(qJ(2));
t185 = qJD(1) * t125;
t171 = t131 * t185;
t106 = qJD(2) * qJ(3) + t171;
t127 = cos(pkin(6));
t184 = qJD(1) * t127;
t113 = t126 * t184;
t76 = t113 + (-pkin(8) * qJD(2) - t106) * t124;
t81 = t126 * t106 + t124 * t184;
t77 = pkin(8) * t183 + t81;
t150 = -t130 * t76 - t134 * t77;
t44 = t96 * pkin(9) - t150;
t203 = t129 * t44;
t227 = -t130 * t77 + t134 * t76;
t43 = -t97 * pkin(9) + t227;
t40 = qJD(4) * pkin(4) + t43;
t16 = t133 * t40 - t203;
t14 = -t123 * pkin(5) - t16;
t242 = t14 * t67;
t118 = -t126 * pkin(3) - pkin(2);
t170 = t135 * t185;
t154 = qJD(3) - t170;
t92 = t118 * qJD(2) + t154;
t71 = -t96 * pkin(4) + t92;
t241 = t71 * t67;
t240 = -t99 * pkin(9) - t249;
t98 = t103 * qJD(4);
t239 = -t98 * pkin(9) - t247;
t237 = t151 * t67;
t236 = t128 * t235;
t198 = t151 * t123;
t234 = -t42 + t198;
t232 = t151 ^ 2 - t67 ^ 2;
t46 = pkin(5) * t151 + t67 * pkin(10);
t212 = t151 * t57;
t229 = t235 * t151;
t36 = t132 * t42;
t228 = -t179 * t235 - t36;
t102 = (qJD(3) + t170) * qJD(2);
t226 = t103 * t102;
t146 = t99 * pkin(4) - t171;
t199 = t133 * t44;
t17 = t129 * t40 + t199;
t15 = t123 * pkin(10) + t17;
t29 = t67 * pkin(5) - pkin(10) * t151 + t71;
t152 = t128 * t15 - t132 * t29;
t225 = t14 * t179 + t151 * t152;
t27 = -t91 * pkin(9) + t227 * qJD(4) - t226;
t143 = t104 * t102;
t28 = -t90 * pkin(9) + t150 * qJD(4) - t143;
t165 = t129 * t27 - t133 * t28;
t3 = t17 * qJD(5) + t165;
t5 = t128 * t29 + t132 * t15;
t224 = t3 * t128 + t14 * t178 + t5 * t151;
t164 = t129 * t28 - t44 * t181;
t2 = (qJD(5) * t40 + t27) * t133 + t164;
t61 = -t104 * pkin(9) - t130 * t108 - t191;
t62 = -t103 * pkin(9) - t149;
t31 = t129 * t62 - t133 * t61;
t217 = t31 * qJD(5) + t239 * t129 - t240 * t133;
t32 = t129 * t61 + t133 * t62;
t72 = t133 * t103 + t129 * t104;
t73 = -t129 * t103 + t133 * t104;
t85 = t103 * pkin(4) + t118;
t33 = t72 * pkin(5) - t73 * pkin(10) + t85;
t47 = -t72 * qJD(5) - t129 * t99 - t133 * t98;
t223 = -(qJD(6) * t29 + t2) * t72 + t14 * t47 + t3 * t73 - (-qJD(6) * t33 + t217) * t235 - t32 * t42;
t222 = -t71 * t151 - t165;
t220 = t97 * pkin(4);
t216 = t32 * qJD(5) + t240 * t129 + t239 * t133;
t215 = t14 * t73;
t214 = t33 * t42;
t211 = t73 * t42;
t207 = qJD(2) * pkin(2);
t205 = t128 * t59;
t195 = t125 * t131;
t136 = qJD(2) ^ 2;
t193 = t125 * t136;
t186 = t124 ^ 2 + t126 ^ 2;
t182 = qJD(2) * t131;
t175 = t73 * t179;
t174 = t131 * t193;
t173 = t135 * t193;
t172 = t125 * t182;
t109 = qJD(2) * t171;
t75 = t91 * pkin(4) + t109;
t167 = -pkin(4) * t123 - t40;
t160 = t186 * t102;
t119 = t129 * pkin(4) + pkin(10);
t158 = qJD(6) * t119 + t220 + t46;
t18 = t129 * t43 + t199;
t157 = pkin(4) * t181 - t18;
t48 = t73 * qJD(5) - t129 * t98 + t133 * t99;
t156 = t48 * pkin(5) - t47 * pkin(10) + t146;
t155 = -t235 * t47 + t211;
t153 = t124 * (-t124 * t106 + t113) - t126 * t81;
t93 = -t124 * t195 + t127 * t126;
t94 = t127 * t124 + t126 * t195;
t64 = -t130 * t94 + t134 * t93;
t65 = t130 * t93 + t134 * t94;
t38 = t129 * t65 - t133 * t64;
t39 = t129 * t64 + t133 * t65;
t148 = t236 * t67 - t228;
t145 = -t128 * t39 - t132 * t194;
t144 = t128 * t194 - t132 * t39;
t19 = t133 * t43 - t203;
t139 = -t119 * t42 + t242 - (-pkin(4) * t180 + t19) * t235;
t120 = -t133 * pkin(4) - pkin(5);
t105 = t154 - t207;
t50 = -qJD(2) * t142 - t65 * qJD(4);
t49 = -qJD(2) * t141 + t64 * qJD(4);
t11 = t42 * pkin(5) - t41 * pkin(10) + t75;
t10 = t132 * t11;
t7 = t39 * qJD(5) + t129 * t49 - t133 * t50;
t6 = -t38 * qJD(5) + t129 * t50 + t133 * t49;
t1 = [0, 0, -t174, -t173, -t126 * t174, t124 * t174, t186 * t173 (-t124 * t93 + t126 * t94) * t102 + (t105 * t131 + (-t153 - t171) * t135) * t125 * qJD(2), 0, 0, 0, 0, 0, t50 * qJD(4) + (-t135 * t91 - t182 * t96) * t125, -t49 * qJD(4) + (-t135 * t90 + t182 * t97) * t125, 0, 0, 0, 0, 0, -t7 * t123 + (-t135 * t42 + t182 * t67) * t125, -t6 * t123 + (-t135 * t41 + t151 * t182) * t125, 0, 0, 0, 0, 0 -(qJD(6) * t144 - t128 * t6 + t132 * t172) * t235 + t145 * t42 + t7 * t57 + t38 * t25 (qJD(6) * t145 + t128 * t172 + t132 * t6) * t235 + t144 * t42 + t7 * t59 + t38 * t24; 0, 0, 0, 0, 0, 0, t154 * qJD(2) * t186 + t160, -t153 * qJD(3) + qJ(3) * t160 + (t153 * t135 + (-t105 - t207) * t131) * t185, t90 * t104 - t97 * t98, -t90 * t103 - t104 * t91 - t98 * t96 - t97 * t99, -t98 * qJD(4), -t99 * qJD(4), 0, t118 * t91 + t92 * t99 + (qJD(2) * t103 + t96) * t171 + t247 * qJD(4), t249 * qJD(4) + t118 * t90 - t92 * t98, t151 * t47 + t41 * t73, -t151 * t48 - t41 * t72 - t47 * t67 - t211, t47 * t123, -t48 * t123, 0, -t216 * t123 + t146 * t67 + t85 * t42 + t71 * t48 + t75 * t72, t217 * t123 + t146 * t151 + t85 * t41 + t71 * t47 + t75 * t73, -t59 * t175 + (t24 * t73 + t47 * t59) * t132 (-t132 * t57 - t205) * t47 + (-t22 - t132 * t25 + (t128 * t57 - t132 * t59) * qJD(6)) * t73, t132 * t155 + t175 * t235 + t24 * t72 + t59 * t48, -t128 * t155 - t25 * t72 - t57 * t48 + t60 * t73, -t235 * t48 + t42 * t72, t10 * t72 + t31 * t25 - t152 * t48 + t216 * t57 + (t214 - t156 * t235 + (-t15 * t72 + t235 * t32 + t215) * qJD(6)) * t132 + t223 * t128, t31 * t24 - t5 * t48 + t216 * t59 + (-t214 - (-qJD(6) * t15 + t11) * t72 - qJD(6) * t215 - (qJD(6) * t32 - t156) * t235) * t128 + t223 * t132; 0, 0, 0, 0, 0, 0, -t186 * t136, t153 * qJD(2) + t109, 0, 0, 0, 0, 0, 0.2e1 * t97 * qJD(4), t110 + (t96 - t169) * qJD(4), 0, 0, 0, 0, 0, t42 + t198, t41 - t197, 0, 0, 0, 0, 0, t148 - t212, -t132 * t235 ^ 2 - t213 - t34; 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t96, -t96 ^ 2 + t97 ^ 2, t110 + (-t96 - t169) * qJD(4), 0, 0, -t92 * t97 - t143, -t92 * t96 + t226, t237, t232, t248, t234, 0, -t67 * t220 + t18 * t123 + (t129 * t167 - t199) * qJD(5) + t222, -t151 * t220 + t19 * t123 + t241 + (qJD(5) * t167 - t27) * t133 - t164, t244, t205 * t235 + t245, t243, t148 + t212, t229, t120 * t25 + t157 * t57 + (t158 * t235 - t3) * t132 + t139 * t128 + t225, t120 * t24 + t132 * t139 + t157 * t59 - t158 * t236 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t232, t248, t234, 0 (-qJD(5) + t123) * t17 + t222, t16 * t123 - t2 + t241, t244, t236 * t59 + t245, t243, -t235 * t236 + t212 + t36, t229, -pkin(5) * t25 - t3 * t132 + (-t128 * t16 + t132 * t46) * t235 - t17 * t57 + t128 * t242 - t208 * pkin(10) + t225, -pkin(5) * t24 - (t128 * t46 + t132 * t16) * t235 - t17 * t59 + t14 * t238 + t228 * pkin(10) + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t57, -t57 ^ 2 + t59 ^ 2, -t235 * t57 + t24, -t235 * t59 - t25, t42, -t128 * t2 - t14 * t59 - t246 * t5 + t10, -t128 * t11 - t132 * t2 + t14 * t57 + t246 * t152;];
tauc_reg  = t1;
