% Calculate minimal parameter regressor of coriolis matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:45
% EndTime: 2019-03-09 01:37:49
% DurationCPUTime: 1.84s
% Computational Cost: add. (776->218), mult. (1715->312), div. (0->0), fcn. (1535->6), ass. (0->180)
t116 = sin(qJ(6));
t225 = 0.2e1 * t116;
t118 = cos(qJ(6));
t117 = sin(qJ(5));
t177 = t117 * qJD(1);
t169 = t118 * t177;
t108 = t116 ^ 2;
t110 = t118 ^ 2;
t94 = t110 - t108;
t123 = -t94 * qJD(5) + t169 * t225;
t112 = sin(pkin(9));
t119 = cos(qJ(5));
t201 = t118 * t119;
t113 = cos(pkin(9));
t206 = t116 * t113;
t54 = t112 * t201 - t206;
t224 = -t54 / 0.2e1;
t207 = t116 * t112;
t55 = t113 * t201 + t207;
t223 = -t55 / 0.2e1;
t220 = t119 * pkin(8);
t221 = t117 * pkin(5);
t82 = -t220 + t221;
t222 = -t82 / 0.2e1;
t217 = t118 * t82;
t114 = pkin(3) + qJ(2);
t115 = pkin(1) + qJ(3);
t61 = t112 * t114 - t113 * t115;
t57 = pkin(7) + t61;
t218 = t117 * t57;
t45 = t116 * t218;
t205 = t116 * t119;
t171 = t57 * t205;
t139 = -t119 * pkin(5) - t117 * pkin(8);
t146 = t112 * t115 + t113 * t114;
t56 = -pkin(4) - t146;
t35 = t139 + t56;
t8 = -t118 * t35 + t171;
t1 = t8 * t117 + (-t45 + t217) * t119;
t219 = t1 * qJD(1);
t170 = t57 * t201;
t9 = t116 * t35 + t170;
t2 = t82 * t205 + (-t9 + t170) * t117;
t216 = t2 * qJD(1);
t109 = t117 ^ 2;
t3 = -t109 * t57 * t116 - t8 * t119;
t215 = t3 * qJD(1);
t204 = t118 * t109;
t4 = -t9 * t119 - t57 * t204;
t214 = t4 * qJD(1);
t202 = t118 * t113;
t52 = t112 * t205 + t202;
t213 = t52 * t119;
t203 = t118 * t112;
t53 = t113 * t205 - t203;
t212 = t53 * t119;
t211 = t54 * t119;
t210 = t55 * t119;
t172 = 0.1e1 / 0.2e1 + t109 / 0.2e1;
t140 = t172 * t118;
t156 = -t207 / 0.2e1;
t10 = (t156 + t223) * t119 - t113 * t140;
t209 = t10 * qJD(1);
t153 = t206 / 0.2e1;
t11 = (t153 + t224) * t119 - t112 * t140;
t208 = t11 * qJD(1);
t151 = -t203 / 0.2e1;
t12 = (t151 + t53 / 0.2e1) * t119 + t172 * t206;
t200 = t12 * qJD(1);
t148 = t202 / 0.2e1;
t155 = t207 / 0.2e1;
t83 = t109 * t207;
t13 = t155 + t83 / 0.2e1 + (t148 + t52 / 0.2e1) * t119;
t199 = t13 * qJD(1);
t152 = t205 / 0.2e1;
t135 = t112 * t152 - t52 / 0.2e1;
t149 = -t202 / 0.2e1;
t19 = (t149 + t135) * t117;
t198 = t19 * qJD(1);
t147 = t201 / 0.2e1;
t133 = t112 * t147 + t224;
t20 = (t153 + t133) * t117;
t197 = t20 * qJD(1);
t142 = t113 * t152;
t134 = t142 - t53 / 0.2e1;
t150 = t203 / 0.2e1;
t22 = (t150 + t134) * t117;
t196 = t22 * qJD(1);
t141 = t113 * t147;
t132 = t141 + t223;
t25 = (t156 + t132) * t117;
t195 = t25 * qJD(1);
t26 = t61 * t112 + t146 * t113;
t194 = t26 * qJD(1);
t27 = -t146 * t112 + t61 * t113;
t193 = t27 * qJD(1);
t30 = -t83 - t213;
t192 = t30 * qJD(1);
t31 = t109 * t206 + t212;
t191 = t31 * qJD(1);
t32 = t109 * t203 + t211;
t190 = t32 * qJD(1);
t33 = -t109 * t202 - t210;
t189 = t33 * qJD(1);
t111 = t119 ^ 2;
t95 = t111 - t109;
t68 = t95 * t116;
t188 = t68 * qJD(1);
t69 = t118 * t111 - t204;
t187 = t69 * qJD(1);
t88 = t112 ^ 2 + t113 ^ 2;
t186 = t88 * qJD(1);
t184 = t95 * qJD(1);
t183 = qJD(5) * t116;
t182 = qJD(5) * t118;
t181 = qJD(6) * t116;
t180 = qJD(6) * t118;
t179 = qJD(6) * t119;
t178 = t115 * qJD(1);
t176 = t117 * qJD(5);
t175 = t117 * qJD(6);
t174 = t119 * qJD(1);
t173 = t119 * qJD(5);
t168 = t116 * t179;
t167 = t118 * t179;
t166 = t116 * t180;
t165 = t116 * t182;
t164 = t117 * t173;
t163 = t112 * t177;
t162 = t112 * t176;
t161 = t113 * t177;
t160 = t113 * t176;
t159 = t117 * t174;
t158 = t118 * t176;
t157 = t112 * t174;
t154 = -t206 / 0.2e1;
t145 = -qJD(6) + t174;
t143 = t116 * t158;
t138 = t145 * t117;
t137 = t113 * qJD(2) + t112 * qJD(3);
t136 = -t112 * qJD(2) + t113 * qJD(3);
t131 = t220 / 0.2e1 - t221 / 0.2e1;
t126 = t222 + t131;
t29 = t126 * t118;
t130 = pkin(5) * t183 - t29 * qJD(1);
t28 = t126 * t116;
t129 = pkin(5) * t182 + t28 * qJD(1);
t58 = (t108 / 0.2e1 - t110 / 0.2e1) * t117;
t128 = -t58 * qJD(1) + t165;
t127 = t118 * t138;
t125 = t116 * qJD(1) * t204 + t58 * qJD(5);
t67 = t94 * t109;
t124 = t67 * qJD(1) + 0.2e1 * t143;
t122 = t116 * t175 - t118 * t173;
t121 = t116 * t173 + t118 * t175;
t120 = -qJD(1) * t56 + t136;
t107 = qJ(2) * qJD(1);
t106 = qJ(2) * qJD(2);
t105 = t113 * qJD(1);
t104 = t112 * qJD(1);
t102 = t176 / 0.2e1;
t98 = t116 * t176;
t93 = t113 * t174;
t62 = (t174 - qJD(6) / 0.2e1) * t117;
t51 = t58 * qJD(6);
t50 = -t93 + t162;
t49 = t157 + t160;
t48 = t112 * t173 + t161;
t47 = -t113 * t173 + t163;
t24 = (t132 + t155) * t117;
t23 = (t134 + t151) * t117;
t21 = (t133 + t154) * t117;
t18 = (t135 + t148) * t117;
t17 = t210 / 0.2e1 + t109 * t148 + t119 * t156 + t149;
t16 = t211 / 0.2e1 + t109 * t150 + t142 + t151;
t15 = -t212 / 0.2e1 + t109 * t154 + t119 * t151 + t153;
t14 = -t213 / 0.2e1 - t83 / 0.2e1 + t141 + t155;
t7 = t45 + t217 / 0.2e1 + t131 * t118;
t6 = t118 * t218 + (-t131 + t222) * t116;
t5 = [0, 0, 0, 0, qJD(2), t106, qJD(2), qJD(3), t115 * qJD(3) + t106, t137, t136, t26 * qJD(2) - t27 * qJD(3), t164, t95 * qJD(5), 0, 0, 0, t137 * t119 + t56 * t176, -t137 * t117 + t56 * t173, -t109 * t166 + t110 * t164, -t67 * qJD(6) - 0.2e1 * t119 * t143, -t69 * qJD(5) + t117 * t168, t68 * qJD(5) + t117 * t167, -t164, -t30 * qJD(2) - t31 * qJD(3) - t1 * qJD(5) - t4 * qJD(6), t32 * qJD(2) + t33 * qJD(3) + t2 * qJD(5) + t3 * qJD(6); 0, 0, 0, 0, qJD(1), t107, qJD(1), 0, t107, t105, -t104, t194, 0, 0, 0, 0, 0, t93, -t161, 0, 0, 0, 0, 0, t23 * qJD(5) + t17 * qJD(6) - t192, t24 * qJD(5) + t15 * qJD(6) + t190; 0, 0, 0, 0, 0, 0, 0, qJD(1), t178, t104, t105, -t193, 0, 0, 0, 0, 0, t157, -t163, 0, 0, 0, 0, 0, t18 * qJD(5) + t16 * qJD(6) - t191, t21 * qJD(5) + t14 * qJD(6) + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t184, t173, -t176, 0, -t57 * t173 + t56 * t177, t56 * t174 + t176 * t57, -t51 + (t110 * t177 + t165) * t119, -0.2e1 * t117 * t166 - t119 * t123, t98 - t187, t158 + t188, -t62, -t219 + t23 * qJD(2) + t18 * qJD(3) + (t116 * t139 - t170) * qJD(5) + t7 * qJD(6), t216 + t24 * qJD(2) + t21 * qJD(3) + (t118 * t139 + t171) * qJD(5) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124, t116 * t138, t127, t102, t17 * qJD(2) + t16 * qJD(3) + t7 * qJD(5) - t9 * qJD(6) - t214, t15 * qJD(2) + t14 * qJD(3) + t6 * qJD(5) + t8 * qJD(6) + t215; 0, 0, 0, 0, -qJD(1), -t107, -qJD(1), 0, -qJD(3) - t107, -t105, t104, -t88 * qJD(3) - t194, 0, 0, 0, 0, 0, t50, t48, 0, 0, 0, 0, 0, t22 * qJD(5) - t10 * qJD(6) + t192, t25 * qJD(5) - t12 * qJD(6) - t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, -t186, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t49, 0, 0, 0, 0, 0, t113 * t122 + t196, t113 * t121 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * qJD(6) + t116 * t160 - t209, t53 * qJD(6) + t113 * t158 - t200; 0, 0, 0, 0, 0, 0, 0, -qJD(1), qJD(2) - t178, -t104, -t105, t88 * qJD(2) + t193, 0, 0, 0, 0, 0, -t49, t47, 0, 0, 0, 0, 0, t19 * qJD(5) - t11 * qJD(6) + t191, t20 * qJD(5) - t13 * qJD(6) - t189; 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, t186, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t50, 0, 0, 0, 0, 0, t112 * t122 + t198, t112 * t121 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 * qJD(6) + t116 * t162 - t208, t52 * qJD(6) + t112 * t158 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t173, 0, 0, 0, 0, 0, -t158 - t168, t98 - t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, -t184, 0, 0, 0, t120 * t117, t120 * t119, -t110 * t159 - t51, t127 * t225, -t167 + t187, t168 - t188, t62, -t22 * qJD(2) - t19 * qJD(3) + t29 * qJD(6) + t219, -t25 * qJD(2) - t20 * qJD(3) - t28 * qJD(6) - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, -t157, 0, 0, 0, 0, 0, -t196, -t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t93, 0, 0, 0, 0, 0, -t198, -t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t94 * qJD(6), 0, 0, 0, -pkin(5) * t181, -pkin(5) * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t123, -t145 * t118, t145 * t116, -t177 / 0.2e1, -pkin(8) * t180 - t130, pkin(8) * t181 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t124 (-t116 * t177 + t182) * t119 (-t169 - t183) * t119, t102, t10 * qJD(2) + t11 * qJD(3) - t29 * qJD(5) + t214, t12 * qJD(2) + t13 * qJD(3) + t28 * qJD(5) - t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t123, t118 * t174, -t116 * t174, t177 / 0.2e1, t130, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
