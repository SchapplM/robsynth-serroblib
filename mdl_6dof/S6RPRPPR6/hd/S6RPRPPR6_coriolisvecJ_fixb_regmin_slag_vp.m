% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:38
% EndTime: 2019-03-09 02:54:45
% DurationCPUTime: 2.20s
% Computational Cost: add. (3293->302), mult. (7390->438), div. (0->0), fcn. (5257->8), ass. (0->168)
t146 = sin(qJ(3));
t148 = cos(qJ(3));
t210 = sin(pkin(9));
t211 = cos(pkin(9));
t156 = t146 * t211 + t148 * t210;
t229 = t156 * qJD(1);
t99 = qJD(6) + t229;
t233 = qJD(6) - t99;
t112 = -t210 * t146 + t211 * t148;
t147 = cos(qJ(6));
t105 = t112 * qJD(1);
t143 = sin(pkin(10));
t144 = cos(pkin(10));
t81 = -t144 * qJD(3) + t143 * t105;
t232 = t147 * t81;
t145 = sin(qJ(6));
t83 = t143 * qJD(3) + t144 * t105;
t163 = t145 * t81 - t147 * t83;
t231 = t163 * t99;
t115 = t147 * t143 + t145 * t144;
t108 = t115 * qJD(6);
t204 = t147 * t144;
t205 = t145 * t143;
t114 = -t204 + t205;
t93 = t105 * qJD(3);
t230 = -t108 * t99 - t114 * t93;
t107 = t114 * qJD(6);
t213 = -t114 * t229 - t107;
t217 = t115 * t93;
t228 = -t213 * t99 - t217;
t94 = qJD(3) * t229;
t17 = -qJD(6) * t163 - t115 * t94;
t100 = t229 ^ 2;
t139 = qJD(1) * qJD(2);
t227 = 0.2e1 * t139;
t149 = -pkin(1) - pkin(7);
t226 = t144 * pkin(8);
t122 = qJD(1) * t149 + qJD(2);
t186 = t148 * qJD(4);
t191 = qJD(3) * t146;
t152 = -t122 * t191 + (qJ(4) * t191 - t186) * qJD(1);
t187 = t146 * qJD(4);
t190 = qJD(3) * t148;
t80 = t122 * t190 + (-qJ(4) * t190 - t187) * qJD(1);
t39 = -t211 * t152 + t210 * t80;
t198 = qJ(4) - t149;
t118 = t198 * t146;
t119 = t198 * t148;
t77 = -t118 * t210 + t211 * t119;
t225 = t39 * t77;
t224 = t77 * t94;
t128 = pkin(3) * t210 + qJ(5);
t223 = pkin(8) + t128;
t184 = qJD(1) * qJD(3);
t180 = t148 * t184;
t197 = pkin(3) * t180 + t139;
t32 = pkin(4) * t93 + qJ(5) * t94 - qJD(5) * t105 + t197;
t40 = t210 * t152 + t211 * t80;
t35 = qJD(3) * qJD(5) + t40;
t9 = t143 * t32 + t144 * t35;
t193 = qJD(1) * t146;
t97 = -qJ(4) * t193 + t146 * t122;
t181 = t211 * t97;
t192 = qJD(1) * t148;
t98 = -qJ(4) * t192 + t148 * t122;
t91 = qJD(3) * pkin(3) + t98;
t54 = t210 * t91 + t181;
t47 = qJD(3) * qJ(5) + t54;
t117 = pkin(3) * t193 + qJD(1) * qJ(2) + qJD(4);
t55 = pkin(4) * t229 - qJ(5) * t105 + t117;
t20 = t143 * t55 + t144 * t47;
t176 = qJD(3) * t210;
t177 = qJD(3) * t211;
t103 = t146 * t176 - t148 * t177;
t104 = -t146 * t177 - t148 * t176;
t185 = pkin(3) * t190 + qJD(2);
t41 = -pkin(4) * t103 - qJ(5) * t104 - qJD(5) * t112 + t185;
t95 = t191 * t198 - t186;
t96 = -qJD(3) * t119 - t187;
t61 = t210 * t95 + t211 * t96;
t15 = t143 * t41 + t144 * t61;
t88 = t210 * t97;
t63 = t211 * t98 - t88;
t64 = pkin(3) * t192 + t105 * pkin(4) + qJ(5) * t229;
t22 = t143 * t64 + t144 * t63;
t199 = t146 * pkin(3) + qJ(2);
t69 = pkin(4) * t156 - t112 * qJ(5) + t199;
t78 = -t118 * t211 - t119 * t210;
t29 = t143 * t69 + t144 * t78;
t188 = qJD(6) * t147;
t222 = -t81 * t188 - t94 * t204;
t221 = t103 * t99;
t42 = t145 * t83 + t232;
t220 = t105 * t42;
t219 = t105 * t163;
t218 = t156 * t93;
t216 = t143 * t94;
t215 = t144 * t94;
t214 = t39 * t112;
t58 = t115 * t229;
t212 = t108 + t58;
t209 = t229 * t143;
t208 = t104 * t143;
t207 = t112 * t143;
t206 = t112 * t144;
t150 = qJD(3) ^ 2;
t203 = t150 * t146;
t202 = t150 * t148;
t151 = qJD(1) ^ 2;
t201 = t151 * qJ(2);
t200 = t151 * t148;
t196 = t146 ^ 2 - t148 ^ 2;
t195 = -t150 - t151;
t194 = qJD(1) * t229;
t189 = qJD(6) * t112;
t183 = 0.2e1 * qJD(1);
t8 = -t143 * t35 + t144 * t32;
t4 = t93 * pkin(5) + pkin(8) * t215 + t8;
t5 = pkin(8) * t216 + t9;
t182 = -t145 * t5 + t147 * t4;
t14 = -t143 * t61 + t144 * t41;
t19 = -t143 * t47 + t144 * t55;
t21 = -t143 * t63 + t144 * t64;
t28 = -t143 * t78 + t144 * t69;
t175 = -t58 * t99 + t230;
t133 = -pkin(3) * t211 - pkin(4);
t172 = t145 * t4 + t147 * t5;
t53 = t211 * t91 - t88;
t60 = t210 * t96 - t211 * t95;
t62 = t210 * t98 + t181;
t171 = -t19 * t103 + t156 * t8;
t170 = t20 * t103 - t156 * t9;
t12 = -pkin(8) * t81 + t20;
t7 = pkin(5) * t229 - t83 * pkin(8) + t19;
t2 = t147 * t12 + t145 * t7;
t169 = t145 * t12 - t147 * t7;
t46 = -qJD(3) * pkin(4) + qJD(5) - t53;
t168 = -t46 * t104 - t214;
t167 = -t143 * t19 + t144 * t20;
t166 = t143 * t83 - t144 * t81;
t18 = pkin(5) * t156 - pkin(8) * t206 + t28;
t23 = -pkin(8) * t207 + t29;
t165 = t145 * t23 - t147 * t18;
t164 = t145 * t18 + t147 * t23;
t162 = -qJD(6) * t83 + t216;
t109 = t223 * t143;
t161 = pkin(8) * t209 - qJD(5) * t144 + qJD(6) * t109 + t22;
t110 = t223 * t144;
t160 = t105 * pkin(5) + qJD(5) * t143 + qJD(6) * t110 + t226 * t229 + t21;
t159 = t115 * t99;
t158 = -t168 - t224;
t157 = t103 * t229 + t112 * t94 - t218;
t16 = t162 * t145 + t222;
t155 = -t128 * t93 - t133 * t94 + (-qJD(5) + t46) * t229;
t153 = -t103 * t54 + t104 * t53 + t156 * t40 - t214;
t120 = -t144 * pkin(5) + t133;
t66 = t114 * t112;
t65 = t115 * t112;
t51 = pkin(5) * t207 + t77;
t36 = -pkin(5) * t209 + t62;
t34 = pkin(5) * t208 + t60;
t27 = t81 * pkin(5) + t46;
t26 = -pkin(5) * t216 + t39;
t25 = t104 * t115 + t188 * t206 - t189 * t205;
t24 = -t104 * t114 - t115 * t189;
t10 = -pkin(8) * t208 + t15;
t6 = -t103 * pkin(5) - t104 * t226 + t14;
t1 = [0, 0, 0, 0, t227, qJ(2) * t227, -0.2e1 * t146 * t180, 0.2e1 * t196 * t184, -t203, -t202, 0, -t149 * t203 + (qJ(2) * t190 + qJD(2) * t146) * t183, -t149 * t202 + (-qJ(2) * t191 + qJD(2) * t148) * t183, t105 * t60 - t229 * t61 - t78 * t93 - t153 - t224, t117 * t185 + t197 * t199 + t40 * t78 - t53 * t60 + t54 * t61 + t225, t14 * t229 + t143 * t158 + t28 * t93 + t60 * t81 + t171, t144 * t158 - t15 * t229 - t29 * t93 + t60 * t83 + t170, -t14 * t83 - t15 * t81 + (-t104 * t19 - t112 * t8 + t28 * t94) * t144 + (-t104 * t20 - t112 * t9 + t29 * t94) * t143, t14 * t19 + t15 * t20 + t28 * t8 + t29 * t9 + t46 * t60 + t225, -t16 * t66 - t163 * t24, -t16 * t65 + t163 * t25 + t17 * t66 - t24 * t42, t103 * t163 + t156 * t16 + t24 * t99 - t66 * t93, t103 * t42 - t156 * t17 - t25 * t99 - t65 * t93, t218 - t221 (-t145 * t10 + t147 * t6) * t99 - t165 * t93 + t182 * t156 + t169 * t103 + t34 * t42 + t51 * t17 + t26 * t65 + t27 * t25 + (-t156 * t2 - t164 * t99) * qJD(6) -(t147 * t10 + t145 * t6) * t99 - t164 * t93 - t172 * t156 + t2 * t103 - t34 * t163 + t51 * t16 - t26 * t66 + t27 * t24 + (t156 * t169 + t165 * t99) * qJD(6); 0, 0, 0, 0, -t151, -t201, 0, 0, 0, 0, 0, t195 * t146, t195 * t148, -t104 * t105 + t157, -qJD(1) * t117 + t153, -t104 * t81 + t143 * t157 - t144 * t194, -t104 * t83 + t143 * t194 + t144 * t157, -t166 * t103 + (t143 * t81 + t144 * t83) * qJD(1) (-qJD(1) * t19 - t170) * t144 + (-qJD(1) * t20 - t171) * t143 + t168, 0, 0, 0, 0, 0, -t104 * t42 - t112 * t17 + t103 * t159 + t114 * t99 * qJD(1) - (-t107 * t99 + t217) * t156, qJD(1) * t159 + t104 * t163 - t112 * t16 - t114 * t221 - t156 * t230; 0, 0, 0, 0, 0, 0, t146 * t200, -t196 * t151, 0, 0, 0, -qJ(2) * t200, t146 * t201 (t54 - t62) * t105 + (-t53 + t63) * t229 + (-t210 * t93 + t211 * t94) * pkin(3), t53 * t62 - t54 * t63 + (-t117 * t192 + t210 * t40 - t211 * t39) * pkin(3), -t19 * t105 + t143 * t155 - t39 * t144 - t21 * t229 - t62 * t81, t20 * t105 + t39 * t143 + t144 * t155 + t22 * t229 - t62 * t83, t21 * t83 + t22 * t81 + (-qJD(5) * t81 - t19 * t229 + t9) * t144 + (qJD(5) * t83 - t20 * t229 - t8) * t143, t39 * t133 - t19 * t21 - t20 * t22 - t46 * t62 + (-t8 * t143 + t9 * t144) * t128 + t167 * qJD(5), t115 * t16 - t163 * t213, -t114 * t16 - t115 * t17 + t163 * t212 - t213 * t42, t219 - t228, t175 + t220, -t99 * t105 (-t147 * t109 - t145 * t110) * t93 + t120 * t17 + t26 * t114 + t169 * t105 - t36 * t42 + (t145 * t161 - t147 * t160) * t99 + t212 * t27 -(-t145 * t109 + t147 * t110) * t93 + t120 * t16 + t26 * t115 + t2 * t105 + t36 * t163 + (t145 * t160 + t147 * t161) * t99 + t213 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 ^ 2 - t100, t105 * t53 + t229 * t54 + t197, -t100 * t143 - t105 * t81 + t144 * t93, -t144 * t100 - t105 * t83 - t143 * t93 (t143 ^ 2 + t144 ^ 2) * t94 + t166 * t229, -t46 * t105 + t9 * t143 + t8 * t144 + t167 * t229, 0, 0, 0, 0, 0, t175 - t220, t219 + t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229 * t83 - t216, -t229 * t81 - t215, -t81 ^ 2 - t83 ^ 2, t19 * t83 + t20 * t81 + t39, 0, 0, 0, 0, 0, t17 - t231, -t99 * t232 + (-t83 * t99 + t162) * t145 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163 * t42, t163 ^ 2 - t42 ^ 2, t42 * t99 + t16, -t17 - t231, t93, t27 * t163 - t233 * t2 + t182, t233 * t169 + t27 * t42 - t172;];
tauc_reg  = t1;
