% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP12_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:22
% EndTime: 2019-12-31 18:57:28
% DurationCPUTime: 1.70s
% Computational Cost: add. (1716->237), mult. (3285->353), div. (0->0), fcn. (2684->4), ass. (0->188)
t135 = sin(qJ(4));
t249 = 0.2e1 * t135;
t136 = sin(qJ(3));
t139 = -pkin(1) - pkin(6);
t211 = t135 * t139;
t162 = pkin(4) - t211;
t137 = cos(qJ(4));
t138 = cos(qJ(3));
t207 = t137 * t138;
t178 = qJ(5) * t207;
t158 = pkin(3) * t136 - pkin(7) * t138;
t95 = qJ(2) + t158;
t83 = t137 * t95;
t47 = t136 * t162 - t178 + t83;
t246 = -t47 / 0.2e1;
t209 = t136 * t139;
t58 = t135 * t209 - t83;
t50 = -t58 - t178;
t180 = t50 / 0.2e1 + t246;
t230 = -qJ(5) - pkin(7);
t97 = t230 * t137;
t248 = t180 * t97;
t131 = t135 ^ 2;
t133 = t137 ^ 2;
t109 = t133 + t131;
t110 = t133 - t131;
t160 = t207 * t249;
t142 = qJD(1) * t160 - qJD(3) * t110;
t210 = t136 * t137;
t231 = t138 * pkin(3);
t233 = t136 * pkin(7);
t98 = t231 + t233;
t93 = t137 * t98;
t48 = qJ(5) * t210 + t138 * t162 + t93;
t245 = t48 / 0.2e1;
t212 = t135 * t138;
t206 = t137 * t139;
t179 = t136 * t206;
t59 = t135 * t95 + t179;
t51 = -qJ(5) * t212 + t59;
t243 = t51 / 0.2e1;
t134 = t138 ^ 2;
t125 = t134 * t137;
t242 = -t125 / 0.2e1;
t241 = -t133 / 0.2e1;
t240 = -t135 / 0.2e1;
t238 = t136 / 0.2e1;
t237 = t137 / 0.2e1;
t236 = -t138 / 0.2e1;
t235 = t135 * pkin(4);
t234 = t136 * pkin(4);
t232 = t137 * pkin(4);
t229 = pkin(4) * qJD(4);
t183 = -t234 / 0.2e1;
t6 = (t183 + t180) * t137;
t228 = qJD(1) * t6;
t225 = t51 * t135;
t175 = t225 / 0.2e1;
t28 = t47 * t210;
t213 = t135 * t136;
t40 = t51 * t213;
t181 = -t28 / 0.2e1 - t40 / 0.2e1;
t184 = -t134 / 0.2e1 - 0.1e1 / 0.2e1;
t8 = t136 * t175 + (pkin(4) * t184 + t238 * t50) * t137 + t181;
t227 = qJD(1) * t8;
t226 = t135 * t97;
t161 = -t139 + t235;
t86 = t161 * t138;
t224 = t86 * t135;
t205 = t138 * t139;
t100 = t137 * t205;
t92 = t135 * t98;
t54 = qJ(5) * t213 + t100 + t92;
t85 = t161 * t136;
t9 = t47 * t48 + t51 * t54 - t85 * t86;
t223 = t9 * qJD(1);
t27 = t47 * t212;
t11 = t212 * t50 - t27;
t222 = qJD(1) * t11;
t185 = pkin(4) * t207;
t12 = t86 * t185 + (-t47 + t50) * t51;
t221 = qJD(1) * t12;
t16 = t47 * t137 + t225;
t15 = t16 * t138;
t220 = qJD(1) * t15;
t219 = qJD(1) * t16;
t41 = -t134 * t211 - t136 * t58;
t218 = qJD(1) * t41;
t42 = -t134 * t206 - t136 * t59;
t217 = qJD(1) * t42;
t132 = t136 ^ 2;
t108 = t132 - t134;
t90 = t108 * t135;
t216 = qJD(1) * t90;
t208 = t137 * t132;
t91 = -t125 + t208;
t215 = qJD(1) * t91;
t10 = -t28 - t40 + (t135 * t54 + t137 * t48) * t138;
t214 = t10 * qJD(1);
t176 = t135 * t205;
t17 = -t58 * t138 + (t93 + t176) * t136;
t204 = t17 * qJD(1);
t18 = t59 * t138 + (-t100 + t92) * t136;
t203 = t18 * qJD(1);
t63 = (-t132 / 0.2e1 + t184) * t135;
t202 = t63 * qJD(1);
t65 = t125 / 0.2e1 + (0.1e1 / 0.2e1 + t132 / 0.2e1) * t137;
t201 = t65 * qJD(1);
t87 = t109 * t134;
t200 = t87 * qJD(1);
t89 = t109 * t138;
t199 = t89 * qJD(1);
t198 = qJD(1) * qJ(2);
t197 = qJD(2) * t136;
t196 = qJD(3) * t135;
t195 = qJD(3) * t137;
t194 = qJD(4) * t135;
t193 = qJD(4) * t137;
t192 = t108 * qJD(1);
t191 = t109 * qJD(3);
t190 = t136 * qJD(1);
t189 = t136 * qJD(3);
t188 = t138 * qJD(1);
t187 = t138 * qJD(3);
t186 = t138 * qJD(4);
t182 = t232 / 0.2e1;
t126 = -pkin(3) - t232;
t177 = t126 * t207;
t174 = qJ(2) * t190;
t173 = qJ(2) * t188;
t172 = t135 * t195;
t171 = t135 * t187;
t170 = t137 * t187;
t169 = t136 * t194;
t168 = t135 * t186;
t167 = t136 * t193;
t166 = t137 * t186;
t165 = t135 * t193;
t164 = t136 * t187;
t163 = t136 * t188;
t159 = -t236 * t97 + t246;
t157 = qJD(3) * t160;
t96 = t230 * t135;
t56 = -t135 * t96 - t137 * t97;
t140 = (t54 * t237 + t48 * t240 + t86 / 0.2e1) * t136 - t27 / 0.2e1;
t4 = -t85 * t236 + t226 / 0.2e1 + (t138 * t243 - t96 / 0.2e1) * t137 + t140;
t57 = (-0.1e1 + t109) * t138 * t136;
t155 = t4 * qJD(1) + t57 * qJD(2);
t1 = t248 + (t245 - t177 / 0.2e1 - t224 / 0.2e1) * pkin(4);
t19 = t126 * t235;
t154 = -qJD(1) * t1 + qJD(3) * t19;
t153 = (-qJD(4) - t190) * t138;
t152 = t233 / 0.2e1 + t231 / 0.2e1;
t151 = (t236 * t96 + t243) * t137;
t146 = t152 * t137;
t53 = -t93 / 0.2e1 - t146;
t150 = pkin(3) * t196 - qJD(1) * t53;
t145 = t152 * t135;
t52 = t92 / 0.2e1 + t145;
t149 = pkin(3) * t195 - qJD(1) * t52;
t79 = (t131 / 0.2e1 + t241) * t138;
t148 = -qJD(1) * t79 + t172;
t147 = t137 * t153;
t144 = qJD(1) * t125 * t135 + qJD(3) * t79;
t88 = t110 * t134;
t143 = qJD(1) * t88 + t157;
t14 = -t209 / 0.2e1 + t151 + (t234 / 0.2e1 + t159) * t135;
t64 = (0.1e1 / 0.2e1 + t241 - t131 / 0.2e1) * t136;
t141 = qJD(1) * t14 - qJD(2) * t64 + qJD(3) * t56;
t124 = t187 / 0.2e1;
t121 = t137 * t190;
t120 = t135 * t189;
t119 = t135 * t190;
t84 = (t190 + qJD(4) / 0.2e1) * t138;
t76 = (t137 * t188 + t196) * pkin(4);
t75 = t79 * qJD(4);
t68 = -t208 / 0.2e1 + t242 + t237;
t67 = (0.1e1 + t109) * t238;
t66 = t240 + (t132 + t134) * t135 / 0.2e1;
t33 = -t176 + t93 / 0.2e1 - t146;
t32 = -t100 - t92 / 0.2e1 + t145;
t29 = pkin(4) * t212;
t13 = t209 / 0.2e1 + t151 + (t183 + t159) * t135;
t7 = pkin(4) * t242 + t182 + (t237 * t50 + t175) * t136 + t181;
t5 = t136 * t182 + t137 * t180;
t3 = -t226 / 0.2e1 + t96 * t237 + (t51 * t237 + t85 / 0.2e1) * t138 + t140;
t2 = pkin(4) * t245 - t248 + (t177 + t224) * pkin(4) / 0.2e1;
t20 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t164, t108 * qJD(3), 0, 0, 0, qJ(2) * t187 + t197, -qJ(2) * t189 + qJD(2) * t138, -t133 * t164 - t134 * t165, -qJD(4) * t88 + t136 * t157, -qJD(3) * t91 - t136 * t168, qJD(3) * t90 - t136 * t166, t164, qJD(3) * t17 + qJD(4) * t42 + t137 * t197, -qJD(3) * t18 - qJD(4) * t41 - t135 * t197, -qJD(2) * t89 - qJD(3) * t10 - qJD(4) * t11 + qJD(5) * t87, qJD(2) * t16 + qJD(3) * t9 + qJD(4) * t12 - qJD(5) * t15; 0, 0, 0, 0, qJD(1), t198, 0, 0, 0, 0, 0, t190, t188, 0, 0, 0, 0, 0, qJD(4) * t68 + t121, qJD(4) * t66 - t119, -t199, qJD(3) * t3 + qJD(4) * t7 + t219; 0, 0, 0, 0, 0, 0, -t163, t192, -t189, -t187, 0, -t139 * t189 + t173, -t139 * t187 - t174, -t75 + (-t133 * t188 - t172) * t136, t136 * t142 - 0.2e1 * t138 * t165, t171 - t215, t170 + t216, t84, t204 + (t135 * t158 - t179) * qJD(3) + t33 * qJD(4), -t203 + (-pkin(7) * t207 + (pkin(3) * t137 + t211) * t136) * qJD(3) + t32 * qJD(4), -t214 + ((t136 * t96 + t54) * t137 + (-t136 * t97 - t48) * t135) * qJD(3) + t5 * qJD(4), t223 + t3 * qJD(2) + (-t126 * t85 + t48 * t96 - t54 * t97) * qJD(3) + t2 * qJD(4) + t13 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -t143, t135 * t153, t147, t124, qJD(2) * t68 + qJD(3) * t33 - qJD(4) * t59 + t217, qJD(2) * t66 + qJD(3) * t32 + qJD(4) * t58 - t218, pkin(4) * t168 + qJD(3) * t5 - t222, qJD(2) * t7 + qJD(3) * t2 - t229 * t51 + t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, qJD(3) * t13 - t220; 0, 0, 0, 0, -qJD(1), -t198, 0, 0, 0, 0, 0, -t190, -t188, 0, 0, 0, 0, 0, -qJD(4) * t65 - t121, -qJD(4) * t63 + t119, t199, qJD(3) * t4 + qJD(4) * t8 - t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t187, 0, 0, 0, 0, 0, -t137 * t189 - t168, t120 - t166, t89 * qJD(3), (t136 * t126 + t138 * t56) * qJD(3) - t29 * qJD(4) + t67 * qJD(5) + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 - t171 - t201, t169 - t170 - t202, 0, -pkin(4) * t167 - qJD(3) * t29 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * qJD(3); 0, 0, 0, 0, 0, 0, t163, -t192, 0, 0, 0, -t173, t174, t133 * t163 - t75, t147 * t249, t167 + t215, -t169 - t216, -t84, qJD(4) * t53 - t204, qJD(4) * t52 + t203, qJD(4) * t6 + t214, -qJD(2) * t4 - qJD(4) * t1 + qJD(5) * t14 - t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t64 - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t110 * qJD(4), 0, 0, 0, -pkin(3) * t194, -pkin(3) * t193, qJD(5) * t109, qJD(4) * t19 + qJD(5) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t142, t121 + t193, -t119 - t194, -t188 / 0.2e1, -pkin(7) * t193 - t150, pkin(7) * t194 - t149, -pkin(4) * t193 + t228, t229 * t97 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t143, (t135 * t188 - t195) * t136, t137 * t163 + t120, t124, qJD(2) * t65 - qJD(3) * t53 - t217, qJD(2) * t63 - qJD(3) * t52 + t218, -qJD(3) * t6 + t222, -qJD(2) * t8 + qJD(3) * t1 - qJD(5) * t185 - t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t202, 0, -t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t142, -t121, t119, t188 / 0.2e1, t150, t149, -t228, -qJD(5) * t235 - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, pkin(4) * t166 - qJD(3) * t14 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, pkin(4) * t194 - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t20;
