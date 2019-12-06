% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:32
% DurationCPUTime: 1.84s
% Computational Cost: add. (2859->219), mult. (7766->308), div. (0->0), fcn. (6071->8), ass. (0->151)
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t141 = sin(pkin(9));
t142 = cos(pkin(9));
t145 = sin(qJ(2));
t148 = cos(qJ(2));
t122 = -t141 * t145 + t142 * t148;
t112 = t122 * qJD(1);
t123 = t141 * t148 + t142 * t145;
t114 = t123 * qJD(1);
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t160 = t112 * t144 + t147 * t114;
t71 = t147 * t112 - t114 * t144;
t26 = t143 * t71 + t146 * t160;
t62 = t146 * t71;
t30 = -t143 * t160 + t62;
t220 = t26 * t30;
t138 = qJD(2) + qJD(4);
t137 = qJD(5) + t138;
t191 = t137 * t30;
t180 = qJD(5) * t143;
t113 = t123 * qJD(2);
t103 = qJD(1) * t113;
t179 = qJD(1) * qJD(2);
t174 = t148 * t179;
t175 = t145 * t179;
t104 = -t141 * t175 + t142 * t174;
t181 = qJD(4) * t147;
t182 = qJD(4) * t144;
t23 = -t144 * t103 + t147 * t104 + t112 * t181 - t114 * t182;
t24 = qJD(4) * t160 + t147 * t103 + t104 * t144;
t6 = qJD(5) * t62 - t143 * t24 + t146 * t23 - t160 * t180;
t216 = t6 - t191;
t152 = -qJD(5) * t26 - t143 * t23 - t146 * t24;
t192 = t137 * t26;
t212 = t152 + t192;
t217 = t26 ^ 2 - t30 ^ 2;
t199 = pkin(7) * t114;
t198 = -qJ(3) - pkin(6);
t131 = t198 * t148;
t128 = qJD(1) * t131;
t117 = t141 * t128;
t130 = t198 * t145;
t127 = qJD(1) * t130;
t194 = qJD(2) * pkin(2);
t121 = t127 + t194;
t73 = t142 * t121 + t117;
t49 = qJD(2) * pkin(3) - t199 + t73;
t200 = pkin(7) * t112;
t188 = t142 * t128;
t74 = t141 * t121 - t188;
t54 = t74 + t200;
t162 = -t144 * t49 - t147 * t54;
t222 = pkin(8) * t71;
t13 = -t162 + t222;
t178 = -pkin(2) * t148 - pkin(1);
t165 = t178 * qJD(1);
t129 = qJD(3) + t165;
t79 = -pkin(3) * t112 + t129;
t40 = -pkin(4) * t71 + t79;
t225 = t13 * t180 - t40 * t30;
t169 = qJD(2) * t198;
t109 = qJD(3) * t148 + t145 * t169;
t91 = t109 * qJD(1);
t110 = -qJD(3) * t145 + t148 * t169;
t92 = t110 * qJD(1);
t55 = -t141 * t91 + t142 * t92;
t38 = -pkin(7) * t104 + t55;
t56 = t141 * t92 + t142 * t91;
t39 = -pkin(7) * t103 + t56;
t155 = -(qJD(4) * t49 + t39) * t147 - t144 * t38 + t54 * t182;
t2 = -pkin(8) * t24 - t155;
t154 = qJD(4) * t162 - t144 * t39 + t147 * t38;
t3 = -pkin(8) * t23 + t154;
t213 = -t143 * t2 + t146 * t3 - t40 * t26;
t189 = t138 * t71;
t223 = t23 - t189;
t221 = t160 * t71;
t190 = t138 * t160;
t219 = -t24 + t190;
t218 = t160 ^ 2 - t71 ^ 2;
t215 = -t79 * t71 + t155;
t214 = (-t13 * t137 - t3) * t143 + t225;
t210 = -0.2e1 * t179;
t202 = pkin(4) * t160;
t209 = pkin(8) * t160;
t134 = pkin(2) * t142 + pkin(3);
t201 = pkin(2) * t141;
t108 = t134 * t144 + t147 * t201;
t77 = -t127 * t141 + t188;
t57 = t77 - t200;
t78 = t142 * t127 + t117;
t58 = t78 - t199;
t208 = -t108 * qJD(4) + t144 * t58 - t147 * t57;
t158 = t134 * t147 - t144 * t201;
t207 = -t158 * qJD(4) + t144 * t57 + t147 * t58;
t204 = qJD(5) - t137;
t203 = -t79 * t160 + t154;
t197 = t208 + t222;
t196 = t207 - t209;
t65 = t142 * t109 + t141 * t110;
t193 = t13 * t146;
t150 = qJD(1) ^ 2;
t187 = t148 * t150;
t149 = qJD(2) ^ 2;
t186 = t149 * t145;
t185 = t149 * t148;
t82 = t141 * t130 - t142 * t131;
t184 = t145 ^ 2 - t148 ^ 2;
t183 = qJD(1) * t145;
t136 = t145 * t194;
t171 = -t144 * t54 + t147 * t49;
t12 = t171 - t209;
t10 = pkin(4) * t138 + t12;
t176 = -pkin(4) * t137 - t10;
t133 = pkin(2) * t175;
t80 = pkin(3) * t103 + t133;
t87 = pkin(3) * t113 + t136;
t86 = pkin(2) * t183 + pkin(3) * t114;
t64 = -t109 * t141 + t142 * t110;
t166 = pkin(1) * t210;
t81 = t142 * t130 + t131 * t141;
t164 = -t10 * t143 - t193;
t159 = t147 * t122 - t123 * t144;
t76 = t122 * t144 + t123 * t147;
t35 = t143 * t76 - t146 * t159;
t36 = t143 * t159 + t146 * t76;
t66 = -pkin(7) * t123 + t81;
t67 = pkin(7) * t122 + t82;
t161 = -t144 * t66 - t147 * t67;
t94 = -pkin(3) * t122 + t178;
t116 = t122 * qJD(2);
t45 = -pkin(7) * t116 + t64;
t46 = -pkin(7) * t113 + t65;
t156 = t144 * t45 + t147 * t46 + t66 * t181 - t67 * t182;
t153 = qJD(4) * t161 - t144 * t46 + t147 * t45;
t107 = pkin(4) + t158;
t50 = -pkin(4) * t159 + t94;
t41 = t86 + t202;
t34 = qJD(4) * t76 + t147 * t113 + t116 * t144;
t33 = qJD(4) * t159 - t113 * t144 + t116 * t147;
t19 = pkin(4) * t34 + t87;
t18 = pkin(4) * t24 + t80;
t17 = pkin(8) * t159 - t161;
t16 = -pkin(8) * t76 - t144 * t67 + t147 * t66;
t9 = qJD(5) * t36 + t143 * t33 + t146 * t34;
t8 = -qJD(5) * t35 - t143 * t34 + t146 * t33;
t5 = -pkin(8) * t33 + t153;
t4 = -pkin(8) * t34 + t156;
t1 = [0, 0, 0, 0.2e1 * t145 * t174, t184 * t210, t185, -t186, 0, -pkin(6) * t185 + t145 * t166, pkin(6) * t186 + t148 * t166, -t103 * t82 - t104 * t81 + t112 * t65 - t113 * t74 - t114 * t64 - t116 * t73 + t122 * t56 - t123 * t55, t55 * t81 + t56 * t82 + t73 * t64 + t74 * t65 + (t129 + t165) * t136, t160 * t33 + t23 * t76, t159 * t23 - t160 * t34 - t24 * t76 + t33 * t71, t33 * t138, -t34 * t138, 0, t138 * t153 - t159 * t80 + t94 * t24 + t79 * t34 - t71 * t87, -t138 * t156 + t160 * t87 + t94 * t23 + t79 * t33 + t80 * t76, t26 * t8 + t36 * t6, t152 * t36 - t26 * t9 + t30 * t8 - t35 * t6, t8 * t137, -t9 * t137, 0, -t19 * t30 - t50 * t152 + t18 * t35 + t40 * t9 + (-t143 * t4 + t146 * t5 + (-t143 * t16 - t146 * t17) * qJD(5)) * t137, t19 * t26 + t50 * t6 + t18 * t36 + t40 * t8 - (t143 * t5 + t146 * t4 + (-t143 * t17 + t146 * t16) * qJD(5)) * t137; 0, 0, 0, -t145 * t187, t184 * t150, 0, 0, 0, t150 * pkin(1) * t145, pkin(1) * t187, (t74 + t77) * t114 + (t73 - t78) * t112 + (-t103 * t141 - t104 * t142) * pkin(2), -t73 * t77 - t74 * t78 + (-t129 * t183 + t141 * t56 + t142 * t55) * pkin(2), -t221, t218, t223, t219, 0, t208 * t138 + t71 * t86 + t203, t207 * t138 - t86 * t160 + t215, -t220, t217, t216, t212, 0, t41 * t30 + (t196 * t143 + t197 * t146) * t137 + ((-t107 * t143 - t108 * t146) * t137 + t164) * qJD(5) + t213, -t41 * t26 + (-t3 + (qJD(5) * t108 - t197) * t137) * t143 + (-qJD(5) * t10 - t2 + (-qJD(5) * t107 + t196) * t137) * t146 + t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 ^ 2 - t114 ^ 2, -t112 * t74 + t114 * t73 + t133, 0, 0, 0, 0, 0, t24 + t190, t23 + t189, 0, 0, 0, 0, 0, -t152 + t192, t6 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, t218, t223, t219, 0, -t138 * t162 + t203, t138 * t171 + t215, -t220, t217, t216, t212, 0, t30 * t202 - (-t12 * t143 - t193) * t137 + (t176 * t143 - t193) * qJD(5) + t213, -t26 * t202 + (t176 * qJD(5) + t12 * t137 - t2) * t146 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t217, t216, t212, 0, t204 * t164 + t213, (-t204 * t10 - t2) * t146 + t214;];
tauc_reg = t1;
