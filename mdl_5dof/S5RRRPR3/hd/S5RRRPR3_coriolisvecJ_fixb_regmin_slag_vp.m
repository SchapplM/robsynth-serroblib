% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:35
% EndTime: 2021-01-15 22:59:43
% DurationCPUTime: 1.40s
% Computational Cost: add. (2216->225), mult. (3910->308), div. (0->0), fcn. (2720->8), ass. (0->160)
t150 = sin(pkin(9));
t151 = cos(pkin(9));
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t119 = t150 * t156 + t151 * t153;
t147 = qJD(1) + qJD(2);
t100 = t119 * t147;
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t187 = qJD(5) * t152;
t196 = t151 * t156;
t179 = t147 * t196;
t197 = t150 * t153;
t98 = t147 * t197 - t179;
t85 = t155 * t98;
t111 = t119 * qJD(3);
t91 = t147 * t111;
t190 = qJD(3) * t153;
t177 = t147 * t190;
t127 = t150 * t177;
t189 = qJD(3) * t156;
t176 = t147 * t189;
t92 = t151 * t176 - t127;
t10 = -qJD(5) * t85 - t100 * t187 - t152 * t91 + t155 * t92;
t146 = qJD(3) + qJD(5);
t42 = -t152 * t100 - t85;
t199 = t42 * t146;
t225 = t10 - t199;
t164 = -t155 * t100 + t152 * t98;
t159 = t164 * qJD(5) - t152 * t92 - t155 * t91;
t200 = t164 * t146;
t224 = t159 - t200;
t223 = qJ(4) + pkin(7);
t222 = t164 * t42;
t221 = t164 ^ 2 - t42 ^ 2;
t217 = t98 * pkin(8);
t154 = sin(qJ(2));
t203 = pkin(1) * qJD(1);
t182 = t154 * t203;
t171 = t223 * t147 + t182;
t88 = t171 * t156;
t201 = t151 * t88;
t87 = t171 * t153;
t78 = qJD(3) * pkin(3) - t87;
t34 = t150 * t78 + t201;
t18 = t34 - t217;
t139 = -t156 * pkin(3) - pkin(2);
t157 = cos(qJ(2));
t181 = t157 * t203;
t97 = t139 * t147 + qJD(4) - t181;
t49 = t98 * pkin(4) + t97;
t220 = t18 * t187 - t49 * t42;
t202 = pkin(1) * qJD(2);
t180 = qJD(1) * t202;
t169 = t157 * t180;
t161 = qJD(4) * t147 + t169;
t163 = qJD(3) * t171;
t47 = -t153 * t163 + t161 * t156;
t48 = -t161 * t153 - t156 * t163;
t12 = -t150 * t47 + t151 * t48;
t4 = -t92 * pkin(8) + t12;
t13 = t150 * t48 + t151 * t47;
t5 = -t91 * pkin(8) + t13;
t219 = -t152 * t5 + t155 * t4 + t49 * t164;
t142 = t156 * qJD(4);
t172 = qJD(3) * t223;
t108 = -t153 * t172 + t142;
t109 = -t153 * qJD(4) - t156 * t172;
t207 = -t150 * t108 + t151 * t109 + t119 * t181;
t118 = -t196 + t197;
t206 = t151 * t108 + t150 * t109 + t118 * t181;
t218 = qJD(5) - t146;
t112 = t118 * qJD(3);
t60 = -t152 * t118 + t155 * t119;
t26 = t60 * qJD(5) + t155 * t111 - t152 * t112;
t134 = t154 * t180;
t110 = pkin(3) * t177 + t134;
t54 = t91 * pkin(4) + t110;
t59 = t155 * t118 + t152 * t119;
t216 = t49 * t26 + t54 * t59;
t25 = -t59 * qJD(5) - t152 * t111 - t155 * t112;
t215 = t49 * t25 + t54 * t60;
t214 = pkin(3) * t150;
t213 = pkin(3) * t153;
t212 = t100 * pkin(8);
t211 = t112 * pkin(8);
t210 = t119 * pkin(8);
t209 = t157 * pkin(1);
t205 = t110 * t118 + t97 * t111;
t204 = t110 * t119 - t97 * t112;
t137 = t154 * pkin(1) + pkin(7);
t193 = -qJ(4) - t137;
t170 = qJD(3) * t193;
t183 = t157 * t202;
t67 = t153 * t170 + t156 * t183 + t142;
t68 = (-qJD(4) - t183) * t153 + t156 * t170;
t30 = t150 * t68 + t151 * t67;
t72 = t150 * t88;
t36 = -t151 * t87 - t72;
t198 = t147 * t153;
t195 = t154 * t156;
t158 = qJD(3) ^ 2;
t194 = t158 * t153;
t143 = t158 * t156;
t116 = t193 * t153;
t144 = t156 * qJ(4);
t117 = t156 * t137 + t144;
t58 = t150 * t116 + t151 * t117;
t126 = -t147 * pkin(2) - t181;
t192 = t126 * t189 + t153 * t134;
t131 = t223 * t153;
t132 = t156 * pkin(7) + t144;
t77 = -t150 * t131 + t151 * t132;
t191 = t153 ^ 2 - t156 ^ 2;
t188 = qJD(3) * t157;
t186 = -qJD(1) - t147;
t185 = -qJD(2) + t147;
t184 = pkin(3) * t198;
t141 = t154 * t202;
t140 = pkin(3) * t190;
t175 = t153 * t188;
t33 = t151 * t78 - t72;
t174 = -t34 * t111 + t33 * t112 - t13 * t118 - t12 * t119;
t84 = t111 * pkin(4) + t140;
t29 = -t150 * t67 + t151 * t68;
t35 = t150 * t87 - t201;
t57 = t151 * t116 - t150 * t117;
t76 = -t151 * t131 - t150 * t132;
t168 = t84 - t182;
t115 = t118 * pkin(8);
t167 = qJD(5) * (-t115 + t77) - t211 - t207;
t107 = t111 * pkin(8);
t166 = -qJD(5) * (t76 - t210) + t107 - t206;
t17 = qJD(3) * pkin(4) - t212 + t33;
t165 = -t152 * t17 - t155 * t18;
t93 = t118 * pkin(4) + t139;
t162 = -t126 * t147 - t169;
t160 = -t154 * t198 + t156 * t188;
t145 = t147 ^ 2;
t138 = -pkin(2) - t209;
t136 = t151 * pkin(3) + pkin(4);
t128 = t139 - t209;
t124 = 0.2e1 * t153 * t176;
t123 = t141 + t140;
t113 = t126 * t190;
t102 = -0.2e1 * t191 * t147 * qJD(3);
t83 = t93 - t209;
t71 = t141 + t84;
t66 = t100 * pkin(4) + t184;
t38 = -t115 + t58;
t37 = t57 - t210;
t22 = t36 - t212;
t21 = t35 + t217;
t20 = t26 * t146;
t19 = t25 * t146;
t15 = -t107 + t30;
t14 = t29 + t211;
t2 = t10 * t60 - t164 * t25;
t1 = -t10 * t59 + t159 * t60 + t164 * t26 + t25 * t42;
t3 = [0, 0, 0, 0, -t147 * t141 - t134, t186 * t183, t124, t102, t143, -t194, 0, t138 * t177 - t137 * t143 + t113 + (t186 * t195 - t175) * t202, t137 * t194 + t138 * t176 - t160 * t202 + t192, t29 * qJD(3) + t123 * t98 + t128 * t91 + t205, -t30 * qJD(3) + t123 * t100 + t128 * t92 + t204, -t29 * t100 - t30 * t98 - t57 * t92 - t58 * t91 + t174, t110 * t128 + t12 * t57 + t97 * t123 + t13 * t58 + t33 * t29 + t34 * t30, t2, t1, t19, -t20, 0, -t71 * t42 - t83 * t159 + (t155 * t14 - t152 * t15 + (-t152 * t37 - t155 * t38) * qJD(5)) * t146 + t216, -t71 * t164 + t83 * t10 - (t152 * t14 + t155 * t15 + (-t152 * t38 + t155 * t37) * qJD(5)) * t146 + t215; 0, 0, 0, 0, t147 * t182 - t134, t185 * t181, t124, t102, t143, -t194, 0, -pkin(2) * t177 - pkin(7) * t143 + t113 + (t185 * t195 + t175) * t203, -pkin(2) * t176 + pkin(7) * t194 + t160 * t203 + t192, -t98 * t182 + t139 * t91 + (t98 * t213 + t207) * qJD(3) + t205, -t100 * t182 + t139 * t92 + (t100 * t213 - t206) * qJD(3) + t204, -t207 * t100 - t206 * t98 - t76 * t92 - t77 * t91 + t174, t110 * t139 + t12 * t76 + t13 * t77 + (-t182 + t140) * t97 + t206 * t34 + t207 * t33, t2, t1, t19, -t20, 0, -t93 * t159 - t168 * t42 + (t166 * t152 - t167 * t155) * t146 + t216, t93 * t10 - t168 * t164 + (t167 * t152 + t166 * t155) * t146 + t215; 0, 0, 0, 0, 0, 0, -t153 * t145 * t156, t191 * t145, 0, 0, 0, t162 * t153, t162 * t156, -t35 * qJD(3) - t97 * t100 - t98 * t184 + t12, t36 * qJD(3) - t100 * t184 + t97 * t98 - t13, (-t33 + t36) * t98 + (t34 + t35) * t100 + (-t150 * t91 - t151 * t92) * pkin(3), -t33 * t35 - t34 * t36 + (t12 * t151 + t13 * t150 - t97 * t198) * pkin(3), t222, t221, t225, t224, 0, t66 * t42 - (-t152 * t22 + t155 * t21) * t146 + ((-t136 * t152 - t155 * t214) * t146 + t165) * qJD(5) + t219, -t155 * t5 - t152 * t4 + t66 * t164 + (t152 * t21 + t155 * t22) * t146 + (-(t136 * t155 - t152 * t214) * t146 - t155 * t17) * qJD(5) + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100 * qJD(3), -t127 + (-t98 + t179) * qJD(3), -t100 ^ 2 - t98 ^ 2, t33 * t100 + t34 * t98 + t110, 0, 0, 0, 0, 0, -t159 - t200, t10 + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t221, t225, t224, 0, t218 * t165 + t219, (-t18 * t146 - t4) * t152 + (-t218 * t17 - t5) * t155 + t220;];
tauc_reg = t3;
