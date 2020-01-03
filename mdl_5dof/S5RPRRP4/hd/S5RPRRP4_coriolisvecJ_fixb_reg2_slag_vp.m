% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:30
% DurationCPUTime: 1.61s
% Computational Cost: add. (3063->252), mult. (8137->338), div. (0->0), fcn. (5649->6), ass. (0->172)
t178 = qJD(3) + qJD(4);
t132 = cos(qJ(4));
t129 = cos(pkin(8));
t181 = t129 * qJD(1);
t117 = -qJD(3) + t181;
t131 = sin(qJ(3));
t198 = qJ(2) * t131;
t170 = t129 * t198;
t133 = cos(qJ(3));
t128 = sin(pkin(8));
t214 = pkin(7) * t128;
t176 = t133 * t214;
t141 = -t170 - t176;
t104 = -t129 * pkin(2) - t128 * pkin(6) - pkin(1);
t92 = t104 * qJD(1) + qJD(2);
t88 = t133 * t92;
t61 = t141 * qJD(1) + t88;
t45 = -t117 * pkin(3) + t61;
t130 = sin(qJ(4));
t188 = qJD(1) * t128;
t166 = t131 * t188;
t168 = qJ(2) * t181;
t74 = t131 * t92 + t133 * t168;
t62 = -pkin(7) * t166 + t74;
t54 = t130 * t62;
t23 = t132 * t45 - t54;
t165 = t133 * t188;
t152 = t132 * t165;
t154 = t130 * t166;
t85 = t152 - t154;
t77 = t85 * qJ(5);
t14 = t23 - t77;
t101 = t130 * t133 + t132 * t131;
t89 = t101 * t128;
t216 = t131 * t133;
t215 = t85 ^ 2;
t213 = t132 * pkin(3);
t142 = qJD(1) * t101;
t82 = t128 * t142;
t212 = t85 * t82;
t93 = pkin(3) * t166 + qJ(2) * t188;
t211 = t93 * t85;
t111 = -qJD(4) + t117;
t13 = -t111 * pkin(4) + t14;
t210 = t13 - t14;
t182 = qJD(4) * t132;
t174 = pkin(3) * t182;
t194 = t132 * t133;
t169 = t128 * t194;
t147 = t178 * t169;
t191 = t178 * t154;
t44 = qJD(1) * t147 - t191;
t209 = -t130 * pkin(3) * t44 - t82 * t174;
t27 = t132 * t61 - t54;
t98 = t133 * t104;
t65 = -t176 + t98 + (-pkin(3) - t198) * t129;
t196 = t128 * t131;
t177 = pkin(7) * t196;
t197 = qJ(2) * t133;
t116 = t129 * t197;
t79 = t131 * t104 + t116;
t75 = -t177 + t79;
t31 = t130 * t65 + t132 * t75;
t195 = t130 * t131;
t100 = t194 - t195;
t208 = (t178 - t181) * t100;
t72 = t178 * t101;
t207 = t129 * t142 - t72;
t56 = t132 * t62;
t50 = t178 * t89;
t43 = qJD(1) * t50;
t206 = t43 * qJ(5);
t153 = t131 * t168;
t73 = t88 - t153;
t205 = t73 * t117;
t204 = t74 * t117;
t203 = t82 * qJ(5);
t202 = t82 * t111;
t201 = t85 * t111;
t180 = qJD(1) * qJD(2);
t163 = t129 * t180;
t184 = qJD(3) * t133;
t200 = t133 * t163 + t92 * t184;
t186 = qJD(2) * t133;
t199 = t104 * t184 + t129 * t186;
t124 = t128 ^ 2;
t134 = qJD(1) ^ 2;
t193 = t134 * t124;
t161 = -pkin(4) * t82 - qJD(5);
t59 = -t161 + t93;
t192 = qJD(5) + t59;
t156 = pkin(3) * t165;
t91 = qJD(3) * t156 + t128 * t180;
t95 = (pkin(3) * t184 + qJD(2)) * t128;
t99 = pkin(3) * t196 + t128 * qJ(2);
t190 = t129 ^ 2 + t124;
t189 = t131 ^ 2 - t133 ^ 2;
t187 = qJD(2) * t131;
t185 = qJD(3) * t131;
t183 = qJD(4) * t130;
t179 = qJD(1) * qJD(3);
t175 = pkin(3) * t183;
t173 = 0.2e1 * t180;
t172 = t92 * t185;
t171 = t128 * t195;
t167 = qJ(2) * t185;
t164 = t129 * t187;
t162 = t124 * t179;
t140 = t141 * qJD(3);
t39 = qJD(1) * t140 + t200;
t40 = -t172 + (-t164 + (-t116 + t177) * qJD(3)) * qJD(1);
t160 = -t130 * t39 + t132 * t40;
t26 = -t130 * t61 - t56;
t30 = -t130 * t75 + t132 * t65;
t159 = -t130 * t40 - t132 * t39 - t45 * t182 + t62 * t183;
t158 = t190 * t134;
t157 = 0.2e1 * qJD(1) * t124;
t155 = t193 * t216;
t32 = pkin(4) * t44 + t91;
t151 = (-qJD(3) - t117) * t188;
t24 = t130 * t45 + t56;
t150 = t162 * t216;
t149 = t190 * t173;
t148 = t93 * t82 + t159;
t146 = qJ(5) * t44 + t159;
t57 = t140 + t199;
t58 = -t164 + (-t116 + (-t104 + t214) * t131) * qJD(3);
t9 = t130 * t58 + t132 * t57 + t65 * t182 - t75 * t183;
t145 = qJD(3) * t128 * (t117 + t181);
t144 = qJ(2) * t184 + t187;
t143 = t129 * t179 + t193;
t139 = -t117 ^ 2 - t193;
t138 = t192 * t82 + t146;
t8 = -t24 * qJD(4) + t160;
t10 = -t31 * qJD(4) - t130 * t57 + t132 * t58;
t137 = t8 + t206;
t136 = (-t56 + (pkin(3) * t111 - t45) * t130) * qJD(4) + t160;
t135 = t72 * t188;
t121 = pkin(4) + t213;
t96 = t111 * t174;
t90 = t169 - t171;
t80 = t82 ^ 2;
t78 = t98 - t170;
t70 = t85 * pkin(4) + t156;
t68 = -t79 * qJD(3) - t164;
t67 = -t129 * t167 + t199;
t66 = pkin(4) * t89 + t99;
t53 = -t144 * t181 - t172;
t52 = -qJD(3) * t153 + t200;
t51 = -t178 * t171 + t147;
t34 = pkin(4) * t51 + t95;
t33 = -t80 + t215;
t29 = -t178 * t152 + t191 - t201;
t28 = -t135 - t202;
t25 = -qJ(5) * t89 + t31;
t22 = -t129 * pkin(4) - t90 * qJ(5) + t30;
t21 = -t207 * t111 - t82 * t188;
t20 = t208 * t111 - t85 * t188;
t19 = t51 * t111 + t44 * t129;
t18 = t50 * t111 + t43 * t129;
t17 = -t77 + t27;
t16 = t26 + t203;
t15 = t24 - t203;
t12 = t44 * t89 + t51 * t82;
t11 = -t43 * t90 - t50 * t85;
t6 = t50 * qJ(5) - t90 * qJD(5) + t10;
t5 = -t51 * qJ(5) - t89 * qJD(5) + t9;
t4 = t43 * t89 - t44 * t90 + t50 * t82 - t51 * t85;
t3 = -t85 * qJD(5) + t137;
t2 = -qJD(5) * t82 - t146;
t1 = t100 * t43 - t101 * t44 - t207 * t85 - t208 * t82;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, qJ(2) * t149, -0.2e1 * t150, 0.2e1 * t189 * t162, t131 * t145, 0.2e1 * t150, t133 * t145, 0, -t68 * t117 - t53 * t129 + t144 * t157, t67 * t117 + t52 * t129 + (-t167 + t186) * t157, (-t52 * t131 - t53 * t133 + (t131 * t73 - t133 * t74) * qJD(3) + (-t131 * t67 - t133 * t68 + (t131 * t78 - t133 * t79) * qJD(3)) * qJD(1)) * t128, qJ(2) * t124 * t173 + t52 * t79 + t53 * t78 + t67 * t74 + t68 * t73, t11, t4, t18, t12, t19, 0, -t10 * t111 - t8 * t129 + t99 * t44 + t93 * t51 + t95 * t82 + t91 * t89, t9 * t111 - t129 * t159 - t99 * t43 - t93 * t50 + t95 * t85 + t91 * t90, -t10 * t85 + t159 * t89 + t23 * t50 - t24 * t51 + t30 * t43 - t31 * t44 - t8 * t90 - t82 * t9, t10 * t23 - t159 * t31 + t24 * t9 + t30 * t8 + t91 * t99 + t93 * t95, t11, t4, t18, t12, t19, 0, -t6 * t111 - t3 * t129 + t32 * t89 + t34 * t82 + t66 * t44 + t59 * t51, t5 * t111 + t2 * t129 + t32 * t90 + t34 * t85 - t66 * t43 - t59 * t50, t13 * t50 - t15 * t51 - t2 * t89 + t22 * t43 - t25 * t44 - t3 * t90 - t5 * t82 - t6 * t85, t13 * t6 + t15 * t5 + t2 * t25 + t22 * t3 + t32 * t66 + t34 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158, -qJ(2) * t158, 0, 0, 0, 0, 0, 0, t139 * t131, t139 * t133, 0, -qJ(2) * t193 + (t53 - t204) * t133 + (t52 + t205) * t131, 0, 0, 0, 0, 0, 0, t21, t20, t1, t8 * t100 - t101 * t159 - t93 * t188 + t207 * t23 + t208 * t24, 0, 0, 0, 0, 0, 0, t21, t20, t1, t3 * t100 + t2 * t101 + t207 * t13 + t208 * t15 - t59 * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t189 * t193, t131 * t151, -t155, t133 * t151, 0, -t204 + (-qJD(3) * t92 - t163) * t131 - t143 * t197, t143 * t198 - t200 - t205, 0, 0, t212, t33, t28, -t212, t29, 0, t26 * t111 - t82 * t156 + t136 - t211, -t27 * t111 - t156 * t85 + t148 + t96, t43 * t213 + (-t23 + t27) * t82 + (t24 + t26 + t175) * t85 + t209, -t23 * t26 - t24 * t27 + (-t93 * t165 - t130 * t159 + t132 * t8 + (-t130 * t23 + t132 * t24) * qJD(4)) * pkin(3), t212, t33, t28, -t212, t29, 0, t16 * t111 - t192 * t85 - t70 * t82 + t136 + t206, -t111 * t17 - t70 * t85 + t138 + t96, t121 * t43 + (-t13 + t17) * t82 + (t15 + t16 + t175) * t85 + t209, t3 * t121 - t13 * t16 - t15 * t17 - t59 * t70 + (t130 * t2 + (-t13 * t130 + t132 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t33, t28, -t212, t29, 0, -t24 * t111 - t211 + t8, -t111 * t23 + t148, 0, 0, t212, t33, t28, -t212, t29, 0, -t15 * t111 + (t161 - t59) * t85 + t137, -pkin(4) * t215 - t111 * t14 + t138, pkin(4) * t43 - t210 * t82, t210 * t15 + (-t59 * t85 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 - t201, -t135 + t202, -t80 - t215, t13 * t85 + t15 * t82 + t32;];
tauc_reg = t7;
