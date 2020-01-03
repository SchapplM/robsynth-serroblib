% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR5
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:58
% EndTime: 2020-01-03 12:04:04
% DurationCPUTime: 2.02s
% Computational Cost: add. (4605->254), mult. (7862->330), div. (0->0), fcn. (5666->8), ass. (0->160)
t160 = cos(qJ(2));
t207 = pkin(1) * qJD(1);
t187 = t160 * t207;
t166 = qJD(3) - t187;
t154 = cos(pkin(9));
t159 = cos(qJ(4));
t195 = t159 * t154;
t153 = sin(pkin(9));
t156 = sin(qJ(4));
t196 = t156 * t153;
t130 = -t195 + t196;
t136 = (-pkin(7) - qJ(3)) * t153;
t147 = t154 * pkin(7);
t137 = t154 * qJ(3) + t147;
t192 = qJD(4) * t159;
t209 = t136 * t192 + qJD(3) * t195 + (-qJD(3) * t153 - qJD(4) * t137) * t156 + t130 * t187;
t131 = t159 * t153 + t156 * t154;
t79 = t156 * t136 + t159 * t137;
t208 = -t79 * qJD(4) - t166 * t131;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t152 = qJD(1) + qJD(2);
t100 = t131 * t152;
t185 = t152 * t196;
t98 = -t152 * t195 + t185;
t165 = -t158 * t100 + t155 * t98;
t183 = t154 * t192;
t134 = t152 * t183;
t184 = qJD(4) * t196;
t93 = t152 * t184 - t134;
t120 = t131 * qJD(4);
t94 = t152 * t120;
t161 = t165 * qJD(5) + t155 * t93 - t158 * t94;
t151 = qJD(4) + qJD(5);
t201 = t165 * t151;
t231 = t161 - t201;
t190 = qJD(5) * t158;
t191 = qJD(5) * t155;
t163 = -t100 * t191 - t155 * t94 - t158 * t93 - t98 * t190;
t63 = -t155 * t100 - t158 * t98;
t200 = t63 * t151;
t230 = t163 - t200;
t215 = t165 ^ 2;
t216 = t63 ^ 2;
t229 = t215 - t216;
t146 = -t154 * pkin(3) - pkin(2);
t97 = t146 * t152 + t166;
t67 = t98 * pkin(4) + t97;
t228 = t67 * t63;
t115 = t120 * pkin(8);
t227 = -t115 + t209;
t119 = -t183 + t184;
t220 = t119 * pkin(8);
t226 = t220 + t208;
t214 = t63 * t165;
t225 = t67 * t165;
t206 = pkin(1) * qJD(2);
t186 = qJD(1) * t206;
t129 = t152 * qJD(3) + t160 * t186;
t193 = t153 ^ 2 + t154 ^ 2;
t224 = t193 * t129;
t157 = sin(qJ(2));
t223 = t154 * t157;
t145 = t157 * pkin(1) + qJ(3);
t121 = (-pkin(7) - t145) * t153;
t122 = t154 * t145 + t147;
t73 = t156 * t121 + t159 * t122;
t162 = t131 * t129;
t188 = t157 * t207;
t133 = t152 * qJ(3) + t188;
t180 = pkin(7) * t152 + t133;
t86 = t180 * t153;
t87 = t180 * t154;
t55 = -t156 * t86 + t159 * t87;
t31 = -t55 * qJD(4) - t162;
t21 = t93 * pkin(8) + t31;
t41 = -t98 * pkin(8) + t55;
t178 = t155 * t21 - t41 * t191;
t197 = t129 * t195 - t86 * t192;
t30 = (-qJD(4) * t87 - t129 * t153) * t156 + t197;
t20 = -t94 * pkin(8) + t30;
t203 = t156 * t87;
t54 = -t159 * t86 - t203;
t40 = -t100 * pkin(8) + t54;
t37 = qJD(4) * pkin(4) + t40;
t4 = (qJD(5) * t37 + t20) * t158 + t178;
t222 = t100 ^ 2;
t221 = pkin(4) * t100;
t219 = t120 * pkin(4);
t218 = t131 * pkin(8);
t217 = t160 * pkin(1);
t78 = t159 * t136 - t156 * t137;
t68 = t78 - t218;
t124 = t130 * pkin(8);
t69 = -t124 + t79;
t32 = -t155 * t69 + t158 * t68;
t213 = t32 * qJD(5) + t226 * t155 + t227 * t158;
t33 = t155 * t68 + t158 * t69;
t212 = -t33 * qJD(5) - t227 * t155 + t226 * t158;
t75 = -t155 * t130 + t158 * t131;
t43 = t75 * qJD(5) - t155 * t119 + t158 * t120;
t74 = t158 * t130 + t155 * t131;
t143 = t157 * t186;
t76 = t94 * pkin(4) + t143;
t211 = t67 * t43 + t76 * t74;
t42 = t158 * t119 + t155 * t120 + t130 * t190 + t131 * t191;
t210 = -t67 * t42 + t76 * t75;
t205 = t100 * t98;
t204 = t155 * t41;
t202 = t158 * t41;
t199 = t97 * t120 + t130 * t143;
t198 = -t97 * t119 + t131 * t143;
t189 = t157 * t206;
t14 = t158 * t37 - t204;
t15 = t155 * t37 + t202;
t179 = -t155 * t20 + t158 * t21;
t5 = -t15 * qJD(5) + t179;
t182 = t14 * t42 - t15 * t43 - t4 * t74 - t5 * t75;
t181 = -pkin(4) * t151 - t37;
t176 = t54 * t119 - t55 * t120 - t30 * t130 - t31 * t131;
t174 = t193 * t160;
t140 = t160 * t206 + qJD(3);
t172 = t193 * t140;
t72 = t159 * t121 - t156 * t122;
t171 = t193 * qJD(3);
t170 = t152 * t189;
t169 = t152 * t188;
t168 = (-qJD(2) + t152) * t207;
t167 = (-qJD(1) - t152) * t206;
t58 = t72 - t218;
t59 = -t124 + t73;
t24 = -t155 * t59 + t158 * t58;
t25 = t155 * t58 + t158 * t59;
t95 = t130 * pkin(4) + t146;
t164 = -t188 + t219;
t48 = t121 * t192 + t140 * t195 + (-qJD(4) * t122 - t140 * t153) * t156;
t49 = -t73 * qJD(4) - t131 * t140;
t138 = t153 * t143;
t135 = t146 - t217;
t132 = -t152 * pkin(2) + t166;
t113 = t120 * qJD(4);
t112 = t119 * qJD(4);
t96 = t98 ^ 2;
t88 = t189 + t219;
t84 = t95 - t217;
t45 = t98 * t120 + t94 * t130;
t44 = -t100 * t119 - t93 * t131;
t39 = t43 * t151;
t38 = t42 * t151;
t35 = t49 + t220;
t34 = -t115 + t48;
t18 = -t100 * t120 + t119 * t98 + t93 * t130 - t131 * t94;
t17 = t158 * t40 - t204;
t16 = -t155 * t40 - t202;
t9 = -t161 * t74 - t43 * t63;
t8 = t163 * t75 + t165 * t42;
t7 = -t25 * qJD(5) - t155 * t34 + t158 * t35;
t6 = t24 * qJD(5) + t155 * t35 + t158 * t34;
t1 = t161 * t75 - t163 * t74 + t165 * t43 - t42 * t63;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 - t170, t160 * t167, 0, 0, 0, 0, 0, 0, 0, 0, t167 * t223, t153 * t170 + t138, t152 * t172 + t224, t133 * t172 + t145 * t224 + (t132 + (-pkin(2) - t217) * qJD(1)) * t189, t44, t18, -t112, t45, -t113, 0, t49 * qJD(4) + t135 * t94 + t189 * t98 + t199, -t48 * qJD(4) + t100 * t189 - t135 * t93 + t198, -t49 * t100 - t48 * t98 + t72 * t93 - t73 * t94 + t176, t30 * t73 + t31 * t72 + t55 * t48 + t54 * t49 + (qJD(1) * t135 + t97) * t189, t8, t1, -t38, t9, -t39, 0, t7 * t151 - t161 * t84 - t63 * t88 + t211, -t6 * t151 + t163 * t84 - t165 * t88 + t210, t161 * t25 - t163 * t24 + t165 * t7 + t6 * t63 + t182, t14 * t7 + t15 * t6 + t5 * t24 + t4 * t25 + t67 * t88 + t76 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 + t169, t160 * t168, 0, 0, 0, 0, 0, 0, 0, 0, t168 * t223, -t153 * t169 + t138, (-t174 * t207 + t171) * t152 + t224, t133 * t171 + qJ(3) * t224 + ((-pkin(2) * qJD(2) - t132) * t157 - t133 * t174) * t207, t44, t18, -t112, t45, -t113, 0, t208 * qJD(4) + t146 * t94 - t98 * t188 + t199, -t209 * qJD(4) - t100 * t188 - t146 * t93 + t198, -t208 * t100 - t209 * t98 + t78 * t93 - t79 * t94 + t176, t30 * t79 + t31 * t78 + t209 * t55 + t208 * t54 + (qJD(2) * t146 - t97) * t188, t8, t1, -t38, t9, -t39, 0, t212 * t151 - t161 * t95 - t164 * t63 + t211, -t213 * t151 + t163 * t95 - t164 * t165 + t210, t161 * t33 - t163 * t32 + t165 * t212 + t213 * t63 + t182, t212 * t14 + t213 * t15 + t164 * t67 + t5 * t32 + t4 * t33 + t76 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 * t152 ^ 2, -t193 * t152 * t133 + t143, 0, 0, 0, 0, 0, 0, 0.2e1 * t100 * qJD(4), t134 + (-t98 - t185) * qJD(4), -t96 - t222, t54 * t100 + t55 * t98 + t143, 0, 0, 0, 0, 0, 0, -t161 - t201, t163 + t200, -t215 - t216, -t14 * t165 - t15 * t63 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t96 + t222, t134 + (t98 - t185) * qJD(4), -t205, 0, 0, -t97 * t100 - t162, t129 * t196 + t97 * t98 + (t54 + t203) * qJD(4) - t197, 0, 0, t214, t229, t230, -t214, t231, 0, t63 * t221 - t16 * t151 + t225 + (t155 * t181 - t202) * qJD(5) + t179, t165 * t221 + t17 * t151 - t228 + (qJD(5) * t181 - t20) * t158 - t178, t14 * t63 - t15 * t165 - t16 * t165 - t17 * t63 + (t155 * t161 - t158 * t163 + (-t155 * t165 + t158 * t63) * qJD(5)) * pkin(4), -t14 * t16 - t15 * t17 + (-t100 * t67 + t155 * t4 + t158 * t5 + (-t14 * t155 + t15 * t158) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t229, t230, -t214, t231, 0, t15 * t151 + t225 + t5, t14 * t151 - t228 - t4, 0, 0;];
tauc_reg = t2;
