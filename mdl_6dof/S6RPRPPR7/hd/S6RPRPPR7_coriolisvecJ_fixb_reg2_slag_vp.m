% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:26
% EndTime: 2019-03-09 02:57:32
% DurationCPUTime: 2.18s
% Computational Cost: add. (4161->328), mult. (8935->410), div. (0->0), fcn. (5938->6), ass. (0->185)
t120 = cos(qJ(6));
t117 = cos(pkin(9));
t121 = cos(qJ(3));
t176 = qJD(1) * t121;
t164 = t117 * t176;
t119 = sin(qJ(3));
t193 = sin(pkin(9));
t161 = t193 * t119;
t99 = qJD(1) * t161;
t85 = -t99 + t164;
t78 = qJD(6) + t85;
t157 = t120 * t78;
t118 = sin(qJ(6));
t89 = t117 * t119 + t193 * t121;
t82 = t89 * qJD(1);
t57 = t120 * qJD(3) + t118 * t82;
t234 = t57 * t157;
t87 = t89 * qJD(3);
t70 = qJD(1) * t87;
t113 = qJD(1) * qJD(2);
t169 = qJD(1) * qJD(3);
t163 = t121 * t169;
t92 = pkin(3) * t163 + t113;
t150 = t70 * qJ(5) + t92;
t131 = -t85 * qJD(5) + t150;
t220 = pkin(4) + pkin(8);
t96 = qJD(3) * t99;
t69 = -t117 * t163 + t96;
t11 = -t220 * t69 + t131;
t177 = qJD(1) * t119;
t122 = -pkin(1) - pkin(7);
t97 = t122 * qJD(1) + qJD(2);
t76 = -qJ(4) * t177 + t119 * t97;
t63 = t193 * t76;
t77 = -qJ(4) * t176 + t121 * t97;
t66 = qJD(3) * pkin(3) + t77;
t35 = t117 * t66 - t63;
t149 = qJD(5) - t35;
t215 = t85 * pkin(5);
t16 = -t220 * qJD(3) + t149 + t215;
t93 = pkin(3) * t177 + qJD(1) * qJ(2) + qJD(4);
t133 = -t85 * qJ(5) + t93;
t20 = t220 * t82 + t133;
t135 = t118 * t20 - t120 * t16;
t170 = t121 * qJD(4);
t175 = qJD(3) * t119;
t53 = -t97 * t175 + (qJ(4) * t175 - t170) * qJD(1);
t171 = t119 * qJD(4);
t174 = qJD(3) * t121;
t54 = t97 * t174 + (-qJ(4) * t174 - t171) * qJD(1);
t23 = -t117 * t53 + t193 * t54;
t14 = -t70 * pkin(5) + t23;
t1 = -t135 * qJD(6) + t120 * t11 + t118 * t14;
t148 = t135 * t78 + t1;
t6 = t118 * t16 + t120 * t20;
t2 = -qJD(6) * t6 - t118 * t11 + t120 * t14;
t233 = t6 * t78 + t2;
t222 = t82 ^ 2;
t81 = t85 ^ 2;
t232 = -t222 - t81;
t231 = -t222 + t81;
t172 = t118 * qJD(3);
t55 = -t120 * t82 + t172;
t162 = t55 * t78;
t173 = qJD(6) * t120;
t33 = qJD(6) * t172 + t118 * t69 - t82 * t173;
t230 = t33 - t162;
t34 = t57 * qJD(6) + t120 * t69;
t229 = t57 * t78 - t34;
t90 = t117 * t121 - t161;
t47 = t90 * t70;
t137 = -t87 * t85 - t47;
t43 = t117 * t77 - t63;
t181 = -qJD(5) + t43;
t190 = qJD(3) * t87;
t228 = qJD(1) * t82 + t190;
t24 = t117 * t54 + t193 * t53;
t197 = t117 * t76;
t36 = t193 * t66 + t197;
t84 = qJD(3) * t161 - t117 * t174;
t227 = -t24 * t89 + t35 * t87 + t36 * t84;
t19 = -qJD(3) * qJD(5) - t24;
t30 = -qJD(3) * pkin(4) + t149;
t32 = -qJD(3) * qJ(5) - t36;
t226 = t19 * t89 - t30 * t87 - t32 * t84;
t105 = -t117 * pkin(3) - pkin(4);
t101 = -pkin(8) + t105;
t216 = t82 * pkin(5);
t18 = -t32 - t216;
t225 = -t101 * t70 + t18 * t78;
t224 = qJD(3) * (t85 + t164) - t96;
t223 = qJD(3) * (-t85 + t164) - t96;
t221 = 0.2e1 * t113;
t217 = t69 * pkin(4);
t182 = qJ(4) - t122;
t94 = t182 * t119;
t95 = t182 * t121;
t51 = t117 * t95 - t193 * t94;
t213 = t23 * t51;
t212 = t23 * t90;
t37 = t82 * pkin(4) + t133;
t206 = t37 * t85;
t205 = t57 * t55;
t203 = t57 * t82;
t202 = t70 * t89;
t201 = t82 * t55;
t200 = t82 * t85;
t196 = t118 * t70;
t195 = t120 * t33;
t62 = t120 * t70;
t194 = t34 * t118;
t191 = qJD(3) * t82;
t123 = qJD(3) ^ 2;
t189 = t123 * t119;
t188 = t123 * t121;
t124 = qJD(1) ^ 2;
t187 = t124 * qJ(2);
t186 = t124 * t121;
t71 = t182 * t175 - t170;
t72 = -qJD(3) * t95 - t171;
t38 = -t117 * t71 + t193 * t72;
t185 = t38 * qJD(3);
t39 = t117 * t72 + t193 * t71;
t184 = t39 * qJD(3);
t183 = t84 * qJD(3);
t106 = t119 * pkin(3) + qJ(2);
t180 = t215 - t181;
t179 = t119 ^ 2 - t121 ^ 2;
t178 = -t123 - t124;
t98 = pkin(3) * t174 + qJD(2);
t168 = 0.2e1 * qJD(1);
t109 = pkin(3) * t176;
t167 = qJD(6) * t118 * t89;
t166 = t89 * t173;
t165 = t119 * t186;
t160 = t82 * qJ(5) + t109;
t159 = t118 * t78;
t156 = qJD(6) * t90 + qJD(1);
t155 = t119 * t163;
t154 = -t1 * t90 + t6 * t87;
t153 = -t135 * t87 - t2 * t90;
t12 = t69 * pkin(5) - t19;
t152 = -qJD(6) * t101 * t78 + t12;
t151 = -t90 * qJ(5) + t106;
t145 = -t118 * t135 - t120 * t6;
t144 = t12 * t89 - t18 * t84;
t143 = -t89 * t33 - t84 * t57;
t142 = -t33 * t90 - t57 * t87;
t141 = -t89 * t34 + t84 * t55;
t140 = t34 * t90 - t55 * t87;
t139 = -t89 * t69 - t84 * t82;
t138 = -t78 * t84 - t202;
t42 = t193 * t77 + t197;
t31 = t220 * t89 + t151;
t40 = t90 * pkin(5) + t51;
t10 = t118 * t40 + t120 * t31;
t9 = -t118 * t31 + t120 * t40;
t134 = qJD(1) * t85 - t183;
t132 = -t78 * t159 - t62;
t52 = -t117 * t94 - t193 * t95;
t130 = t87 * qJ(5) - t90 * qJD(5) + t98;
t129 = t42 * qJD(3) - t23;
t128 = -t157 * t78 + t196;
t127 = t90 * t69 + t87 * t82 + t85 * t84 + t202;
t126 = -t137 - t139;
t125 = t38 * t85 - t39 * t82 - t51 * t70 + t52 * t69 + t212;
t111 = qJ(2) * t221;
t103 = t193 * pkin(3) + qJ(5);
t46 = t89 * pkin(4) + t151;
t45 = -t191 + t70;
t44 = t85 * pkin(4) + t160;
t41 = -t89 * pkin(5) + t52;
t29 = t120 * t34;
t28 = -t84 * pkin(4) + t130;
t27 = t220 * t85 + t160;
t25 = t42 - t216;
t22 = t84 * pkin(5) + t39;
t21 = -t87 * pkin(5) + t38;
t17 = t131 - t217;
t15 = -t220 * t84 + t130;
t8 = t118 * t25 + t120 * t27;
t7 = -t118 * t27 + t120 * t25;
t4 = -t10 * qJD(6) - t118 * t15 + t120 * t21;
t3 = t9 * qJD(6) + t118 * t21 + t120 * t15;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t111, -0.2e1 * t155, 0.2e1 * t179 * t169, -t189, 0.2e1 * t155, -t188, 0, -t122 * t189 + (qJ(2) * t174 + qJD(2) * t119) * t168, -t122 * t188 + (-qJ(2) * t175 + qJD(2) * t121) * t168, 0, t111, t137, t127, -t190, t139, t183, 0, -t106 * t69 + t98 * t82 - t93 * t84 + t92 * t89 - t185, -t106 * t70 + t98 * t85 - t93 * t87 + t92 * t90 - t184, t125 + t227, t92 * t106 + t24 * t52 - t35 * t38 + t36 * t39 + t93 * t98 + t213, 0, t190, -t183, t137, t127, t139, t125 + t226, -t17 * t89 - t28 * t82 + t37 * t84 + t46 * t69 + t185, -t17 * t90 - t28 * t85 + t37 * t87 + t46 * t70 + t184, t17 * t46 - t19 * t52 + t37 * t28 + t30 * t38 - t32 * t39 + t213, t118 * t143 + t166 * t57 (t118 * t55 - t120 * t57) * t84 + (-t194 - t195 + (-t118 * t57 - t120 * t55) * qJD(6)) * t89, t118 * t138 + t166 * t78 + t142, t120 * t141 + t167 * t55, t120 * t138 - t167 * t78 - t140, -t78 * t87 - t47, -t120 * t144 + t167 * t18 + t22 * t55 + t41 * t34 + t4 * t78 - t9 * t70 - t153, t10 * t70 + t118 * t144 + t166 * t18 + t22 * t57 - t3 * t78 - t41 * t33 + t154, -t10 * t34 - t3 * t55 + t9 * t33 - t4 * t57 + t145 * t84 + (t1 * t120 - t2 * t118 + (-t118 * t6 + t120 * t135) * qJD(6)) * t89, t1 * t10 + t12 * t41 - t135 * t4 + t18 * t22 + t2 * t9 + t6 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t187, 0, 0, 0, 0, 0, 0, t178 * t119, t178 * t121, 0, -t187, 0, 0, 0, 0, 0, 0, -t228, -t134, t126, -t93 * qJD(1) - t212 - t227, 0, 0, 0, 0, 0, 0, t126, t228, t134, -t37 * qJD(1) - t212 - t226, 0, 0, 0, 0, 0, 0, t90 * t62 + (t118 * t156 + t120 * t87) * t78 - t141, -t90 * t196 + (-t118 * t87 + t120 * t156) * t78 + t143 (t156 * t55 + t142) * t120 + (-t156 * t57 + t140) * t118 (-t156 * t6 + t153) * t120 + (-t135 * t156 + t154) * t118 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t179 * t124, 0, -t165, 0, 0, -qJ(2) * t186, t119 * t187, 0, 0, t200, t231, 0, -t200, -t223, 0, -t109 * t82 - t93 * t85 + t129, t43 * qJD(3) - t109 * t85 + t93 * t82 - t24 (t36 - t42) * t85 + (-t35 + t43) * t82 + (t117 * t70 + t193 * t69) * pkin(3), t35 * t42 - t36 * t43 + (-t117 * t23 - t176 * t93 + t193 * t24) * pkin(3), 0, t45, t223, t200, t231, -t200, t103 * t69 - t105 * t70 + (-t32 - t42) * t85 + (t30 + t181) * t82, t44 * t82 - t129 + t206, -t37 * t82 + t44 * t85 + (0.2e1 * qJD(5) - t43) * qJD(3) + t24, -t19 * t103 + t23 * t105 + t181 * t32 - t30 * t42 - t37 * t44, -t159 * t57 - t195, -t29 - t234 + (t33 + t162) * t118, t132 + t203, t157 * t55 + t194, t128 - t201, t78 * t82, t103 * t34 + t152 * t118 + t225 * t120 - t135 * t82 + t180 * t55 - t7 * t78, -t103 * t33 - t225 * t118 + t152 * t120 + t180 * t57 - t6 * t82 + t8 * t78, t8 * t55 + t7 * t57 + (t101 * t33 - t6 * t85 - t2 + (-t101 * t55 - t6) * qJD(6)) * t120 + (-t101 * t34 - t135 * t85 - t1 + (t101 * t57 - t135) * qJD(6)) * t118, t12 * t103 + t135 * t7 - t6 * t8 + t180 * t18 + (-qJD(6) * t145 + t1 * t118 + t2 * t120) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -0.2e1 * t191, t232, t35 * t85 + t36 * t82 + t92, 0, 0, 0, 0, 0, 0, t232, -t224, t191 + t70, -t217 - t32 * t82 + (-qJD(5) - t30) * t85 + t150, 0, 0, 0, 0, 0, 0, t128 + t201, t118 * t78 ^ 2 + t203 + t62, -t230 * t118 + t234 - t29, -t233 * t118 + t148 * t120 + t18 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t200, -t81 - t123, t32 * qJD(3) + t206 + t23, 0, 0, 0, 0, 0, 0, -qJD(3) * t55 + t132, -qJD(3) * t57 + t128, t229 * t118 + t230 * t120, -t18 * qJD(3) + t148 * t118 + t233 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t55 ^ 2 + t57 ^ 2, -t230, -t205, t229, -t70, -t18 * t57 + t233, t18 * t55 - t148, 0, 0;];
tauc_reg  = t5;
