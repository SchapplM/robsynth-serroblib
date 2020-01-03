% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:40
% DurationCPUTime: 1.71s
% Computational Cost: add. (3020->274), mult. (7457->381), div. (0->0), fcn. (5133->6), ass. (0->160)
t207 = cos(qJ(3));
t161 = t207 * qJD(3);
t153 = pkin(2) * t161;
t128 = sin(qJ(2));
t208 = -pkin(7) - pkin(6);
t111 = t208 * t128;
t104 = qJD(1) * t111;
t130 = cos(qJ(2));
t112 = t208 * t130;
t106 = qJD(1) * t112;
t127 = sin(qJ(3));
t94 = t127 * t106;
t67 = t207 * t104 + t94;
t213 = -t153 + t67;
t171 = qJD(1) * qJD(2);
t212 = -0.2e1 * t171;
t129 = cos(qJ(4));
t126 = sin(qJ(4));
t162 = qJD(1) * t207;
t175 = qJD(1) * t128;
t91 = t127 * t175 - t130 * t162;
t187 = t126 * t91;
t211 = -qJ(5) * t187 + t129 * qJD(5);
t210 = t207 * t111 + t127 * t112;
t170 = qJD(2) + qJD(3);
t140 = -t127 * t128 + t207 * t130;
t138 = t140 * qJD(3);
t70 = qJD(2) * t140 + t138;
t134 = t70 * qJD(1);
t181 = t127 * t130;
t93 = -qJD(1) * t181 - t128 * t162;
t139 = -t126 * t170 + t129 * t93;
t36 = -qJD(4) * t139 + t126 * t134;
t209 = t139 ^ 2;
t206 = t129 * pkin(4);
t191 = qJD(2) * pkin(2);
t96 = t104 + t191;
t64 = t207 * t96 + t94;
t51 = -t170 * pkin(3) - t64;
t205 = t51 * t70;
t154 = t129 * t170;
t74 = -t126 * t93 - t154;
t88 = qJD(4) + t91;
t204 = t74 * t88;
t203 = t139 * t88;
t202 = t88 * t93;
t201 = t93 * t91;
t200 = -qJ(5) - pkin(8);
t121 = -pkin(2) * t130 - pkin(1);
t110 = t121 * qJD(1);
t49 = pkin(3) * t91 + pkin(8) * t93 + t110;
t95 = t207 * t106;
t65 = t127 * t96 - t95;
t52 = t170 * pkin(8) + t65;
t21 = -t126 * t52 + t129 * t49;
t14 = qJ(5) * t139 + t21;
t13 = pkin(4) * t88 + t14;
t199 = t13 - t14;
t61 = -pkin(3) * t93 + pkin(8) * t91;
t198 = t126 * t61 + t129 * t64;
t50 = pkin(2) * t175 + t61;
t197 = t126 * t50 + t129 * t67;
t119 = t127 * pkin(2) + pkin(8);
t177 = -qJ(5) - t119;
t155 = qJD(4) * t177;
t196 = t126 * t155 + t129 * t153 - t197 + t211;
t123 = t129 * qJ(5);
t148 = -t93 * pkin(4) + t91 * t123;
t48 = t129 * t50;
t195 = t129 * t155 - t148 - t48 + (-qJD(5) + t213) * t126;
t102 = t207 * t128 + t181;
t63 = -pkin(3) * t140 - pkin(8) * t102 + t121;
t79 = t127 * t111 - t207 * t112;
t72 = t129 * t79;
t194 = t126 * t63 + t72;
t158 = qJD(4) * t200;
t193 = t126 * t158 - t198 + t211;
t159 = -t126 * t64 + t129 * t61;
t192 = -t126 * qJD(5) + t129 * t158 - t148 - t159;
t173 = qJD(4) * t126;
t35 = -qJD(4) * t154 - t129 * t134 - t93 * t173;
t190 = t126 * t35;
t71 = t170 * t102;
t60 = t71 * qJD(1);
t189 = t126 * t60;
t188 = t126 * t70;
t186 = t129 * t60;
t185 = t129 * t139;
t157 = t129 * t88;
t184 = t129 * t91;
t174 = qJD(3) * t127;
t166 = qJD(2) * t208;
t151 = qJD(1) * t166;
t97 = t128 * t151;
t98 = t130 * t151;
t30 = -t106 * t161 + t127 * t97 + t96 * t174 - t207 * t98;
t183 = t30 * t129;
t182 = t102 * t126;
t132 = qJD(1) ^ 2;
t180 = t130 * t132;
t131 = qJD(2) ^ 2;
t179 = t131 * t128;
t178 = t131 * t130;
t176 = t128 ^ 2 - t130 ^ 2;
t172 = qJD(4) * t129;
t167 = t128 * t191;
t33 = pkin(3) * t71 - pkin(8) * t70 + t167;
t105 = t128 * t166;
t107 = t130 * t166;
t38 = t210 * qJD(3) + t207 * t105 + t127 * t107;
t169 = t126 * t33 + t129 * t38 + t63 * t172;
t45 = t51 * t173;
t164 = t102 * t172;
t163 = t21 * t93 + t45;
t160 = t128 * t171;
t156 = pkin(1) * t212;
t22 = t126 * t49 + t129 * t52;
t152 = t30 * t126 + t51 * t172 - t22 * t93;
t120 = -t207 * pkin(2) - pkin(3);
t66 = t127 * t104 - t95;
t150 = pkin(2) * t174 - t66;
t149 = (t173 + t187) * pkin(4);
t147 = -t119 * t60 + t51 * t91;
t15 = -qJ(5) * t74 + t22;
t146 = -t126 * t15 - t129 * t13;
t145 = -qJ(5) * t70 - qJD(5) * t102;
t9 = pkin(4) * t36 + t30;
t144 = t110 * t93 - t30;
t143 = -t88 * t172 - t189;
t142 = t88 * t173 - t186;
t26 = t60 * pkin(3) + (-pkin(8) * t138 + (t128 * pkin(2) - pkin(8) * t140) * qJD(2)) * qJD(1);
t29 = t106 * t174 + t127 * t98 + t96 * t161 + t207 * t97;
t141 = t126 * t26 + t129 * t29 + t49 * t172 - t52 * t173;
t24 = t129 * t26;
t137 = -t22 * qJD(4) - t126 * t29 + t24;
t1 = pkin(4) * t60 + qJ(5) * t35 + qJD(5) * t139 + t137;
t3 = -qJ(5) * t36 - qJD(5) * t74 + t141;
t136 = t146 * qJD(4) - t1 * t126 + t3 * t129 - t13 * t184 - t15 * t187;
t135 = t110 * t91 - t29;
t39 = t79 * qJD(3) + t127 * t105 - t207 * t107;
t109 = pkin(8) * t129 + t123;
t108 = t200 * t126;
t100 = t119 * t129 + t123;
t99 = t177 * t126;
t73 = t74 ^ 2;
t58 = t129 * t63;
t42 = -t91 ^ 2 + t93 ^ 2;
t41 = (-qJD(1) * t102 - t93) * t170;
t40 = t91 * t170 + t134;
t34 = t74 * pkin(4) + qJD(5) + t51;
t32 = t129 * t33;
t25 = -qJ(5) * t182 + t194;
t19 = -pkin(4) * t140 - t102 * t123 - t126 * t79 + t58;
t11 = -t139 * t93 + t88 * t157 + t189;
t10 = -t88 ^ 2 * t126 - t74 * t93 + t186;
t8 = -t139 * t157 - t190;
t6 = -qJ(5) * t164 + (-qJD(4) * t79 + t145) * t126 + t169;
t5 = (-t35 - t204) * t129 + (-t36 + t203) * t126;
t4 = pkin(4) * t71 - t126 * t38 + t32 + t145 * t129 + (-t72 + (qJ(5) * t102 - t63) * t126) * qJD(4);
t2 = [0, 0, 0, 0.2e1 * t130 * t160, t176 * t212, t178, -t179, 0, -pkin(6) * t178 + t128 * t156, pkin(6) * t179 + t130 * t156, t102 * t134 - t93 * t70, -t102 * t60 + t134 * t140 - t70 * t91 + t93 * t71, t70 * t170, -t71 * t170, 0, t121 * t60 + t110 * t71 - t39 * t170 + (-qJD(1) * t140 + t91) * t167, pkin(2) * t102 * t160 + t110 * t70 + t121 * t134 - t93 * t167 - t38 * t170, -t70 * t185 + (-t129 * t35 + t139 * t173) * t102, (t126 * t139 - t129 * t74) * t70 + (t190 - t129 * t36 + (t126 * t74 + t185) * qJD(4)) * t102, -t102 * t142 - t139 * t71 + t140 * t35 + t70 * t157, t102 * t143 + t140 * t36 - t88 * t188 - t71 * t74, -t140 * t60 + t71 * t88, (-t172 * t79 + t32) * t88 + t58 * t60 - (-t172 * t52 + t24) * t140 + t21 * t71 + t39 * t74 - t210 * t36 + t51 * t164 + ((-qJD(4) * t63 - t38) * t88 - t79 * t60 - (-qJD(4) * t49 - t29) * t140 + t30 * t102 + t205) * t126, -(-t173 * t79 + t169) * t88 - t194 * t60 + t141 * t140 - t22 * t71 - t39 * t139 + t210 * t35 + t129 * t205 + (-t45 + t183) * t102, t19 * t35 - t25 * t36 + t4 * t139 - t6 * t74 + t146 * t70 + (-t1 * t129 - t126 * t3 + (t126 * t13 - t129 * t15) * qJD(4)) * t102, t3 * t25 + t15 * t6 + t1 * t19 + t13 * t4 + t9 * (pkin(4) * t182 - t210) + t34 * ((t164 + t188) * pkin(4) + t39); 0, 0, 0, -t128 * t180, t176 * t132, 0, 0, 0, t132 * pkin(1) * t128, pkin(1) * t180, -t201, t42, t40, t41, 0, t66 * t170 + (-t170 * t174 - t91 * t175) * pkin(2) + t144, t67 * t170 + (-t170 * t161 + t93 * t175) * pkin(2) + t135, t8, t5, t11, t10, t202, t120 * t36 - t48 * t88 + t150 * t74 + (-qJD(4) * t119 * t88 - t30) * t129 + (t213 * t88 + t147) * t126 + t163, -t120 * t35 + (t119 * t173 + t197) * t88 - t150 * t139 + (-t153 * t88 + t147) * t129 + t152, -t100 * t36 + t139 * t195 - t196 * t74 + t35 * t99 + t136, t3 * t100 + t1 * t99 + t9 * (t120 - t206) + (t95 + (pkin(2) * qJD(3) - t104) * t127 + t149) * t34 + t196 * t15 + t195 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, t42, t40, t41, 0, t65 * t170 + t144, t64 * t170 + t135, t8, t5, t11, t10, t202, -pkin(3) * t36 + pkin(8) * t143 - t159 * t88 + t51 * t187 - t65 * t74 + t163 - t183, pkin(3) * t35 + t142 * pkin(8) + t139 * t65 + t51 * t184 + t198 * t88 + t152, t108 * t35 - t109 * t36 + t139 * t192 - t193 * t74 + t136, t3 * t109 + t1 * t108 + t9 * (-pkin(3) - t206) + (t149 - t65) * t34 + t193 * t15 + t192 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139 * t74, -t73 + t209, -t35 + t204, -t36 - t203, t60, t139 * t51 + t22 * t88 + t137, t21 * t88 + t51 * t74 - t141, pkin(4) * t35 - t199 * t74, t199 * t15 + (t139 * t34 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 - t209, -t13 * t139 + t15 * t74 + t9;];
tauc_reg = t2;
