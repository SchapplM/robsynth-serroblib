% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:18
% EndTime: 2019-03-08 20:12:25
% DurationCPUTime: 2.23s
% Computational Cost: add. (3761->309), mult. (9673->424), div. (0->0), fcn. (7647->10), ass. (0->154)
t118 = sin(pkin(11));
t120 = cos(pkin(11));
t123 = sin(qJ(4));
t199 = cos(qJ(4));
t136 = -t123 * t118 + t199 * t120;
t193 = pkin(8) + qJ(3);
t106 = t193 * t118;
t107 = t193 * t120;
t71 = t199 * t106 + t123 * t107;
t48 = t136 * qJD(3) - qJD(4) * t71;
t119 = sin(pkin(6));
t126 = cos(qJ(2));
t178 = t119 * t126;
t131 = t136 * t178;
t74 = qJD(1) * t131;
t220 = -t48 + t74;
t124 = sin(qJ(2));
t172 = qJD(1) * t119;
t157 = t124 * t172;
t104 = qJD(2) * qJ(3) + t157;
t121 = cos(pkin(6));
t171 = qJD(1) * t121;
t110 = t120 * t171;
t187 = pkin(8) * qJD(2);
t68 = t110 + (-t104 - t187) * t118;
t80 = t120 * t104 + t118 * t171;
t69 = t120 * t187 + t80;
t33 = t123 * t68 + t199 * t69;
t219 = t33 * qJD(4);
t218 = t136 * qJD(2);
t102 = t199 * t118 + t123 * t120;
t95 = qJD(2) * t102;
t216 = t95 * qJD(4);
t89 = qJD(5) - t218;
t96 = t136 * qJD(4);
t129 = qJD(2) * t96;
t215 = qJD(5) * qJD(4) + t129;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t168 = qJD(5) * t125;
t169 = qJD(5) * t122;
t156 = t126 * t172;
t100 = (qJD(3) + t156) * qJD(2);
t209 = t136 * t100;
t210 = -t123 * t69 + t199 * t68;
t18 = qJD(4) * t210 + t209;
t29 = qJD(4) * pkin(9) + t33;
t113 = -t120 * pkin(3) - pkin(2);
t146 = qJD(3) - t156;
t88 = t113 * qJD(2) + t146;
t37 = -pkin(4) * t218 - t95 * pkin(9) + t88;
t108 = qJD(2) * t157;
t97 = t102 * qJD(4);
t87 = qJD(2) * t97;
t44 = t87 * pkin(4) - pkin(9) * t129 + t108;
t154 = t122 * t18 - t125 * t44 + t29 * t168 + t37 * t169;
t203 = t87 * pkin(5);
t2 = t154 - t203;
t12 = t122 * t37 + t125 * t29;
t7 = t89 * qJ(6) + t12;
t204 = t7 * t89;
t214 = -t2 + t204;
t132 = t102 * t178;
t72 = -t123 * t106 + t199 * t107;
t189 = -qJD(1) * t132 + t102 * qJD(3) + t72 * qJD(4);
t82 = t122 * t87;
t213 = -t89 * t168 - t82;
t60 = t97 * pkin(4) - t96 * pkin(9);
t61 = -pkin(4) * t136 - t102 * pkin(9) + t113;
t212 = -t61 * t168 + t72 * t169 + t220 * t125 + (t157 - t60) * t122;
t84 = t125 * t87;
t211 = -t89 * t169 + t84;
t28 = -qJD(4) * pkin(4) - t210;
t75 = -t125 * qJD(4) + t122 * t95;
t77 = t122 * qJD(4) + t125 * t95;
t15 = t75 * pkin(5) - t77 * qJ(6) + t28;
t205 = pkin(9) * t87;
t208 = t15 * t89 - t205;
t207 = t77 ^ 2;
t206 = t89 ^ 2;
t202 = t97 * qJ(6) - qJD(6) * t136 - t212;
t188 = t122 * t61 + t125 * t72;
t50 = t122 * t74 - t125 * t157;
t201 = -t97 * pkin(5) + t188 * qJD(5) + t122 * t48 - t125 * t60 - t50;
t147 = pkin(5) * t122 - qJ(6) * t125;
t148 = t125 * pkin(5) + t122 * qJ(6);
t200 = -t147 * t96 - (t148 * qJD(5) - qJD(6) * t125) * t102 - t189;
t198 = t12 * t89;
t197 = t75 * t218;
t196 = t77 * t75;
t155 = t77 * t89;
t195 = t77 * t95;
t194 = t95 * t75;
t192 = t122 * qJD(6) - t89 * t147 + t33;
t59 = t95 * pkin(4) - pkin(9) * t218;
t191 = t122 * t59 + t125 * t210;
t46 = t215 * t122 + t95 * t168;
t190 = -t122 * t46 - t75 * t168;
t186 = pkin(9) * qJD(5);
t185 = qJD(2) * pkin(2);
t184 = t122 * t89;
t183 = t122 * t96;
t182 = t125 * t96;
t45 = -t215 * t125 + t95 * t169;
t181 = t45 * t122;
t180 = t87 * qJ(6);
t179 = t119 * t124;
t127 = qJD(2) ^ 2;
t177 = t119 * t127;
t11 = -t122 * t29 + t125 * t37;
t174 = qJD(6) - t11;
t173 = t118 ^ 2 + t120 ^ 2;
t170 = qJD(2) * t124;
t166 = t89 * t186;
t163 = t124 * t177;
t162 = t126 * t177;
t158 = t119 * t170;
t19 = t102 * t100 + t219;
t152 = t173 * t100;
t3 = t46 * pkin(5) + t45 * qJ(6) - t77 * qJD(6) + t19;
t151 = -t3 - t166;
t138 = t122 * t44 + t125 * t18 + t37 * t168 - t29 * t169;
t1 = t89 * qJD(6) + t138 + t180;
t6 = -t89 * pkin(5) + t174;
t150 = t89 * t6 + t1;
t149 = -t122 * t7 + t125 * t6;
t145 = t118 * (-t118 * t104 + t110) - t120 * t80;
t144 = -t125 * t218 * t89 - t213;
t143 = t184 * t218 + t211;
t142 = t15 * t77 + t154;
t91 = -t118 * t179 + t121 * t120;
t92 = t121 * t118 + t120 * t179;
t140 = -t123 * t92 + t199 * t91;
t53 = t123 * t91 + t199 * t92;
t41 = t122 * t53 + t125 * t178;
t42 = -t122 * t178 + t125 * t53;
t139 = t89 * t28 - t205;
t128 = qJD(2) * t131;
t30 = t140 * qJD(4) + t128;
t10 = t42 * qJD(5) + t122 * t30 - t125 * t158;
t31 = qJD(2) * t132 + t53 * qJD(4);
t135 = -t10 * t89 - t140 * t46 + t31 * t75 - t41 * t87;
t9 = -t41 * qJD(5) + t122 * t158 + t125 * t30;
t130 = t140 * t45 + t31 * t77 - t42 * t87 - t9 * t89;
t105 = -pkin(4) - t148;
t103 = t146 - t185;
t40 = t77 * pkin(5) + t75 * qJ(6);
t34 = t147 * t102 + t71;
t25 = pkin(5) * t136 + t122 * t72 - t125 * t61;
t24 = -qJ(6) * t136 + t188;
t22 = t75 * t89 - t45;
t14 = -t95 * pkin(5) + t122 * t210 - t125 * t59;
t13 = t95 * qJ(6) + t191;
t4 = [0, 0, -t163, -t162, -t120 * t163, t118 * t163, t173 * t162 (-t118 * t91 + t120 * t92) * t100 + (t103 * t124 + (-t145 - t157) * t126) * t119 * qJD(2), 0, 0, 0, 0, 0, -t31 * qJD(4) + (-t126 * t87 - t170 * t218) * t119, t95 * t158 + (-t128 - t30) * qJD(4), 0, 0, 0, 0, 0, t135, t130, t135, t10 * t77 - t41 * t45 - t42 * t46 - t9 * t75, -t130, t1 * t42 + t6 * t10 - t140 * t3 + t15 * t31 + t2 * t41 + t7 * t9; 0, 0, 0, 0, 0, 0, qJD(2) * t146 * t173 + t152, -t145 * qJD(3) + qJ(3) * t152 + (t145 * t126 + (-t103 - t185) * t124) * t172, t102 * t129 + t95 * t96, -t102 * t87 + t129 * t136 + t218 * t96 - t95 * t97, t96 * qJD(4), -t97 * qJD(4), 0, -t189 * qJD(4) + t113 * t87 + t88 * t97, t220 * qJD(4) + t102 * t108 + t113 * t129 - t95 * t157 + t88 * t96, t77 * t182 + (-t45 * t125 - t77 * t169) * t102 (-t122 * t77 - t125 * t75) * t96 + (t181 - t125 * t46 + (t122 * t75 - t125 * t77) * qJD(5)) * t102, t211 * t102 + t136 * t45 + t89 * t182 + t77 * t97, t213 * t102 + t136 * t46 - t89 * t183 - t75 * t97, -t136 * t87 + t89 * t97, t154 * t136 + t11 * t97 + t71 * t46 + t50 * t89 + t189 * t75 + ((-qJD(5) * t72 + t60) * t89 + t61 * t87 + t28 * qJD(5) * t102) * t125 + ((-qJD(5) * t61 - t48) * t89 - t72 * t87 + t19 * t102 + t28 * t96) * t122, -t188 * t87 + t138 * t136 - t12 * t97 - t71 * t45 + t28 * t182 + t212 * t89 + t189 * t77 + (t19 * t125 - t169 * t28) * t102, t15 * t183 + t2 * t136 - t25 * t87 + t34 * t46 - t6 * t97 - t201 * t89 - t200 * t75 + (t3 * t122 + t15 * t168) * t102, -t24 * t46 - t25 * t45 + t149 * t96 + t201 * t77 - t202 * t75 + (-t1 * t122 + t2 * t125 + (-t122 * t6 - t125 * t7) * qJD(5)) * t102, -t15 * t182 - t1 * t136 + t24 * t87 + t34 * t45 + t7 * t97 + t202 * t89 + t200 * t77 + (-t3 * t125 + t15 * t169) * t102, t1 * t24 - t200 * t15 + t2 * t25 + t201 * t6 + t202 * t7 + t3 * t34; 0, 0, 0, 0, 0, 0, -t173 * t127, t145 * qJD(2) + t108, 0, 0, 0, 0, 0, 0.2e1 * t216, 0.2e1 * t218 * qJD(4), 0, 0, 0, 0, 0, t143 - t194, -t125 * t206 - t195 - t82, -t184 * t89 - t194 + t84 (t45 + t197) * t125 + t122 * t155 + t190, t144 + t195, t150 * t122 + t214 * t125 - t15 * t95; 0, 0, 0, 0, 0, 0, 0, 0, -t95 * t218, -t218 ^ 2 + t95 ^ 2, 0, 0, 0, -t88 * t95 - t19 + t219, -t218 * t88 - t209, t125 * t155 - t181 (-t45 + t197) * t125 - t77 * t184 + t190, t144 - t195, t143 + t194, -t89 * t95, -pkin(4) * t46 - t11 * t95 - t33 * t75 + (-t19 + (-t59 - t186) * t89) * t125 + (t210 * t89 + t139) * t122, pkin(4) * t45 + t191 * t89 + t12 * t95 - t33 * t77 + (t19 + t166) * t122 + t139 * t125, t105 * t46 + t122 * t208 + t151 * t125 + t14 * t89 - t192 * t75 + t6 * t95, t13 * t75 - t14 * t77 + ((qJD(5) * t77 - t46) * pkin(9) + t150) * t125 + ((qJD(5) * t75 - t45) * pkin(9) - t214) * t122, t105 * t45 + t151 * t122 - t125 * t208 - t13 * t89 + t192 * t77 - t7 * t95, t3 * t105 - t7 * t13 - t6 * t14 - t192 * t15 + (qJD(5) * t149 + t1 * t125 + t2 * t122) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t75 ^ 2 + t207, t22, -t46 + t155, t87, -t28 * t77 - t154 + t198, t11 * t89 + t28 * t75 - t138, -t40 * t75 - t142 + t198 + 0.2e1 * t203, pkin(5) * t45 - t46 * qJ(6) + (-t12 + t7) * t77 + (t6 - t174) * t75, 0.2e1 * t180 - t15 * t75 + t40 * t77 + (0.2e1 * qJD(6) - t11) * t89 + t138, -t2 * pkin(5) + t1 * qJ(6) - t6 * t12 - t15 * t40 + t174 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196 - t216, t22, -t206 - t207, t142 - t203 - t204;];
tauc_reg  = t4;
