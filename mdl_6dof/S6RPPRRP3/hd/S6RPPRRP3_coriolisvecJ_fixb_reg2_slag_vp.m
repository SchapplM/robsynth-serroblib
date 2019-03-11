% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:48
% EndTime: 2019-03-09 02:03:55
% DurationCPUTime: 2.17s
% Computational Cost: add. (2979->317), mult. (5992->409), div. (0->0), fcn. (3373->6), ass. (0->173)
t95 = cos(qJ(5));
t167 = qJD(5) * t95;
t93 = sin(qJ(5));
t168 = qJD(5) * t93;
t96 = cos(qJ(4));
t163 = t96 * qJD(2);
t82 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t70 = t82 * qJD(1) + qJD(3);
t94 = sin(qJ(4));
t51 = t94 * t70 + t163;
t44 = qJD(4) * pkin(8) + t51;
t218 = -t94 * qJD(2) + t70 * t96;
t45 = t218 * qJD(4);
t84 = sin(pkin(9)) * pkin(1) + qJ(3);
t67 = pkin(4) * t94 - pkin(8) * t96 + t84;
t54 = t67 * qJD(1);
t133 = pkin(4) * t96 + pkin(8) * t94;
t71 = t133 * qJD(4) + qJD(3);
t59 = t71 * qJD(1);
t141 = t44 * t167 + t54 * t168 + t93 * t45 - t95 * t59;
t17 = t44 * t95 + t54 * t93;
t173 = qJD(1) * t94;
t83 = qJD(5) + t173;
t111 = t17 * t83 - t141;
t164 = t95 * qJD(4);
t187 = t93 * t94;
t172 = qJD(1) * t96;
t73 = t93 * t172 - t164;
t227 = ((t73 + t164) * t96 - t83 * t187) * qJD(1) - t83 * t168;
t137 = qJD(5) * t94 + qJD(1);
t152 = t95 * t172;
t183 = t95 * t96;
t188 = t83 * t93;
t166 = qJD(5) * t96;
t147 = t93 * t166;
t109 = t94 * t164 + t147;
t48 = t109 * qJD(1) - qJD(5) * t164;
t171 = qJD(4) * t93;
t75 = t152 + t171;
t226 = ((-t75 + t152) * t94 + t83 * t183) * qJD(4) - t137 * t188 - t48 * t96;
t119 = t75 * t83;
t120 = t73 * t83;
t161 = qJD(1) * qJD(4);
t142 = t93 * t161;
t80 = t94 * t142;
t49 = t75 * qJD(5) - t80;
t225 = (t49 + t119) * t93 + (t48 + t120) * t95;
t14 = qJ(6) * t83 + t17;
t85 = t96 * t161;
t136 = pkin(5) * t85;
t2 = -t136 + t141;
t224 = -t14 * t83 + t2;
t169 = qJD(4) * t96;
t39 = t94 * t48;
t177 = t75 * t169 - t39;
t185 = t94 * t95;
t223 = t82 * t185 + t93 * t67;
t222 = t137 * t95;
t220 = t49 - t119;
t90 = t96 ^ 2;
t81 = t95 * t90 * t161;
t217 = -t109 * t83 + t81;
t16 = -t44 * t93 + t54 * t95;
t162 = qJD(6) - t16;
t13 = -pkin(5) * t83 + t162;
t128 = t13 * t93 + t14 * t95;
t46 = t51 * qJD(4);
t5 = pkin(5) * t49 + qJ(6) * t48 - qJD(6) * t75 + t46;
t216 = qJD(4) * t128 - t5;
t126 = t16 * t93 - t17 * t95;
t213 = qJD(4) * t126 + t46;
t170 = qJD(4) * t94;
t146 = t95 * t166;
t186 = t93 * t96;
t178 = t75 * t146 - t48 * t186;
t189 = t75 * t93;
t191 = t73 * t95;
t212 = -(t189 + t191) * t170 + t178;
t196 = t46 * t96;
t210 = (t218 * t94 - t51 * t96) * qJD(4) - t45 * t94 + t196;
t209 = qJD(5) * t223 - t71 * t95;
t208 = t75 ^ 2;
t207 = 0.2e1 * qJD(3);
t206 = pkin(8) * t75;
t205 = t5 * t93;
t204 = t5 * t95;
t43 = -qJD(4) * pkin(4) - t218;
t15 = pkin(5) * t73 - qJ(6) * t75 + t43;
t202 = t15 * t75;
t200 = t43 * t93;
t199 = t43 * t95;
t198 = t46 * t93;
t197 = t46 * t95;
t194 = t49 * t95;
t190 = t75 * t73;
t38 = t94 * t49;
t184 = t95 * t67;
t97 = qJD(4) ^ 2;
t182 = t97 * t94;
t181 = t97 * t96;
t130 = pkin(5) * t93 - qJ(6) * t95;
t180 = t163 + (-t130 * qJD(1) + t70) * t94 - t130 * qJD(5) + t93 * qJD(6);
t63 = t73 * t169;
t179 = t38 + t63;
t77 = t133 * qJD(1);
t24 = t218 * t95 + t93 * t77;
t176 = t94 ^ 2 - t90;
t98 = qJD(1) ^ 2;
t175 = -t97 - t98;
t79 = qJD(1) * t84;
t174 = qJD(1) * t90;
t165 = t79 * qJD(1);
t160 = pkin(8) * t188;
t159 = pkin(8) * t83 * t95;
t157 = t83 * t185;
t155 = t96 * t98 * t94;
t154 = t96 * t82 * t164 + t67 * t167 + t93 * t71;
t153 = pkin(8) * t169;
t151 = t93 * t170;
t150 = t93 * t169;
t62 = t73 * t170;
t149 = t82 * t168;
t145 = t83 * t172;
t144 = t73 ^ 2 - t208;
t143 = t82 * t93 - pkin(5);
t139 = qJD(5) * t73 - t48;
t138 = t83 * t146;
t135 = t94 * t85;
t134 = qJ(6) * t85;
t131 = pkin(5) * t95 + qJ(6) * t93;
t129 = t13 * t95 - t14 * t93;
t127 = t16 * t95 + t17 * t93;
t23 = -t218 * t93 + t77 * t95;
t122 = t139 * pkin(8);
t118 = t79 * t207;
t32 = t49 * t183;
t117 = t73 * t147 - t32;
t116 = t83 * t151 - t138;
t114 = t130 - t82;
t113 = -t15 * t94 + t153;
t112 = t43 * t94 - t153;
t110 = -t54 * t167 + t44 * t168 - t95 * t45 - t93 * t59;
t108 = t93 * t120 - t194;
t107 = -t90 * t142 + t116;
t1 = qJD(6) * t83 - t110 + t134;
t106 = t129 * qJD(5) + t1 * t95 + t2 * t93;
t105 = -t127 * qJD(5) - t110 * t95 + t141 * t93;
t103 = t49 * t186 + (t146 - t151) * t73;
t102 = qJD(4) * t15 + t106;
t101 = qJD(4) * t43 + t105;
t100 = t62 + (-t49 - t80) * t96 + (-t150 - t222) * t83;
t99 = -t49 * t185 - t95 * t63 + t75 * t222 + (t137 * t73 + t177) * t93;
t78 = -pkin(4) - t131;
t55 = (t83 + t173) * t169;
t47 = t114 * t96;
t37 = pkin(8) * t194;
t34 = pkin(5) * t75 + qJ(6) * t73;
t28 = -t82 * t187 + t184;
t27 = t143 * t94 - t184;
t26 = qJ(6) * t94 + t223;
t22 = t120 - t48;
t21 = -pkin(5) * t172 - t23;
t20 = qJ(6) * t172 + t24;
t19 = t83 * t167 + (t157 + (-t75 + t171) * t96) * qJD(1);
t18 = (t131 * qJD(5) - qJD(6) * t95) * t96 - t114 * t170;
t12 = -t82 * t150 - t209;
t11 = -t94 * t149 + t154;
t10 = t95 * t119 - t48 * t93;
t9 = -t109 * t75 - t48 * t183;
t8 = t143 * t169 + t209;
t7 = qJ(6) * t169 + (qJD(6) - t149) * t94 + t154;
t6 = t177 + t217;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t207, t118, -0.2e1 * t135, 0.2e1 * t176 * t161, -t182, 0.2e1 * t135, -t181, 0, t79 * t169 - t82 * t182 + (t84 * t169 + t94 * t207) * qJD(1), -t79 * t170 - t82 * t181 + (-t84 * t170 + t96 * t207) * qJD(1), t210, -t210 * t82 + t118, t9, t117 - t212, t6, t103, -t38 + (-t174 * t93 - t73 * t96) * qJD(4) + t116, t55, t12 * t83 + (-t141 + (t73 * t82 - t200) * qJD(4)) * t94 + (t43 * t167 + t198 - t49 * t82 + (qJD(1) * t28 + t16) * qJD(4)) * t96, -t11 * t83 + (t110 + (t75 * t82 - t199) * qJD(4)) * t94 + (-t43 * t168 + t197 + t48 * t82 + (-qJD(1) * t223 - t17) * qJD(4)) * t96, -t11 * t73 - t12 * t75 + t28 * t48 - t223 * t49 + t127 * t170 + (qJD(5) * t126 + t110 * t93 + t141 * t95) * t96, t11 * t17 + t12 * t16 - t28 * t141 - t223 * t110 + (t170 * t43 - t196) * t82, t9, t6 (-t168 * t73 + t194) * t96 + t212, t55, -t107 + t179, t103, t18 * t73 + t47 * t49 - t8 * t83 + (-t15 * t171 - t2) * t94 + (t15 * t167 + t205 + (-qJD(1) * t27 - t13) * qJD(4)) * t96, -t26 * t49 - t27 * t48 - t7 * t73 + t75 * t8 - t129 * t170 + (-qJD(5) * t128 - t1 * t93 + t2 * t95) * t96, -t18 * t75 + t47 * t48 + t7 * t83 + (t15 * t164 + t1) * t94 + (t15 * t168 - t204 + (qJD(1) * t26 + t14) * qJD(4)) * t96, t1 * t26 + t13 * t8 + t14 * t7 + t15 * t18 + t2 * t27 + t47 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t182, 0, t45 * t96 + t46 * t94 + (-t218 * t96 - t51 * t94) * qJD(4), 0, 0, 0, 0, 0, 0, t107 + t179, t177 - t217, -t32 + (t166 * t75 + t62) * t95 + (t139 * t96 - t170 * t75) * t93, t101 * t96 + t213 * t94, 0, 0, 0, 0, 0, 0, -t138 + (t83 * t94 - t174) * t171 + t179 (-t189 + t191) * t170 + t117 + t178, -t83 * t147 + t39 + t81 + (-t75 * t96 - t157) * qJD(4), t102 * t96 - t216 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t165, 0, 0, 0, 0, 0, 0, t175 * t94, t175 * t96, 0, -t210 - t165, 0, 0, 0, 0, 0, 0, t100, -t226, t99, -t127 * qJD(1) + t101 * t94 - t213 * t96, 0, 0, 0, 0, 0, 0, t100, t99, t226, t129 * qJD(1) + t102 * t94 + t216 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t176 * t98, 0, -t155, 0, 0, -t96 * t165, t94 * t165, 0, 0, t10, -t225, t19, t108, t227, -t145, -pkin(4) * t49 - t23 * t83 - t197 - t51 * t73 + (-t159 + t200) * qJD(5) + (t112 * t93 - t16 * t96) * qJD(1), pkin(4) * t48 + t24 * t83 + t198 - t51 * t75 + (t160 + t199) * qJD(5) + (t112 * t95 + t17 * t96) * qJD(1), t23 * t75 + t24 * t73 - t37 + (-t16 * t173 - t110 + (-t16 + t206) * qJD(5)) * t95 + (t122 - t111) * t93, -pkin(4) * t46 + pkin(8) * t105 - t16 * t23 - t17 * t24 - t43 * t51, t10, t19, t225, -t145, -t227, t108, t21 * t83 + t49 * t78 - t204 - t180 * t73 + (t15 * t93 - t159) * qJD(5) + (-t113 * t93 + t13 * t96) * qJD(1), t20 * t73 - t21 * t75 - t37 + (t13 * t173 + t1 + (t13 + t206) * qJD(5)) * t95 + (t122 + t224) * t93, -t20 * t83 + t48 * t78 - t205 + t180 * t75 + (-t15 * t95 - t160) * qJD(5) + (t113 * t95 - t14 * t96) * qJD(1), pkin(8) * t106 - t13 * t21 - t14 * t20 - t15 * t180 + t5 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, -t144, t22, -t190, -t220, t85, -t43 * t75 + t111, t16 * t83 + t43 * t73 + t110, 0, 0, t190, t22, t144, t85, t220, -t190, -t34 * t73 + t111 + 0.2e1 * t136 - t202, pkin(5) * t48 - qJ(6) * t49 + (t14 - t17) * t75 + (t13 - t162) * t73, 0.2e1 * t134 - t15 * t73 + t34 * t75 + (0.2e1 * qJD(6) - t16) * t83 - t110, -pkin(5) * t2 + qJ(6) * t1 - t13 * t17 + t14 * t162 - t15 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 + t190, t22, -t83 ^ 2 - t208, t202 + t224;];
tauc_reg  = t3;
