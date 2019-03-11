% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:19
% EndTime: 2019-03-09 02:06:24
% DurationCPUTime: 1.74s
% Computational Cost: add. (2532->294), mult. (4889->398), div. (0->0), fcn. (2800->6), ass. (0->153)
t98 = cos(qJ(4));
t165 = qJD(1) * t98;
t80 = qJD(5) + t165;
t160 = qJD(4) * t98;
t95 = sin(qJ(5));
t156 = t95 * qJD(4);
t96 = sin(qJ(4));
t166 = qJD(1) * t96;
t97 = cos(qJ(5));
t66 = t166 * t97 - t156;
t139 = t66 * t160;
t162 = qJD(4) * t96;
t92 = sin(pkin(9));
t141 = t92 * t162;
t177 = t97 * t98;
t180 = t95 * t98;
t93 = cos(pkin(9));
t62 = t180 * t92 + t93 * t97;
t174 = -t62 * qJD(5) - t97 * t141 - (t177 * t93 + t92 * t95) * qJD(1);
t152 = qJD(1) * qJD(4);
t133 = t98 * t152;
t161 = qJD(4) * t97;
t65 = t166 * t95 + t161;
t159 = qJD(5) * t65;
t43 = -t133 * t97 + t159;
t63 = t177 * t92 - t93 * t95;
t206 = (t43 * t92 + (qJD(4) * t63 + t66 * t93) * qJD(1)) * t96 - t92 * t139 - t174 * t80;
t154 = qJ(2) * qJD(1);
t99 = -pkin(1) - pkin(2);
t78 = qJD(1) * t99 + qJD(2);
t57 = t154 * t93 + t78 * t92;
t47 = -qJD(1) * pkin(7) + t57;
t205 = qJD(3) * t98 - t47 * t96;
t158 = qJD(5) * t95;
t138 = t96 * t158;
t108 = t160 * t97 - t138;
t179 = t96 * t97;
t157 = qJD(5) * t97;
t137 = t96 * t157;
t107 = t156 * t98 + t137;
t44 = qJD(1) * t107 - qJD(5) * t156;
t204 = -t108 * t65 - t179 * t44;
t90 = t96 ^ 2;
t134 = t90 * t152;
t114 = t134 * t97 + t138 * t80;
t148 = t80 * t177;
t184 = t66 * t96;
t189 = t43 * t98;
t201 = qJD(4) * (t148 + t184) - t114 + t189;
t200 = t66 ^ 2;
t153 = qJD(1) * qJD(2);
t132 = t93 * t153;
t89 = t96 * qJD(3);
t21 = qJD(4) * t89 + t132 * t96 + t160 * t47;
t3 = -pkin(5) * t44 - qJ(6) * t43 + qJD(6) * t66 + t21;
t199 = t3 * t95;
t198 = t3 * t97;
t32 = t47 * t98 + t89;
t27 = qJD(4) * pkin(8) + t32;
t124 = pkin(4) * t98 + pkin(8) * t96;
t56 = -t154 * t92 + t93 * t78;
t46 = qJD(1) * pkin(3) - t56;
t28 = qJD(1) * t124 + t46;
t9 = t27 * t97 + t28 * t95;
t7 = qJ(6) * t80 + t9;
t197 = t7 * t80;
t196 = t80 * t9;
t195 = t21 * t95;
t194 = t21 * t97;
t26 = -qJD(4) * pkin(4) - t205;
t193 = t26 * t95;
t192 = t26 * t97;
t190 = t43 * t95;
t188 = t44 * t98;
t187 = t65 * t80;
t186 = t66 * t65;
t185 = t66 * t80;
t183 = t80 * t95;
t182 = t80 * t97;
t181 = t92 * t96;
t130 = -qJ(2) * t92 + t93 * t99;
t67 = pkin(3) - t130;
t49 = t124 + t67;
t178 = t97 * t49;
t119 = pkin(5) * t95 - qJ(6) * t97;
t176 = t95 * qJD(6) - t119 * t80 + t32;
t123 = -pkin(4) * t96 + pkin(8) * t98;
t69 = t123 * qJD(1);
t175 = t205 * t97 + t69 * t95;
t173 = t63 * qJD(5) - t95 * t141 - (t180 * t93 - t92 * t97) * qJD(1);
t171 = qJ(2) * t93 + t92 * t99;
t68 = -pkin(7) + t171;
t172 = t177 * t68 + t49 * t95;
t170 = -t98 ^ 2 + t90;
t100 = qJD(4) ^ 2;
t169 = t100 * t96;
t168 = t100 * t98;
t8 = -t27 * t95 + t28 * t97;
t167 = qJD(6) - t8;
t164 = qJD(2) * t93;
t101 = qJD(1) ^ 2;
t155 = t100 + t101;
t151 = pkin(8) * t183;
t150 = pkin(8) * t182;
t149 = t80 * t180;
t20 = qJD(4) * t205 + t132 * t98;
t61 = t92 * qJD(2) + qJD(4) * t123;
t48 = t61 * qJD(1);
t146 = -t157 * t28 - t20 * t97 - t48 * t95;
t145 = -t137 * t66 - t139 * t95 + t190 * t96;
t142 = t98 * t164;
t144 = t142 * t97 + t157 * t49 + t61 * t95;
t143 = pkin(8) * t161;
t140 = t92 * t160;
t136 = t80 * t157;
t135 = 0.2e1 * t153;
t131 = t96 * t152;
t129 = -t157 * t27 - t158 * t28 - t20 * t95 + t48 * t97;
t127 = t92 * t135;
t126 = 0.2e1 * t131;
t125 = qJ(6) * t131;
t6 = -pkin(5) * t80 + t167;
t122 = t6 * t97 - t7 * t95;
t121 = t6 * t95 + t7 * t97;
t120 = -pkin(5) * t97 - qJ(6) * t95;
t118 = -t205 * t95 + t69 * t97;
t117 = t56 * t92 - t57 * t93;
t112 = -t119 + t68;
t11 = -pkin(5) * t65 + qJ(6) * t66 + t26;
t111 = -t11 * t66 - t129;
t110 = t27 * t158 + t146;
t109 = -t100 * t68 + t127;
t105 = qJD(4) * (-qJD(1) * t67 - t164 - t46);
t1 = qJD(6) * t80 - t110 - t125;
t79 = pkin(5) * t131;
t2 = t79 - t129;
t104 = qJD(5) * t122 + t1 * t97 + t2 * t95;
t75 = t95 * t134;
t103 = -t107 * t80 - t162 * t65 + t188 + t75;
t102 = -t44 * t181 - t65 * t140 - t173 * t80 + (qJD(4) * t62 + t65 * t93) * t166;
t73 = t95 * pkin(8) * t131;
t71 = -pkin(4) + t120;
t34 = -pkin(5) * t66 - qJ(6) * t65;
t29 = t112 * t96;
t19 = t43 - t187;
t15 = -t178 + (t68 * t95 - pkin(5)) * t98;
t14 = qJ(6) * t98 + t172;
t13 = pkin(5) * t166 - t118;
t12 = -qJ(6) * t166 + t175;
t10 = t112 * t160 + (qJD(5) * t120 + qJD(6) * t97 + t164) * t96;
t5 = pkin(5) * t162 + (qJD(5) * t68 * t98 - t61) * t97 + (qJD(5) * t49 - t162 * t68 + t142) * t95;
t4 = (-t158 * t68 + qJD(6)) * t98 + (-t68 * t97 - qJ(6)) * t162 + t144;
t16 = [0, 0, 0, 0, t135, qJ(2) * t135, t127, 0.2e1 * t132 ((-t130 * t92 + t171 * t93) * qJD(1) - t117) * qJD(2), t98 * t126, -0.2e1 * t170 * t152, -t168, t169, 0, t105 * t96 + t109 * t98, t105 * t98 - t109 * t96, t108 * t66 - t179 * t43, t145 + t204, t189 + (-t148 + t184) * qJD(4) + t114, t96 * t136 + t188 - t75 + (-t65 * t96 + t149) * qJD(4) (-t80 - t165) * t162 (-t158 * t49 + t97 * t61) * t80 + ((-t157 * t68 - t164 * t95) * t80 + (-t65 * t68 - t193) * qJD(4) + t129) * t98 + (-t65 * t164 - t26 * t157 - t195 - t68 * t44 + (t68 * t183 - (-t180 * t68 + t178) * qJD(1) - t8) * qJD(4)) * t96, -t144 * t80 + ((t68 * t80 + t27) * t158 + (-t66 * t68 - t192) * qJD(4) + t146) * t98 + (-t66 * t164 + t26 * t158 - t194 + t68 * t43 + (qJD(1) * t172 + t182 * t68 + t9) * qJD(4)) * t96, -t10 * t65 - t29 * t44 - t5 * t80 + (-t11 * t156 - t2) * t98 + (-t11 * t157 - t199 + (qJD(1) * t15 + t6) * qJD(4)) * t96, t14 * t44 + t15 * t43 + t4 * t65 - t5 * t66 - t122 * t160 + (qJD(5) * t121 + t1 * t95 - t2 * t97) * t96, t10 * t66 - t29 * t43 + t4 * t80 + (t11 * t161 + t1) * t98 + (-t11 * t158 + t198 + (-qJD(1) * t14 - t7) * qJD(4)) * t96, t1 * t14 + t10 * t11 + t15 * t2 + t29 * t3 + t4 * t7 + t5 * t6; 0, 0, 0, 0, -t101, -t101 * qJ(2), -t92 * t101, -t93 * t101, t117 * qJD(1), 0, 0, 0, 0, 0, -t155 * t92 * t98 + t126 * t93, 0.2e1 * t133 * t93 + t155 * t181, 0, 0, 0, 0, 0, t102, t206, t102, -t173 * t66 + t174 * t65 + t62 * t43 + t63 * t44, -t206, t3 * t181 + t1 * t63 + t2 * t62 + t174 * t7 + t173 * t6 + (-t166 * t93 + t140) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, -t168, 0, 0, 0, 0, 0, t103, -t201, t103, t145 - t204, t201 (qJD(4) * t121 - t3) * t98 + (qJD(4) * t11 + t104) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96 * t101 * t98, t170 * t101, 0, 0, 0, qJD(4) * t32 + t166 * t46 - t21 (t46 - t164) * t165, -t182 * t66 + t190 (t43 + t187) * t97 + (t44 + t185) * t95, t136 + (t148 + (-t66 - t156) * t96) * qJD(1), -t80 * t158 + (-t149 + (t65 - t161) * t96) * qJD(1), t80 * t166, t73 + pkin(4) * t44 - t194 - t118 * t80 + t32 * t65 + (-t150 + t193) * qJD(5) + (t180 * t26 + t8 * t96) * qJD(1), -pkin(4) * t43 + t195 + t175 * t80 + t32 * t66 + (t151 + t192) * qJD(5) + (t26 * t177 + (-t9 + t143) * t96) * qJD(1), t13 * t80 - t198 - t44 * t71 + t73 + t176 * t65 + (t11 * t95 - t150) * qJD(5) + (t11 * t180 - t6 * t96) * qJD(1), -t12 * t65 + t13 * t66 + (t1 + t80 * t6 + (-qJD(5) * t66 + t44) * pkin(8)) * t97 + (t2 - t197 + (t43 - t159) * pkin(8)) * t95, -t12 * t80 - t199 - t43 * t71 - t176 * t66 + (-t11 * t97 - t151) * qJD(5) + (-t11 * t177 + (t7 - t143) * t96) * qJD(1), pkin(8) * t104 - t11 * t176 - t7 * t12 - t6 * t13 + t3 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, -t65 ^ 2 + t200, t19, t44 - t185, -t131, t26 * t66 + t129 + t196, -t26 * t65 + t8 * t80 + t110, t34 * t65 - t111 + t196 - 0.2e1 * t79, -pkin(5) * t43 + t44 * qJ(6) + (-t7 + t9) * t66 + (-t6 + t167) * t65, -0.2e1 * t125 + t11 * t65 - t34 * t66 + (0.2e1 * qJD(6) - t8) * t80 - t110, -t2 * pkin(5) + t1 * qJ(6) - t11 * t34 + t167 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131 + t186, t19, -t80 ^ 2 - t200, t111 + t79 - t197;];
tauc_reg  = t16;
