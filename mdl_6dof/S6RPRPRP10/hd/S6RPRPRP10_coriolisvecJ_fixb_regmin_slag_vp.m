% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:32
% EndTime: 2019-03-09 03:32:38
% DurationCPUTime: 1.98s
% Computational Cost: add. (2206->316), mult. (4260->403), div. (0->0), fcn. (2205->4), ass. (0->163)
t100 = sin(qJ(3));
t101 = cos(qJ(5));
t102 = cos(qJ(3));
t160 = qJD(1) * t102;
t85 = qJD(5) + t160;
t118 = t100 * (t85 + t160);
t130 = qJD(5) * t102 + qJD(1);
t153 = qJD(1) * qJD(3);
t138 = t102 * t153;
t161 = qJD(1) * t100;
t139 = t101 * t161;
t99 = sin(qJ(5));
t169 = qJD(5) * t99;
t30 = qJD(3) * t169 - qJD(5) * t139 - t99 * t138;
t158 = qJD(3) * t101;
t61 = t99 * t161 + t158;
t213 = t130 * t101 * t85 + (t102 * t61 - t99 * t118) * qJD(3) - t100 * t30;
t104 = -pkin(1) - pkin(7);
t79 = t104 * qJD(1) + qJD(2);
t67 = t100 * t79;
t95 = qJD(3) * qJ(4);
t53 = -t67 - t95;
t179 = t102 * t53;
t157 = qJD(3) * t102;
t93 = qJD(3) * qJD(4);
t49 = -t79 * t157 - t93;
t183 = qJD(3) * pkin(3);
t132 = -qJD(4) + t183;
t68 = t102 * t79;
t51 = -t132 - t68;
t211 = ((-t51 + t68) * t100 + t179) * qJD(3) + t49 * t100;
t171 = qJ(4) * t102;
t120 = pkin(8) * t100 - t171;
t209 = t100 * pkin(3) + qJ(2);
t55 = t120 + t209;
t189 = pkin(4) - t104;
t72 = t189 * t102;
t186 = t101 * t55 + t99 * t72;
t59 = qJD(3) * t99 - t139;
t193 = t59 * t85;
t208 = t30 - t193;
t191 = t61 * t85;
t170 = qJD(5) * t61;
t31 = -t101 * t138 + t170;
t207 = -t31 + t191;
t206 = qJD(4) - t68;
t103 = -pkin(3) - pkin(8);
t156 = qJD(3) * t103;
t200 = t85 ^ 2;
t137 = t100 * t153;
t75 = t99 * t137;
t204 = qJD(3) * t61 + t200 * t101 - t75;
t154 = qJD(5) * t101;
t117 = (qJD(3) * pkin(8) - qJD(4)) * t102;
t94 = qJD(1) * qJD(2);
t151 = pkin(3) * t138 + qJ(4) * t137 + t94;
t20 = qJD(1) * t117 + t151;
t163 = pkin(4) * t160 + t206;
t27 = t163 + t156;
t159 = qJD(3) * t100;
t63 = t79 * t159;
t40 = -pkin(4) * t137 + t63;
t152 = -t101 * t20 - t27 * t154 - t99 * t40;
t185 = pkin(3) * t161 + qJD(1) * qJ(2);
t36 = t120 * qJD(1) + t185;
t112 = t36 * t169 + t152;
t129 = qJ(6) * t137;
t1 = qJD(6) * t85 - t112 - t129;
t9 = t101 * t27 - t36 * t99;
t172 = qJD(6) - t9;
t6 = -pkin(5) * t85 + t172;
t10 = t101 * t36 + t27 * t99;
t7 = qJ(6) * t85 + t10;
t126 = t101 * t7 + t6 * t99;
t135 = -t101 * t40 + t36 * t154 + t27 * t169 + t99 * t20;
t83 = pkin(5) * t137;
t2 = t83 + t135;
t202 = qJD(5) * t126 + t1 * t99 - t2 * t101;
t201 = t61 ^ 2;
t199 = 0.2e1 * t94;
t33 = -pkin(4) * t138 - t49;
t3 = pkin(5) * t31 + qJ(6) * t30 - qJD(6) * t61 + t33;
t198 = t3 * t99;
t197 = t7 * t85;
t196 = t10 * t85;
t195 = t3 * t101;
t194 = t33 * t99;
t192 = t61 * t59;
t190 = t61 * t99;
t123 = pkin(5) * t101 + qJ(6) * t99;
t115 = -pkin(4) - t123;
t110 = t115 * t102;
t188 = qJD(1) * t110 - t123 * qJD(5) + t101 * qJD(6) - t206;
t65 = pkin(3) * t160 + qJ(4) * t161;
t45 = pkin(8) * t160 + t65;
t46 = -pkin(4) * t161 + t67;
t187 = t101 * t45 + t99 * t46;
t97 = t100 ^ 2;
t98 = t102 ^ 2;
t184 = t97 - t98;
t181 = t100 * t59;
t177 = t102 * t99;
t176 = t103 * t85;
t175 = t30 * t101;
t174 = t33 * t101;
t168 = t101 * t102;
t105 = qJD(3) ^ 2;
t167 = t105 * t100;
t166 = t105 * t102;
t106 = qJD(1) ^ 2;
t165 = t106 * qJ(2);
t38 = t46 + t95;
t11 = pkin(5) * t59 - qJ(6) * t61 + t38;
t164 = t11 * qJD(3);
t162 = t105 + t106;
t155 = qJD(4) * t102;
t150 = 0.2e1 * qJD(1);
t149 = t85 * t177;
t148 = t99 * t176;
t147 = t85 * t169;
t146 = t101 * t176;
t145 = pkin(3) * t157 + qJ(4) * t159 + qJD(2);
t144 = t99 * t156;
t143 = t85 * t154;
t142 = t102 * t106 * t100;
t141 = t100 * t158;
t140 = t101 * t156;
t52 = -qJ(4) * t160 + t185;
t70 = -t171 + t209;
t134 = qJD(1) * t70 + t52;
t128 = t102 * t137;
t125 = t101 * t6 - t7 * t99;
t124 = pkin(5) * t99 - qJ(6) * t101;
t122 = t101 * t46 - t45 * t99;
t121 = t101 * t72 - t55 * t99;
t119 = -qJD(1) * t97 + t102 * t85;
t114 = t11 * t61 + t135;
t29 = -qJD(1) * t155 + t151;
t43 = t145 - t155;
t113 = -qJD(1) * t43 + t104 * t105 - t29;
t32 = t117 + t145;
t57 = t189 * t159;
t111 = t101 * t32 + t72 * t154 - t55 * t169 - t99 * t57;
t108 = t100 * t31 + t59 * t157 + t101 * t128 + (t130 * t99 + t141) * t85;
t107 = -t147 - qJD(3) * t59 + (-t141 - t149) * qJD(1);
t88 = t100 * t104;
t82 = t104 * t157;
t74 = t162 * t102;
t73 = t162 * t100;
t71 = -pkin(4) * t100 + t88;
t69 = qJ(4) + t124;
t58 = -pkin(4) * t157 + t82;
t41 = t52 * t160;
t28 = t115 * t100 + t88;
t21 = pkin(5) * t61 + qJ(6) * t59;
t16 = -pkin(5) * t102 - t121;
t15 = qJ(6) * t102 + t186;
t13 = pkin(5) * t161 - t122;
t12 = -qJ(6) * t161 + t187;
t8 = t82 + (t124 * qJD(5) - qJD(6) * t99) * t100 + qJD(3) * t110;
t5 = pkin(5) * t159 + t186 * qJD(5) + t101 * t57 + t99 * t32;
t4 = -qJ(6) * t159 + qJD(6) * t102 + t111;
t14 = [0, 0, 0, 0, t199, qJ(2) * t199, -0.2e1 * t128, 0.2e1 * t184 * t153, -t167, -t166, 0, -t104 * t167 + (qJ(2) * t157 + qJD(2) * t100) * t150, -t104 * t166 + (-qJ(2) * t159 + qJD(2) * t102) * t150, t211, t113 * t100 - t134 * t157, t113 * t102 + t134 * t159, -t104 * t211 + t29 * t70 + t52 * t43, t157 * t190 + (t61 * t154 - t30 * t99) * t100 (t101 * t61 - t59 * t99) * t157 + (-t175 - t99 * t31 + (-t101 * t59 - t190) * qJD(5)) * t100, t100 * t143 - t30 * t102 + (-t100 * t61 + t119 * t99) * qJD(3), -t100 * t147 - t31 * t102 + (t119 * t101 + t181) * qJD(3), -qJD(3) * t118, -t135 * t102 + t58 * t59 + t71 * t31 + (-qJD(5) * t72 - t32) * t85 * t99 + ((-qJD(5) * t55 - t57) * t85 - t38 * t157) * t101 + (t38 * t169 - t174 + (-t121 * qJD(1) - t9) * qJD(3)) * t100, -t111 * t85 + t58 * t61 - t71 * t30 + ((qJD(3) * t38 + qJD(5) * t36) * t99 + t152) * t102 + (t38 * t154 + t194 + (t186 * qJD(1) + t10) * qJD(3)) * t100, t28 * t31 - t5 * t85 + t8 * t59 + (-t11 * t158 - t2) * t102 + (t11 * t169 - t195 + (qJD(1) * t16 + t6) * qJD(3)) * t100, -t15 * t31 - t16 * t30 - t4 * t59 + t5 * t61 + t126 * t157 + (qJD(5) * t125 + t1 * t101 + t2 * t99) * t100, t28 * t30 + t4 * t85 - t8 * t61 + (-t99 * t164 + t1) * t102 + (-t11 * t154 - t198 + (-qJD(1) * t15 - t7) * qJD(3)) * t100, t1 * t15 + t11 * t8 + t16 * t2 + t28 * t3 + t4 * t7 + t5 * t6; 0, 0, 0, 0, -t106, -t165, 0, 0, 0, 0, 0, -t73, -t74, 0, t73, t74, -t52 * qJD(1) - t211, 0, 0, 0, 0, 0, t108, t213, t108 (-t59 * t159 - qJD(1) * t61 + (t31 - t170) * t102) * t99 + (-t61 * t159 + qJD(1) * t59 + (qJD(5) * t59 - t30) * t102) * t101, -t213, -t126 * qJD(1) + (-qJD(3) * t125 + t3) * t100 + (t164 - t202) * t102; 0, 0, 0, 0, 0, 0, t142, -t184 * t106, 0, 0, 0, -t102 * t165, t100 * t165 ((-t53 - t95) * t102 + (t132 + t51) * t100) * qJD(1), t65 * t161 + t41, 0.2e1 * t93 + (-t100 * t52 + t102 * t65) * qJD(1), -t49 * qJ(4) - t53 * qJD(4) - t52 * t65 + (t179 + (-t51 - t183) * t100) * t79, -t191 * t99 - t175 (t30 + t193) * t99 + (-t31 - t191) * t101, -t147 + (-t149 + (t61 - t158) * t100) * qJD(1), -t143 + t75 + (-t85 * t168 - t181) * qJD(1), t85 * t161, qJ(4) * t31 + t194 - t122 * t85 + t163 * t59 + (t101 * t38 - t148) * qJD(5) + (t38 * t168 + (t9 - t140) * t100) * qJD(1), -qJ(4) * t30 + t174 + t187 * t85 + t163 * t61 + (-t38 * t99 - t146) * qJD(5) + (-t38 * t177 + (-t10 + t144) * t100) * qJD(1), t13 * t85 + t198 + t69 * t31 - t188 * t59 + (t101 * t11 - t148) * qJD(5) + (t11 * t168 + (-t6 - t140) * t100) * qJD(1), t12 * t59 - t13 * t61 + (-t6 * t160 - t103 * t31 - t1 + (t103 * t61 - t6) * qJD(5)) * t99 + (-t7 * t160 + t103 * t30 + t2 + (-t103 * t59 - t7) * qJD(5)) * t101, -t195 - t12 * t85 + t69 * t30 + t188 * t61 + (t11 * t99 + t146) * qJD(5) + (t11 * t177 + (t7 - t144) * t100) * qJD(1), t202 * t103 - t188 * t11 - t7 * t12 - t6 * t13 + t3 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t106 * t98 - t105, qJD(3) * t53 + t41 + t63, 0, 0, 0, 0, 0, t107, -t204, t107, t208 * t101 + t207 * t99, t204, -t164 + (t6 * t85 + t1) * t99 + (-t2 + t197) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, -t59 ^ 2 + t201, -t208, t207, -t137, -t38 * t61 - t135 + t196, t38 * t59 + t85 * t9 + t112, -t21 * t59 - t114 + t196 - 0.2e1 * t83, pkin(5) * t30 - t31 * qJ(6) + (-t10 + t7) * t61 + (t6 - t172) * t59, -0.2e1 * t129 - t11 * t59 + t21 * t61 + (0.2e1 * qJD(6) - t9) * t85 - t112, -t2 * pkin(5) + t1 * qJ(6) - t6 * t10 - t11 * t21 + t172 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 + t192, -t208, -t200 - t201, t114 + t83 - t197;];
tauc_reg  = t14;
