% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP1
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:07
% DurationCPUTime: 1.96s
% Computational Cost: add. (2843->226), mult. (7495->314), div. (0->0), fcn. (5391->6), ass. (0->154)
t118 = qJD(2) + qJD(3);
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t126 = cos(qJ(2));
t177 = qJD(1) * t126;
t123 = sin(qJ(2));
t178 = qJD(1) * t123;
t212 = -t122 * t178 + t125 * t177;
t64 = t212 * t118;
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t87 = -t122 * t177 - t125 * t178;
t147 = -t121 * t212 + t124 * t87;
t206 = t147 ^ 2;
t58 = t121 * t87 + t124 * t212;
t54 = t58 ^ 2;
t10 = -t54 + t206;
t205 = pkin(6) + pkin(7);
t105 = t205 * t126;
t101 = qJD(1) * t105;
t88 = t122 * t101;
t188 = qJD(2) * pkin(2);
t104 = t205 * t123;
t99 = qJD(1) * t104;
t94 = -t99 + t188;
t153 = t125 * t94 - t88;
t81 = t87 * pkin(8);
t41 = t153 + t81;
t31 = pkin(3) * t118 + t41;
t92 = t125 * t101;
t146 = -t122 * t94 - t92;
t201 = pkin(8) * t212;
t42 = -t146 + t201;
t37 = t124 * t42;
t149 = -t121 * t31 - t37;
t98 = t122 * t126 + t123 * t125;
t137 = t98 * qJD(2);
t131 = -qJD(3) * t98 - t137;
t130 = t131 * qJD(1);
t176 = qJD(3) * t122;
t171 = qJD(2) * t205;
t151 = qJD(1) * t171;
t96 = t126 * t151;
t158 = -t101 * t176 - t122 * t96;
t95 = t123 * t151;
t208 = t125 * (qJD(3) * t94 - t95);
t20 = pkin(8) * t130 + t158 + t208;
t159 = t122 * t95 - t125 * t96;
t132 = qJD(3) * t146 + t159;
t21 = -pkin(8) * t64 + t132;
t165 = -t121 * t20 + t124 * t21;
t135 = qJD(4) * t149 + t165;
t114 = -t126 * pkin(2) - pkin(1);
t103 = t114 * qJD(1);
t71 = -pkin(3) * t212 + t103;
t197 = t71 * t147;
t214 = t135 + t197;
t175 = qJD(4) * t121;
t164 = t121 * t21 - t42 * t175;
t136 = -(qJD(4) * t31 + t20) * t124 - t164;
t196 = t71 * t58;
t213 = t136 - t196;
t189 = qJ(5) * t58;
t198 = t147 * t58;
t48 = t147 * qJ(5);
t117 = qJD(4) + t118;
t174 = qJD(4) * t124;
t141 = t121 * t130 + t124 * t64 + t174 * t212 + t87 * t175;
t5 = -t117 * t58 + t141;
t134 = qJD(4) * t147 - t121 * t64 + t124 * t130;
t6 = -t117 * t147 + t134;
t173 = qJD(2) * qJD(1);
t211 = -0.2e1 * t173;
t113 = pkin(2) * t125 + pkin(3);
t183 = t122 * t124;
t157 = t122 * t99 - t92;
t44 = t157 - t201;
t190 = -t125 * t99 - t88;
t45 = t81 + t190;
t210 = t121 * t45 - t124 * t44 - t113 * t175 + (-t122 * t174 + (-t121 * t125 - t183) * qJD(3)) * pkin(2);
t185 = t121 * t122;
t209 = -t113 * t174 - (-t122 * t175 + (t124 * t125 - t185) * qJD(3)) * pkin(2) + t121 * t44 + t124 * t45;
t207 = -pkin(8) * t98 - t122 * t105;
t35 = t121 * t42;
t163 = t124 * t31 - t35;
t8 = t163 + t48;
t7 = pkin(4) * t117 + t8;
t204 = t7 - t8;
t203 = pkin(3) * t87;
t202 = pkin(4) * t147;
t199 = t123 * pkin(2);
t195 = t87 * t212;
t194 = t124 * t41 - t35;
t192 = -t48 - t209;
t191 = t189 + t210;
t187 = t103 * t87;
t186 = t104 * t125;
t128 = qJD(1) ^ 2;
t182 = t126 * t128;
t127 = qJD(2) ^ 2;
t181 = t127 * t123;
t180 = t127 * t126;
t179 = t123 ^ 2 - t126 ^ 2;
t100 = t123 * t171;
t102 = t126 * t171;
t172 = -qJD(3) * t186 - t125 * t100 - t122 * t102;
t116 = t123 * t188;
t115 = pkin(2) * t178;
t168 = -pkin(2) * t118 - t94;
t72 = t115 - t203;
t167 = -pkin(3) * t117 - t31;
t166 = t123 * t173;
t162 = -t121 * t41 - t37;
t154 = t100 * t122 - t125 * t102;
t152 = pkin(1) * t211;
t9 = -t149 + t189;
t150 = -t147 * t9 + t58 * t7;
t49 = -t186 + t207;
t144 = t122 * t123 - t125 * t126;
t145 = t122 * t104 - t125 * t105;
t50 = -pkin(8) * t144 - t145;
t148 = -t121 * t49 - t124 * t50;
t143 = -t103 * t212 - t158;
t27 = -pkin(8) * t137 + t207 * qJD(3) + t172;
t70 = t118 * t144;
t28 = pkin(8) * t70 + qJD(3) * t145 + t154;
t142 = t121 * t28 + t124 * t27 + t49 * t174 - t50 * t175;
t140 = t103 * t98;
t139 = t114 * t98;
t138 = t124 * t144;
t69 = -t121 * t144 + t124 * t98;
t74 = pkin(3) * t144 + t114;
t133 = qJD(4) * t148 - t121 * t27 + t124 * t28;
t61 = -pkin(3) * t131 + t116;
t46 = pkin(2) * t166 - pkin(3) * t130;
t129 = -pkin(4) * t134 + t46;
t112 = pkin(3) * t124 + pkin(4);
t82 = pkin(2) * t183 + t113 * t121;
t79 = -pkin(2) * t185 + t113 * t124 + pkin(4);
t68 = t121 * t98 + t138;
t43 = -t212 ^ 2 + t87 ^ 2;
t34 = -t118 * t87 + t130;
t29 = -pkin(4) * t58 + qJD(5) + t71;
t23 = qJD(4) * t69 - t121 * t70 - t124 * t131;
t22 = qJD(4) * t138 - t121 * t131 + t124 * t70 + t98 * t175;
t16 = -qJ(5) * t68 - t148;
t15 = -qJ(5) * t69 - t121 * t50 + t124 * t49;
t12 = t48 + t194;
t11 = t162 - t189;
t4 = qJ(5) * t22 - qJD(5) * t69 + t133;
t3 = -qJ(5) * t23 - qJD(5) * t68 + t142;
t2 = -qJ(5) * t141 + qJD(5) * t147 + t135;
t1 = qJ(5) * t134 + qJD(5) * t58 - t136;
t13 = [0, 0, 0, 0.2e1 * t126 * t166, t179 * t211, t180, -t181, 0, -pkin(6) * t180 + t123 * t152, pkin(6) * t181 + t126 * t152, t64 * t98 + t70 * t87, t98 * t130 - t87 * t131 - t144 * t64 - t212 * t70, -t70 * t118, t131 * t118, 0, t154 * t118 + (t118 * t145 + t140) * qJD(3) + (-t199 * t212 + t140) * qJD(2) + (qJD(3) * t139 + (t144 * t199 + t139) * qJD(2)) * qJD(1), t114 * t64 - t103 * t70 - (-t105 * t176 + t172) * t118 + (qJD(1) * t98 - t87) * t116, t141 * t69 + t147 * t22, t134 * t69 - t141 * t68 + t147 * t23 - t22 * t58, -t22 * t117, -t23 * t117, 0, t117 * t133 - t134 * t74 + t71 * t23 + t46 * t68 - t58 * t61, -t117 * t142 + t141 * t74 - t147 * t61 - t71 * t22 + t46 * t69, -t1 * t68 + t134 * t16 - t141 * t15 + t147 * t4 - t2 * t69 + t22 * t7 - t23 * t9 + t3 * t58, t1 * t16 + t9 * t3 + t2 * t15 + t7 * t4 + t129 * (t68 * pkin(4) + t74) + t29 * (t23 * pkin(4) + t61); 0, 0, 0, -t123 * t182, t179 * t128, 0, 0, 0, t128 * pkin(1) * t123, pkin(1) * t182, t195, t43, 0, t34, 0, t212 * t115 + t187 - t157 * t118 + (t168 * t122 - t92) * qJD(3) + t159, t87 * t115 + t190 * t118 + (t168 * qJD(3) + t95) * t125 + t143, t198, t10, t5, t6, 0, t210 * t117 + t58 * t72 + t214, t209 * t117 + t147 * t72 + t213, t134 * t82 - t141 * t79 + t147 * t191 + t192 * t58 + t150, t1 * t82 + t2 * t79 - t29 * (t72 - t202) + t192 * t9 + t191 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t43, 0, t34, 0, -t118 * t146 + t132 + t187, t118 * t153 + t143 - t208, t198, t10, t5, t6, 0, -t58 * t203 + t197 - t162 * t117 + (t121 * t167 - t37) * qJD(4) + t165, -t147 * t203 - t196 + t194 * t117 + (qJD(4) * t167 - t20) * t124 - t164, -t11 * t147 - t112 * t141 - t12 * t58 + (t121 * t134 + (-t121 * t147 + t124 * t58) * qJD(4)) * pkin(3) + t150, t29 * t202 - t7 * t11 + t2 * t112 - t9 * t12 + (t1 * t121 + t29 * t87 + (-t121 * t7 + t124 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t10, t5, t6, 0, -t117 * t149 + t214, t117 * t163 + t213, -pkin(4) * t141 + t204 * t58, t204 * t9 + (t147 * t29 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 - t206, -t147 * t7 - t58 * t9 + t129;];
tauc_reg = t13;
