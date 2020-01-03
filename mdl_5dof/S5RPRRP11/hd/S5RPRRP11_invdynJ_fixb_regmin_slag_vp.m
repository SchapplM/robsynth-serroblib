% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:35
% EndTime: 2019-12-31 18:54:42
% DurationCPUTime: 2.64s
% Computational Cost: add. (3306->336), mult. (7762->415), div. (0->0), fcn. (5757->10), ass. (0->159)
t111 = sin(pkin(8));
t112 = cos(pkin(8));
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t223 = -t111 * t115 + t118 * t112;
t226 = t223 * qJD(1);
t110 = pkin(8) + qJ(3);
t103 = cos(t110);
t200 = g(3) * t103;
t102 = sin(t110);
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t153 = g(1) * t119 + g(2) * t116;
t216 = t153 * t102;
t127 = -t200 + t216;
t194 = pkin(6) + qJ(2);
t91 = t194 * t111;
t84 = qJD(1) * t91;
t92 = t194 * t112;
t85 = qJD(1) * t92;
t54 = -t115 * t84 + t118 * t85;
t224 = t54 * qJD(3);
t169 = qJD(1) * qJD(2);
t211 = t194 * qJDD(1) + t169;
t65 = t211 * t111;
t66 = t211 * t112;
t143 = -t115 * t66 - t118 * t65 - t224;
t22 = -qJDD(3) * pkin(3) - t143;
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t80 = t223 * qJD(3);
t131 = qJD(1) * t80;
t83 = t111 * t118 + t112 * t115;
t132 = t83 * qJDD(1);
t122 = t132 + t131;
t170 = t117 * qJD(3);
t172 = qJD(4) * t114;
t215 = t83 * qJD(1);
t26 = -qJD(4) * t170 - t114 * qJDD(3) - t117 * t122 + t172 * t215;
t171 = qJD(4) * t117;
t173 = qJD(3) * t114;
t27 = -t117 * qJDD(3) + t114 * t132 + t171 * t215 + (qJD(4) + t226) * t173;
t61 = t117 * t215 + t173;
t3 = pkin(4) * t27 + qJ(5) * t26 - qJD(5) * t61 + t22;
t225 = -t3 + t127;
t69 = qJD(4) - t226;
t182 = qJDD(1) * pkin(1);
t214 = g(1) * t116 - g(2) * t119;
t138 = -qJDD(2) + t182 + t214;
t222 = qJD(3) * t215;
t146 = -t115 * t65 + t118 * t66;
t219 = -t115 * t85 - t118 * t84;
t21 = qJDD(3) * pkin(7) + qJD(3) * t219 + t146;
t97 = t112 * pkin(2) + pkin(1);
t218 = -pkin(7) * t83 - t97;
t167 = t112 * qJDD(1);
t168 = t111 * qJDD(1);
t148 = t115 * t168 - t118 * t167;
t81 = t83 * qJD(3);
t52 = qJD(1) * t81 + t148;
t23 = t52 * pkin(3) - pkin(7) * t131 + t218 * qJDD(1) + qJDD(2);
t89 = -t97 * qJD(1) + qJD(2);
t33 = -pkin(3) * t226 - pkin(7) * t215 + t89;
t48 = qJD(3) * pkin(7) + t54;
t161 = t114 * t21 - t117 * t23 + t48 * t171 + t33 * t172;
t46 = qJDD(4) + t52;
t208 = pkin(4) * t46;
t2 = qJDD(5) + t161 - t208;
t13 = t114 * t33 + t117 * t48;
t9 = qJ(5) * t69 + t13;
t205 = t69 * t9;
t220 = -t2 + t205;
t217 = t214 * t102;
t213 = qJ(2) * qJDD(1);
t47 = -qJD(3) * pkin(3) - t219;
t59 = t114 * t215 - t170;
t14 = pkin(4) * t59 - qJ(5) * t61 + t47;
t207 = pkin(7) * t46;
t212 = t14 * t69 - t207;
t201 = g(3) * t102;
t128 = -t153 * t103 - t201;
t210 = t61 ^ 2;
t209 = t69 ^ 2;
t199 = t13 * t69;
t198 = t59 * t226;
t197 = t61 * t59;
t162 = t61 * t69;
t196 = t61 * t215;
t195 = t215 * t59;
t193 = -t114 * t27 - t59 * t171;
t49 = pkin(3) * t215 - pkin(7) * t226;
t192 = t114 * t49 + t117 * t219;
t51 = -pkin(3) * t223 + t218;
t58 = -t115 * t91 + t118 * t92;
t191 = t114 * t51 + t117 * t58;
t147 = pkin(4) * t114 - qJ(5) * t117;
t190 = -t114 * qJD(5) + t69 * t147 - t54;
t178 = t117 * t119;
t179 = t116 * t117;
t189 = (g(1) * t178 + g(2) * t179) * t102;
t188 = pkin(7) * qJD(4);
t38 = t114 * t46;
t187 = t114 * t69;
t39 = t117 * t46;
t185 = t26 * t114;
t184 = t46 * qJ(5);
t183 = qJD(4) * t83;
t180 = t114 * t116;
t176 = t119 * t114;
t12 = -t114 * t48 + t117 * t33;
t175 = qJD(5) - t12;
t174 = t111 ^ 2 + t112 ^ 2;
t166 = t69 * t188;
t165 = t83 * t172;
t164 = t83 * t171;
t159 = t174 * qJD(1) ^ 2;
t158 = 0.2e1 * t174;
t70 = t103 * t180 + t178;
t72 = t103 * t176 - t179;
t157 = -g(1) * t70 + g(2) * t72;
t71 = t103 * t179 - t176;
t73 = t103 * t178 + t180;
t156 = g(1) * t71 - g(2) * t73;
t155 = t14 * t80 + t3 * t83;
t136 = t114 * t23 + t117 * t21 + t33 * t171 - t48 * t172;
t1 = qJD(5) * t69 + t136 + t184;
t8 = -pkin(4) * t69 + t175;
t154 = t69 * t8 + t1;
t151 = -t114 * t9 + t117 * t8;
t150 = t22 * t83 + t47 * t80;
t149 = t46 * t83 + t69 * t80;
t57 = t115 * t92 + t118 * t91;
t144 = t38 + (-t117 * t226 + t171) * t69;
t142 = -t69 * t172 + t187 * t226 + t39;
t141 = pkin(3) * t103 + t102 * pkin(7) + t97;
t140 = t166 + t200;
t139 = pkin(4) * t117 + qJ(5) * t114 + pkin(3);
t137 = t69 * t47 - t207;
t34 = qJD(2) * t223 - qJD(3) * t57;
t50 = pkin(3) * t81 - pkin(7) * t80;
t135 = t114 * t50 + t117 * t34 + t51 * t171 - t58 * t172;
t130 = t138 + t182;
t129 = g(1) * t72 + g(2) * t70 + t114 * t201 - t161;
t126 = t158 * t169 - t153;
t125 = t14 * t61 + qJDD(5) - t129;
t123 = -g(1) * t73 - g(2) * t71 - t117 * t201 + t136;
t35 = qJD(2) * t83 + qJD(3) * t58;
t88 = -t97 * qJDD(1) + qJDD(2);
t29 = pkin(4) * t61 + qJ(5) * t59;
t28 = t147 * t83 + t57;
t20 = pkin(4) * t223 + t114 * t58 - t117 * t51;
t19 = -qJ(5) * t223 + t191;
t11 = -pkin(4) * t215 + t114 * t219 - t117 * t49;
t10 = qJ(5) * t215 + t192;
t7 = t59 * t69 - t26;
t6 = (pkin(4) * t80 + qJ(5) * t183) * t114 + (-qJ(5) * t80 + (pkin(4) * qJD(4) - qJD(5)) * t83) * t117 + t35;
t5 = -t81 * pkin(4) + t191 * qJD(4) + t114 * t34 - t117 * t50;
t4 = qJ(5) * t81 - qJD(5) * t223 + t135;
t15 = [qJDD(1), t214, t153, t130 * t112, -t130 * t111, t158 * t213 + t126, t138 * pkin(1) + (t174 * t213 + t126) * qJ(2), t122 * t83 + t215 * t80, t122 * t223 - t215 * t81 + t226 * t80 - t83 * t52, qJD(3) * t80 + qJDD(3) * t83, -qJD(3) * t81 + qJDD(3) * t223, 0, -t35 * qJD(3) - t57 * qJDD(3) + t103 * t214 - t223 * t88 - t97 * t52 + t89 * t81, -t58 * qJDD(3) + t89 * t80 + t88 * t83 - t217 - t97 * t132 + (-t226 * t97 - t34) * qJD(3), -t61 * t165 + (-t26 * t83 + t61 * t80) * t117, (-t114 * t61 - t117 * t59) * t80 + (t185 - t117 * t27 + (t114 * t59 - t117 * t61) * qJD(4)) * t83, t117 * t149 - t165 * t69 + t223 * t26 + t61 * t81, -t114 * t149 - t164 * t69 + t223 * t27 - t59 * t81, -t223 * t46 + t69 * t81, t161 * t223 + t12 * t81 + t35 * t59 + t57 * t27 + ((-qJD(4) * t58 + t50) * t69 + t51 * t46 + t47 * t183) * t117 + ((-qJD(4) * t51 - t34) * t69 - t58 * t46 + t150) * t114 + t156, t150 * t117 - t13 * t81 - t135 * t69 + t136 * t223 - t47 * t165 - t191 * t46 - t57 * t26 + t35 * t61 + t157, t114 * t155 + t14 * t164 + t2 * t223 - t20 * t46 + t28 * t27 - t5 * t69 + t6 * t59 - t8 * t81 + t156, -t19 * t27 - t20 * t26 - t4 * t59 + t5 * t61 + t151 * t80 + t217 + (-t1 * t114 + t2 * t117 + (-t114 * t8 - t117 * t9) * qJD(4)) * t83, -t1 * t223 - t117 * t155 + t14 * t165 + t19 * t46 + t28 * t26 + t4 * t69 - t6 * t61 + t9 * t81 - t157, t1 * t19 + t9 * t4 + t3 * t28 + t14 * t6 + t2 * t20 + t8 * t5 - g(1) * (-t71 * pkin(4) - t70 * qJ(5)) - g(2) * (t73 * pkin(4) + t72 * qJ(5)) + (-g(1) * t194 - g(2) * t141) * t119 + (g(1) * t141 - g(2) * t194) * t116; 0, 0, 0, -t167, t168, -t159, -qJ(2) * t159 - t138, 0, 0, 0, 0, 0, t148 + 0.2e1 * t222, 0.2e1 * t226 * qJD(3) + t132, 0, 0, 0, 0, 0, t142 - t195, -t117 * t209 - t196 - t38, -t187 * t69 - t195 + t39, (t26 + t198) * t117 + t114 * t162 + t193, t144 + t196, t154 * t114 + t220 * t117 - t14 * t215 - t214; 0, 0, 0, 0, 0, 0, 0, -t215 * t226, t215 ^ 2 - t226 ^ 2, t132, -t148, qJDD(3), -t215 * t89 + t127 + t143 + t224, -t226 * t89 - t128 - t146, t117 * t162 - t185, (-t26 + t198) * t117 - t61 * t187 + t193, t144 - t196, t142 + t195, -t69 * t215, -pkin(3) * t27 - t12 * t215 - t54 * t59 + (-t200 - t22 + (-t49 - t188) * t69) * t117 + (t219 * t69 + t137) * t114 + t189, pkin(3) * t26 + t192 * t69 + t13 * t215 - t54 * t61 + t137 * t117 + (t140 + t22 - t216) * t114, t11 * t69 - t139 * t27 + t8 * t215 + t190 * t59 + (-t140 - t3) * t117 + t212 * t114 + t189, t10 * t59 - t11 * t61 + ((qJD(4) * t61 - t27) * pkin(7) + t154) * t117 + ((qJD(4) * t59 - t26) * pkin(7) - t220) * t114 + t128, -t10 * t69 - t26 * t139 - t215 * t9 - t190 * t61 - t212 * t117 + (-t166 + t225) * t114, -t9 * t10 - t8 * t11 + t190 * t14 + (qJD(4) * t151 + t1 * t117 + t2 * t114 + t128) * pkin(7) + t225 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t59 ^ 2 + t210, t7, t162 - t27, t46, -t47 * t61 + t129 + t199, t12 * t69 + t47 * t59 - t123, -t29 * t59 - t125 + t199 + 0.2e1 * t208, pkin(4) * t26 - t27 * qJ(5) + (-t13 + t9) * t61 + (t8 - t175) * t59, 0.2e1 * t184 - t14 * t59 + t29 * t61 + (0.2e1 * qJD(5) - t12) * t69 + t123, t1 * qJ(5) - t2 * pkin(4) - t14 * t29 - t8 * t13 - g(1) * (-pkin(4) * t72 + qJ(5) * t73) - g(2) * (-pkin(4) * t70 + qJ(5) * t71) + t175 * t9 + t147 * t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t148 + t197 - t222, t7, -t209 - t210, t125 - t205 - t208;];
tau_reg = t15;
