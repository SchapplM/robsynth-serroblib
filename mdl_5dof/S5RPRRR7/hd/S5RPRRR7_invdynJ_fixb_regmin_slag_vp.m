% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:12
% EndTime: 2019-12-31 19:04:18
% DurationCPUTime: 2.55s
% Computational Cost: add. (2171->317), mult. (4596->455), div. (0->0), fcn. (3202->14), ass. (0->171)
t140 = cos(qJ(3));
t194 = t140 * qJD(1);
t112 = -qJD(4) + t194;
t136 = sin(qJ(3));
t128 = qJ(1) + pkin(9);
t119 = sin(t128);
t120 = cos(t128);
t165 = g(1) * t120 + g(2) * t119;
t149 = -g(3) * t140 + t136 * t165;
t132 = sin(pkin(9));
t113 = pkin(1) * t132 + pkin(6);
t103 = t113 * qJD(1);
t202 = qJD(3) * t140;
t99 = t113 * qJDD(1);
t240 = -qJD(2) * qJD(3) - t99;
t188 = -t103 * t202 + t136 * t240;
t24 = -qJDD(3) * pkin(3) - t140 * qJDD(2) - t188;
t242 = qJD(4) * pkin(7) * t112 + t149 - t24;
t134 = sin(qJ(5));
t138 = cos(qJ(5));
t135 = sin(qJ(4));
t139 = cos(qJ(4));
t195 = t139 * qJD(3);
t204 = qJD(1) * t136;
t87 = t135 * t204 - t195;
t196 = t135 * qJD(3);
t89 = t139 * t204 + t196;
t162 = t134 * t87 - t138 * t89;
t39 = t134 * t89 + t138 * t87;
t241 = t162 * t39;
t64 = qJD(2) * t140 - t103 * t136;
t239 = t64 * qJD(3);
t200 = qJD(4) * t136;
t238 = -qJD(1) * t200 + qJDD(3);
t237 = t162 ^ 2 - t39 ^ 2;
t110 = -qJD(5) + t112;
t197 = qJD(5) * t138;
t198 = qJD(5) * t134;
t179 = t140 * t195;
t190 = t136 * qJDD(1);
t33 = qJD(1) * t179 + qJD(4) * t195 + t135 * t238 + t139 * t190;
t34 = ((qJD(4) + t194) * qJD(3) + t190) * t135 - t238 * t139;
t6 = -t134 * t34 + t138 * t33 - t197 * t87 - t198 * t89;
t236 = -t110 * t39 + t6;
t131 = qJ(4) + qJ(5);
t127 = cos(t131);
t65 = qJD(2) * t136 + t103 * t140;
t55 = qJD(3) * pkin(7) + t65;
t133 = cos(pkin(9));
t114 = -pkin(1) * t133 - pkin(2);
t80 = -pkin(3) * t140 - pkin(7) * t136 + t114;
t56 = t80 * qJD(1);
t20 = t135 * t56 + t139 * t55;
t15 = -pkin(8) * t87 + t20;
t13 = t15 * t198;
t227 = g(3) * t136;
t54 = -qJD(3) * pkin(3) - t64;
t35 = pkin(4) * t87 + t54;
t126 = sin(t131);
t211 = t127 * t140;
t50 = -t119 * t211 + t120 * t126;
t52 = t119 * t126 + t120 * t211;
t235 = g(1) * t52 - g(2) * t50 + t127 * t227 + t35 * t39 + t13;
t23 = qJDD(3) * pkin(7) + t136 * qJDD(2) + t140 * t99 + t239;
t166 = pkin(3) * t136 - pkin(7) * t140;
t96 = t166 * qJD(3);
t36 = qJD(1) * t96 + qJDD(1) * t80;
t30 = t139 * t36;
t123 = t140 * qJDD(1);
t192 = qJD(1) * qJD(3);
t176 = t136 * t192;
t83 = qJDD(4) - t123 + t176;
t2 = t83 * pkin(4) - t33 * pkin(8) - qJD(4) * t20 - t135 * t23 + t30;
t199 = qJD(4) * t139;
t189 = t135 * t36 + t139 * t23 + t199 * t56;
t201 = qJD(4) * t135;
t156 = -t201 * t55 + t189;
t3 = -pkin(8) * t34 + t156;
t185 = -t134 * t3 + t138 * t2;
t212 = t126 * t140;
t49 = t119 * t212 + t120 * t127;
t19 = -t135 * t55 + t139 * t56;
t14 = -pkin(8) * t89 + t19;
t11 = -pkin(4) * t112 + t14;
t218 = t138 * t15;
t5 = t11 * t134 + t218;
t51 = t119 * t127 - t120 * t212;
t234 = -g(1) * t51 + g(2) * t49 - qJD(5) * t5 + t126 * t227 + t35 * t162 + t185;
t147 = qJD(5) * t162 - t134 * t33 - t138 * t34;
t233 = t110 * t162 + t147;
t91 = t134 * t139 + t135 * t138;
t66 = t91 * t136;
t153 = -t135 * t200 + t179;
t231 = qJD(4) + qJD(5);
t208 = t136 * t139;
t230 = -t112 * t153 + t208 * t83;
t228 = pkin(7) + pkin(8);
t180 = t140 * t196;
t210 = t135 * t136;
t17 = -t198 * t210 + (t208 * t231 + t180) * t138 + t153 * t134;
t79 = qJDD(5) + t83;
t225 = t110 * t17 - t66 * t79;
t90 = t134 * t135 - t138 * t139;
t155 = t90 * t140;
t224 = qJD(1) * t155 - t231 * t90;
t223 = (-t194 + t231) * t91;
t93 = t166 * qJD(1);
t222 = t135 * t93 + t139 * t64;
t221 = t135 * t96 + t199 * t80;
t220 = t113 * t136 * t196 + t139 * t96;
t207 = t139 * t140;
t92 = t113 * t207;
t219 = t135 * t80 + t92;
t217 = t33 * t135;
t216 = t87 * t112;
t215 = t89 * t112;
t214 = qJD(3) * t87;
t213 = t112 * t139;
t209 = t135 * t140;
t206 = qJDD(2) - g(3);
t129 = t136 ^ 2;
t205 = -t140 ^ 2 + t129;
t104 = qJD(1) * t114;
t203 = qJD(3) * t136;
t184 = qJD(4) * t228;
t183 = t135 * t194;
t182 = t112 * t196;
t178 = -t140 * t6 - t162 * t203;
t175 = qJD(5) * t11 + t3;
t173 = -t33 * t140 + t203 * t89;
t172 = -qJD(4) * t56 - t23;
t171 = t112 * t113 + t55;
t169 = -t65 + (-t183 + t201) * pkin(4);
t105 = t228 * t135;
t168 = pkin(8) * t183 - qJD(5) * t105 - t135 * t184 - t222;
t106 = t228 * t139;
t161 = pkin(4) * t136 - pkin(8) * t207;
t77 = t139 * t93;
t167 = qJD(1) * t161 + qJD(5) * t106 - t135 * t64 + t139 * t184 + t77;
t137 = sin(qJ(1));
t141 = cos(qJ(1));
t164 = g(1) * t137 - g(2) * t141;
t16 = -qJD(3) * t155 - t231 * t66;
t67 = t90 * t136;
t163 = -t110 * t16 - t67 * t79;
t159 = -t140 * t147 - t203 * t39;
t158 = t112 * t199 - t135 * t83;
t154 = -qJD(1) * t104 + t165;
t152 = t136 * t199 + t180;
t151 = -pkin(7) * t83 - t112 * t54;
t150 = 0.2e1 * qJD(3) * t104 - qJDD(3) * t113;
t142 = qJD(3) ^ 2;
t146 = g(1) * t119 - g(2) * t120 - 0.2e1 * qJDD(1) * t114 - t113 * t142;
t143 = qJD(1) ^ 2;
t118 = -pkin(4) * t139 - pkin(3);
t98 = qJDD(3) * t140 - t136 * t142;
t97 = qJDD(3) * t136 + t140 * t142;
t72 = (pkin(4) * t135 + t113) * t136;
t69 = t139 * t80;
t62 = t119 * t135 + t120 * t207;
t61 = t119 * t139 - t120 * t209;
t60 = -t119 * t207 + t120 * t135;
t59 = t119 * t209 + t120 * t139;
t45 = pkin(4) * t152 + t113 * t202;
t32 = -pkin(8) * t210 + t219;
t25 = -pkin(8) * t208 + t69 + (-t113 * t135 - pkin(4)) * t140;
t10 = (-t136 * t195 - t140 * t201) * t113 - t152 * pkin(8) + t221;
t9 = pkin(4) * t34 + t24;
t8 = t161 * qJD(3) + (-t92 + (pkin(8) * t136 - t80) * t135) * qJD(4) + t220;
t4 = t11 * t138 - t134 * t15;
t1 = [qJDD(1), t164, g(1) * t141 + g(2) * t137, (t164 + (t132 ^ 2 + t133 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t129 + 0.2e1 * t140 * t176, 0.2e1 * t123 * t136 - 0.2e1 * t192 * t205, t97, t98, 0, t136 * t150 + t140 * t146, -t136 * t146 + t140 * t150, t89 * t179 + (t33 * t139 - t201 * t89) * t136, (-t135 * t89 - t139 * t87) * t202 + (-t217 - t139 * t34 + (t135 * t87 - t139 * t89) * qJD(4)) * t136, t173 + t230, (t34 + t182) * t140 + (t158 - t214) * t136, -t112 * t203 - t140 * t83, -(-t201 * t80 + t220) * t112 + t69 * t83 - g(1) * t60 - g(2) * t62 + (t113 * t214 - t30 + t171 * t199 + (qJD(3) * t54 - t113 * t83 - t172) * t135) * t140 + (qJD(3) * t19 + t113 * t34 + t24 * t135 + t199 * t54) * t136, t221 * t112 - t219 * t83 - g(1) * t59 - g(2) * t61 + (-t171 * t201 + (t113 * t89 + t139 * t54) * qJD(3) + t189) * t140 + (-t54 * t201 + t113 * t33 + t24 * t139 + (-t113 * t213 - t20) * qJD(3)) * t136, -t16 * t162 - t6 * t67, -t147 * t67 - t16 * t39 + t162 * t17 - t6 * t66, t163 + t178, t159 + t225, -t110 * t203 - t140 * t79, -(-t134 * t10 + t138 * t8) * t110 + (-t134 * t32 + t138 * t25) * t79 - t185 * t140 + t4 * t203 + t45 * t39 - t72 * t147 + t9 * t66 + t35 * t17 - g(1) * t50 - g(2) * t52 + (-(-t134 * t25 - t138 * t32) * t110 + t5 * t140) * qJD(5), -t5 * t203 - g(1) * t49 - g(2) * t51 - t13 * t140 + t35 * t16 - t45 * t162 + t72 * t6 - t9 * t67 + ((-qJD(5) * t32 + t8) * t110 - t25 * t79 + t2 * t140) * t134 + ((qJD(5) * t25 + t10) * t110 - t32 * t79 + t175 * t140) * t138; 0, 0, 0, t206, 0, 0, 0, 0, 0, t98, -t97, 0, 0, 0, 0, 0, (-t34 + t182) * t140 + (t158 + t214) * t136, t173 - t230, 0, 0, 0, 0, 0, -t159 + t225, -t163 + t178; 0, 0, 0, 0, -t136 * t143 * t140, t205 * t143, t190, t123, qJDD(3), t65 * qJD(3) + t136 * t154 + t140 * t206 + t188, t239 + (qJD(3) * t103 - t206) * t136 + (t154 + t240) * t140, -t213 * t89 + t217, (t33 + t216) * t139 + (-t34 + t215) * t135, (t112 * t207 - t136 * t89) * qJD(1) - t158, t112 * t201 + t139 * t83 + (-t112 * t209 + t136 * t87) * qJD(1), t112 * t204, -t19 * t204 - pkin(3) * t34 + t77 * t112 - t65 * t87 + (-t112 * t64 + t151) * t135 + t242 * t139, -pkin(3) * t33 - t222 * t112 - t135 * t242 + t151 * t139 + t20 * t204 - t65 * t89, -t162 * t224 + t6 * t91, t147 * t91 + t162 * t223 - t224 * t39 - t6 * t90, -t110 * t224 + t162 * t204 + t91 * t79, t110 * t223 + t204 * t39 - t90 * t79, t110 * t204, (-t105 * t138 - t106 * t134) * t79 - t118 * t147 + t9 * t90 - t4 * t204 + t169 * t39 + t223 * t35 + (t134 * t168 + t138 * t167) * t110 + t149 * t127, -(-t105 * t134 + t106 * t138) * t79 + t118 * t6 + t9 * t91 + t5 * t204 - t169 * t162 + t224 * t35 + (-t134 * t167 + t138 * t168) * t110 - t149 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t87, -t87 ^ 2 + t89 ^ 2, t33 - t216, -t215 - t34, t83, -t55 * t199 - g(1) * t61 + g(2) * t59 - t20 * t112 - t54 * t89 + t30 + (t172 + t227) * t135, g(1) * t62 - g(2) * t60 + g(3) * t208 - t112 * t19 + t54 * t87 - t156, -t241, t237, t236, t233, t79, (-t134 * t14 - t218) * t110 + (t110 * t198 + t138 * t79 - t39 * t89) * pkin(4) + t234, (t110 * t15 - t2) * t134 + (-t110 * t14 - t175) * t138 + (t110 * t197 - t134 * t79 + t162 * t89) * pkin(4) + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t237, t236, t233, t79, -t5 * t110 + t234, -t4 * t110 - t134 * t2 - t138 * t175 + t235;];
tau_reg = t1;
