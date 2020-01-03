% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR5
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:54:15
% EndTime: 2020-01-03 11:54:18
% DurationCPUTime: 1.72s
% Computational Cost: add. (3465->279), mult. (6082->351), div. (0->0), fcn. (3792->16), ass. (0->187)
t150 = qJ(1) + pkin(9);
t138 = qJ(3) + t150;
t127 = sin(t138);
t128 = cos(t138);
t190 = -g(2) * t128 - g(3) * t127;
t159 = sin(qJ(3));
t163 = cos(qJ(3));
t156 = cos(pkin(9));
t129 = t156 * pkin(1) + pkin(2);
t107 = t129 * qJD(1);
t155 = sin(pkin(9));
t237 = pkin(1) * qJDD(1);
t211 = t155 * t237;
t260 = qJD(3) * t107 + t211;
t257 = pkin(1) * t155;
t214 = qJD(1) * t257;
t262 = -qJD(3) * t214 + t129 * qJDD(1);
t195 = t260 * t159 - t262 * t163;
t147 = qJDD(1) + qJDD(3);
t253 = t147 * pkin(3);
t44 = t195 - t253;
t182 = t190 - t44;
t196 = t262 * t159 + t260 * t163;
t43 = t147 * pkin(7) + t196;
t263 = qJD(2) * qJD(4) + t43;
t119 = g(3) * t128;
t259 = g(2) * t127 - t119;
t161 = cos(qJ(5));
t162 = cos(qJ(4));
t217 = qJD(5) * t161;
t219 = qJD(4) * t162;
t261 = -t161 * t219 - t162 * t217;
t149 = qJD(1) + qJD(3);
t157 = sin(qJ(5));
t158 = sin(qJ(4));
t231 = t161 * t158;
t86 = t157 * t162 + t231;
t74 = t86 * t149;
t212 = t158 * qJDD(2) + t263 * t162;
t220 = qJD(4) * t158;
t69 = t159 * t107 + t163 * t214;
t61 = t149 * pkin(7) + t69;
t16 = -t61 * t220 + t212;
t137 = t162 * qJDD(2);
t216 = t158 * qJD(2);
t243 = t162 * t61;
t53 = t216 + t243;
t17 = -t53 * qJD(4) - t158 * t43 + t137;
t139 = t162 * qJD(2);
t245 = t158 * t61;
t52 = t139 - t245;
t171 = -t17 * t158 + t16 * t162 + (-t158 * t53 - t162 * t52) * qJD(4);
t79 = t163 * t129 - t159 * t257;
t80 = t159 * t129 + t163 * t257;
t148 = qJD(4) + qJD(5);
t204 = pkin(8) * t149 + t61;
t51 = t204 * t162 + t216;
t13 = qJDD(4) * pkin(4) + t137 + (-pkin(8) * t147 - t43) * t158 - t51 * qJD(4);
t207 = t149 * t220;
t229 = t162 * t147;
t180 = t207 - t229;
t14 = -pkin(8) * t180 + t16;
t218 = qJD(5) * t157;
t50 = -t204 * t158 + t139;
t49 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t49 + t14) * t161 + t157 * t13 - t51 * t218;
t165 = -pkin(8) - pkin(7);
t78 = pkin(7) + t80;
t258 = -pkin(8) - t78;
t256 = g(1) * t162;
t252 = t149 * pkin(3);
t251 = t162 * pkin(4);
t233 = t157 * t158;
t210 = t149 * t233;
t230 = t161 * t162;
t72 = -t149 * t230 + t210;
t250 = t74 * t72;
t111 = t165 * t158;
t143 = t162 * pkin(8);
t112 = t162 * pkin(7) + t143;
t65 = t161 * t111 - t157 * t112;
t68 = t163 * t107 - t159 * t214;
t85 = -t230 + t233;
t208 = qJD(4) * t165;
t93 = t158 * t208;
t94 = t162 * t208;
t249 = t65 * qJD(5) + t157 * t94 + t161 * t93 + t85 * t68;
t66 = t157 * t111 + t161 * t112;
t248 = -t66 * qJD(5) - t157 * t93 + t161 * t94 + t86 * t68;
t232 = t158 * t147;
t187 = t157 * t232 - t161 * t229;
t58 = t148 * t86;
t31 = t149 * t58 + t187;
t185 = t148 * t233;
t57 = t185 + t261;
t247 = -t86 * t31 + t57 * t72;
t246 = t157 * t51;
t244 = t161 * t51;
t242 = t68 * t149;
t241 = t69 * t149;
t70 = t79 * qJD(3);
t240 = t70 * t149;
t71 = t80 * qJD(3);
t239 = t71 * t149;
t133 = pkin(3) + t251;
t238 = t127 * t133 + t128 * t165;
t154 = qJ(4) + qJ(5);
t140 = sin(t154);
t236 = t127 * t140;
t235 = t128 * t140;
t234 = t149 * t158;
t227 = t128 * pkin(3) + t127 * pkin(7);
t134 = sin(t150);
t160 = sin(qJ(1));
t225 = t160 * pkin(1) + pkin(2) * t134;
t135 = cos(t150);
t164 = cos(qJ(1));
t224 = t164 * pkin(1) + pkin(2) * t135;
t151 = t158 ^ 2;
t152 = t162 ^ 2;
t223 = t151 - t152;
t222 = t151 + t152;
t213 = pkin(4) * t220;
t145 = t149 ^ 2;
t209 = t158 * t145 * t162;
t203 = qJD(4) * t258;
t27 = pkin(4) * t180 + t44;
t54 = -t133 * t149 - t68;
t201 = g(2) * t235 + g(3) * t236 + t27 * t86 - t54 * t57;
t199 = -t127 * t165 + t128 * t133;
t198 = t222 * t147;
t60 = -t68 - t252;
t197 = -t182 * t158 + t60 * t219;
t194 = -t147 * t231 + t261 * t149 - t157 * t229;
t192 = t162 * t207;
t191 = -t69 + t213;
t77 = -pkin(3) - t79;
t188 = -g(2) * t164 - g(3) * t160;
t30 = t149 * t185 + t194;
t186 = -t85 * t30 + t74 * t58;
t18 = t161 * t49 - t246;
t19 = t157 * t49 + t244;
t4 = -t19 * qJD(5) + t161 * t13 - t157 * t14;
t184 = t18 * t57 - t19 * t58 - t3 * t85 - t4 * t86 - t259;
t146 = qJDD(4) + qJDD(5);
t45 = t86 * t146 - t57 * t148;
t62 = t258 * t158;
t63 = t162 * t78 + t143;
t33 = -t157 * t63 + t161 * t62;
t34 = t157 * t62 + t161 * t63;
t183 = t52 * t158 - t53 * t162;
t181 = -t196 + t259;
t166 = qJD(4) ^ 2;
t179 = pkin(7) * t166 - t241 - t253;
t178 = -t149 * t60 + t259;
t177 = t147 * t77 + t166 * t78 + t239;
t176 = -pkin(7) * qJDD(4) + (t68 - t252) * qJD(4);
t175 = -qJDD(4) * t78 + (t149 * t77 - t70) * qJD(4);
t174 = t190 - t195;
t141 = cos(t154);
t173 = t190 * t141 + t27 * t85 + t54 * t58;
t170 = -t259 + t171;
t169 = g(1) * t140 + t259 * t141 + t54 * t72 - t3;
t168 = -g(1) * t141 + g(2) * t236 - g(3) * t235 - t54 * t74 + t4;
t153 = qJDD(2) - g(1);
t117 = t127 * pkin(3);
t98 = qJDD(4) * t162 - t166 * t158;
t97 = qJDD(4) * t158 + t166 * t162;
t76 = t152 * t147 - 0.2e1 * t192;
t75 = t151 * t147 + 0.2e1 * t192;
t67 = t77 - t251;
t64 = t71 + t213;
t59 = -0.2e1 * t223 * t149 * qJD(4) + 0.2e1 * t158 * t229;
t55 = t60 * t220;
t46 = -t85 * t146 - t58 * t148;
t42 = -t158 * t70 + t162 * t203;
t41 = t158 * t203 + t162 * t70;
t32 = -t72 ^ 2 + t74 ^ 2;
t22 = -t194 + (-t210 + t72) * t148;
t21 = t161 * t50 - t246;
t20 = -t157 * t50 - t244;
t11 = t31 * t85 + t72 * t58;
t10 = -t30 * t86 - t74 * t57;
t7 = -t34 * qJD(5) - t157 * t41 + t161 * t42;
t6 = t33 * qJD(5) + t157 * t42 + t161 * t41;
t5 = -t186 + t247;
t1 = [0, 0, 0, 0, 0, qJDD(1), t188, g(2) * t160 - g(3) * t164, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -g(2) * t135 - g(3) * t134 + 0.2e1 * t156 * t237, g(2) * t134 - g(3) * t135 - 0.2e1 * t211, 0, (t188 + (t155 ^ 2 + t156 ^ 2) * t237) * pkin(1), 0, 0, 0, 0, 0, t147, t79 * t147 + t174 - t239, -t80 * t147 + t181 - t240, 0, -g(2) * t224 - g(3) * t225 - t195 * t79 + t196 * t80 - t68 * t71 + t69 * t70, t75, t59, t97, t76, t98, 0, t55 + t175 * t158 + (-t177 + t182) * t162, t158 * t177 + t162 * t175 + t197, t78 * t198 + t222 * t240 + t170, t44 * t77 + t60 * t71 - g(2) * (t224 + t227) - g(3) * (-t128 * pkin(7) + t117 + t225) - t183 * t70 + t171 * t78, t10, t5, t45, t11, t46, 0, t33 * t146 + t7 * t148 + t67 * t31 + t64 * t72 + t173, -t34 * t146 - t6 * t148 - t67 * t30 + t64 * t74 + t201, t33 * t30 - t34 * t31 - t6 * t72 - t7 * t74 + t184, t3 * t34 + t19 * t6 + t4 * t33 + t18 * t7 + t27 * t67 + t54 * t64 - g(2) * (t199 + t224) - g(3) * (t225 + t238); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, 0, 0, 0, 0, 0, t98, -t97, 0, -t183 * qJD(4) + t16 * t158 + t17 * t162 - g(1), 0, 0, 0, 0, 0, 0, t46, -t45, t186 + t247, -t18 * t58 - t19 * t57 + t3 * t86 - t4 * t85 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t174 + t241, t181 + t242, 0, 0, t75, t59, t97, t76, t98, 0, t55 + t176 * t158 + (-t179 + t182) * t162, t158 * t179 + t162 * t176 + t197, pkin(7) * t198 - t222 * t242 + t170, -t44 * pkin(3) - t60 * t69 - g(2) * t227 - g(3) * t117 + t183 * t68 + (t171 + t119) * pkin(7), t10, t5, t45, t11, t46, 0, -t133 * t31 + t65 * t146 + t248 * t148 + t191 * t72 + t173, t133 * t30 - t66 * t146 - t249 * t148 + t191 * t74 + t201, -t248 * t74 - t249 * t72 + t65 * t30 - t66 * t31 + t184, -g(2) * t199 - g(3) * t238 - t27 * t133 + t248 * t18 + t249 * t19 + t191 * t54 + t3 * t66 + t4 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t223 * t145, t232, t209, t229, qJDD(4), -t256 + t137 + (t53 - t243) * qJD(4) + (t178 - t263) * t158, g(1) * t158 + (t52 + t245) * qJD(4) + t178 * t162 - t212, 0, 0, t250, t32, t22, -t250, -t187, t146, -t20 * t148 + (t146 * t161 - t148 * t218 - t72 * t234) * pkin(4) + t168, t21 * t148 + (-t146 * t157 - t148 * t217 - t74 * t234) * pkin(4) + t169, (t19 + t20) * t74 + (-t18 + t21) * t72 + (-t157 * t31 + t161 * t30 + (t157 * t74 - t161 * t72) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t256 + t157 * t3 + t161 * t4 + (-t157 * t18 + t161 * t19) * qJD(5) + (-t149 * t54 + t259) * t158) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, t32, t22, -t250, -t187, t146, t19 * t148 + t168, t18 * t148 + t169, 0, 0;];
tau_reg = t1;
