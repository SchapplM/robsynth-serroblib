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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:49:02
% EndTime: 2022-01-20 09:49:09
% DurationCPUTime: 2.00s
% Computational Cost: add. (3465->282), mult. (6082->356), div. (0->0), fcn. (3792->16), ass. (0->189)
t146 = qJDD(1) + qJDD(3);
t158 = sin(qJ(3));
t162 = cos(qJ(3));
t155 = cos(pkin(9));
t129 = pkin(1) * t155 + pkin(2);
t108 = t129 * qJD(1);
t154 = sin(pkin(9));
t239 = pkin(1) * qJDD(1);
t211 = t154 * t239;
t261 = qJD(3) * t108 + t211;
t258 = pkin(1) * t154;
t214 = qJD(1) * t258;
t263 = -qJD(3) * t214 + t129 * qJDD(1);
t196 = t263 * t158 + t261 * t162;
t43 = pkin(7) * t146 + t196;
t265 = qJD(2) * qJD(4) + t43;
t149 = qJ(1) + pkin(9);
t138 = qJ(3) + t149;
t128 = cos(t138);
t254 = g(2) * t128;
t195 = t261 * t158 - t263 * t162;
t257 = pkin(3) * t146;
t44 = t195 - t257;
t264 = t44 + t254;
t160 = cos(qJ(5));
t161 = cos(qJ(4));
t218 = qJD(5) * t160;
t220 = qJD(4) * t161;
t262 = -t160 * t220 - t161 * t218;
t148 = qJD(1) + qJD(3);
t156 = sin(qJ(5));
t157 = sin(qJ(4));
t230 = t157 * t160;
t86 = t156 * t161 + t230;
t74 = t86 * t148;
t212 = t157 * qJDD(2) + t265 * t161;
t221 = qJD(4) * t157;
t69 = t108 * t158 + t162 * t214;
t61 = pkin(7) * t148 + t69;
t16 = -t61 * t221 + t212;
t137 = t161 * qJDD(2);
t217 = t157 * qJD(2);
t241 = t161 * t61;
t53 = t217 + t241;
t17 = -qJD(4) * t53 - t157 * t43 + t137;
t139 = t161 * qJD(2);
t243 = t157 * t61;
t52 = t139 - t243;
t170 = -t17 * t157 + t16 * t161 + (-t157 * t53 - t161 * t52) * qJD(4);
t127 = sin(t138);
t260 = g(1) * t128 + g(2) * t127;
t79 = t129 * t162 - t158 * t258;
t80 = t158 * t129 + t162 * t258;
t147 = qJD(4) + qJD(5);
t204 = pkin(8) * t148 + t61;
t51 = t204 * t161 + t217;
t13 = qJDD(4) * pkin(4) + t137 + (-pkin(8) * t146 - t43) * t157 - t51 * qJD(4);
t207 = t148 * t221;
t228 = t161 * t146;
t178 = t207 - t228;
t14 = -pkin(8) * t178 + t16;
t219 = qJD(5) * t156;
t50 = -t204 * t157 + t139;
t49 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t49 + t14) * t160 + t156 * t13 - t51 * t219;
t164 = -pkin(8) - pkin(7);
t78 = pkin(7) + t80;
t259 = -pkin(8) - t78;
t256 = pkin(3) * t148;
t255 = pkin(4) * t161;
t120 = g(1) * t127;
t253 = g(3) * t161;
t232 = t156 * t157;
t210 = t148 * t232;
t229 = t160 * t161;
t72 = -t148 * t229 + t210;
t252 = t74 * t72;
t111 = t164 * t157;
t142 = t161 * pkin(8);
t112 = pkin(7) * t161 + t142;
t65 = t111 * t160 - t112 * t156;
t68 = t162 * t108 - t158 * t214;
t85 = -t229 + t232;
t208 = qJD(4) * t164;
t94 = t157 * t208;
t95 = t161 * t208;
t251 = qJD(5) * t65 + t156 * t95 + t160 * t94 + t85 * t68;
t66 = t111 * t156 + t112 * t160;
t250 = -qJD(5) * t66 - t156 * t94 + t160 * t95 + t86 * t68;
t231 = t157 * t146;
t187 = t156 * t231 - t160 * t228;
t58 = t147 * t86;
t31 = t148 * t58 + t187;
t185 = t147 * t232;
t57 = t185 + t262;
t249 = -t86 * t31 + t57 * t72;
t248 = t148 * t68;
t247 = t148 * t69;
t70 = t79 * qJD(3);
t246 = t148 * t70;
t71 = t80 * qJD(3);
t245 = t148 * t71;
t244 = t156 * t51;
t242 = t160 * t51;
t60 = -t68 - t256;
t240 = t161 * t120 + t60 * t221;
t153 = qJ(4) + qJ(5);
t140 = sin(t153);
t238 = t127 * t140;
t141 = cos(t153);
t237 = t127 * t141;
t236 = t128 * t140;
t235 = t128 * t141;
t233 = t148 * t157;
t227 = t128 * pkin(3) + t127 * pkin(7);
t135 = cos(t149);
t163 = cos(qJ(1));
t225 = t163 * pkin(1) + pkin(2) * t135;
t150 = t157 ^ 2;
t151 = t161 ^ 2;
t224 = t150 - t151;
t223 = t150 + t151;
t215 = t264 * t157 + t60 * t220;
t213 = pkin(4) * t221;
t144 = t148 ^ 2;
t209 = t157 * t144 * t161;
t133 = pkin(3) + t255;
t202 = qJD(4) * t259;
t201 = -pkin(3) * t127 + t128 * pkin(7);
t198 = -t127 * t164 + t128 * t133;
t197 = t223 * t146;
t194 = -t146 * t230 + t262 * t148 - t156 * t228;
t192 = t161 * t207;
t191 = -t69 + t213;
t77 = -pkin(3) - t79;
t134 = sin(t149);
t159 = sin(qJ(1));
t190 = -pkin(1) * t159 - pkin(2) * t134;
t188 = g(1) * t159 - g(2) * t163;
t30 = t148 * t185 + t194;
t186 = -t30 * t85 + t58 * t74;
t18 = t160 * t49 - t244;
t19 = t156 * t49 + t242;
t4 = -qJD(5) * t19 + t160 * t13 - t156 * t14;
t184 = t18 * t57 - t19 * t58 - t3 * t85 - t4 * t86 - t260;
t145 = qJDD(4) + qJDD(5);
t45 = t145 * t86 - t147 * t57;
t62 = t259 * t157;
t63 = t161 * t78 + t142;
t33 = -t156 * t63 + t160 * t62;
t34 = t156 * t62 + t160 * t63;
t183 = t157 * t52 - t161 * t53;
t182 = -t127 * t133 - t128 * t164;
t181 = -t196 + t260;
t27 = pkin(4) * t178 + t44;
t54 = -t133 * t148 - t68;
t180 = -g(1) * t238 + g(2) * t236 + t27 * t86 - t54 * t57;
t179 = g(1) * t237 - g(2) * t235 + t27 * t85 + t54 * t58;
t165 = qJD(4) ^ 2;
t177 = pkin(7) * t165 - t247 - t257;
t176 = -t148 * t60 + t260;
t175 = t120 - t195 - t254;
t174 = t146 * t77 + t165 * t78 + t245;
t173 = -pkin(7) * qJDD(4) + (t68 - t256) * qJD(4);
t172 = -qJDD(4) * t78 + (t148 * t77 - t70) * qJD(4);
t169 = -t260 + t170;
t168 = g(1) * t235 + g(2) * t237 + g(3) * t140 + t54 * t72 - t3;
t167 = g(1) * t236 + g(2) * t238 - g(3) * t141 - t54 * t74 + t4;
t152 = qJDD(2) - g(3);
t99 = qJDD(4) * t161 - t157 * t165;
t98 = qJDD(4) * t157 + t161 * t165;
t76 = t146 * t151 - 0.2e1 * t192;
t75 = t146 * t150 + 0.2e1 * t192;
t67 = t77 - t255;
t64 = t71 + t213;
t59 = -0.2e1 * t224 * t148 * qJD(4) + 0.2e1 * t157 * t228;
t46 = -t145 * t85 - t147 * t58;
t42 = -t157 * t70 + t161 * t202;
t41 = t157 * t202 + t161 * t70;
t32 = -t72 ^ 2 + t74 ^ 2;
t22 = -t194 + (-t210 + t72) * t147;
t21 = t160 * t50 - t244;
t20 = -t156 * t50 - t242;
t11 = t31 * t85 + t58 * t72;
t10 = -t30 * t86 - t57 * t74;
t7 = -qJD(5) * t34 - t156 * t41 + t160 * t42;
t6 = qJD(5) * t33 + t156 * t42 + t160 * t41;
t5 = -t186 + t249;
t1 = [0, 0, 0, 0, 0, qJDD(1), t188, g(1) * t163 + g(2) * t159, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t134 - g(2) * t135 + 0.2e1 * t155 * t239, g(1) * t135 + g(2) * t134 - 0.2e1 * t211, 0, (t188 + (t154 ^ 2 + t155 ^ 2) * t239) * pkin(1), 0, 0, 0, 0, 0, t146, t146 * t79 + t175 - t245, -t146 * t80 + t181 - t246, 0, -g(1) * t190 - g(2) * t225 - t195 * t79 + t196 * t80 - t68 * t71 + t69 * t70, t75, t59, t98, t76, t99, 0, t172 * t157 + (-t174 - t264) * t161 + t240, t172 * t161 + (t174 - t120) * t157 + t215, t197 * t78 + t223 * t246 + t169, t44 * t77 + t60 * t71 - g(1) * (t190 + t201) - g(2) * (t225 + t227) - t183 * t70 + t170 * t78, t10, t5, t45, t11, t46, 0, t145 * t33 + t147 * t7 + t31 * t67 + t64 * t72 + t179, -t145 * t34 - t147 * t6 - t30 * t67 + t64 * t74 + t180, t30 * t33 - t31 * t34 - t6 * t72 - t7 * t74 + t184, t3 * t34 + t19 * t6 + t4 * t33 + t18 * t7 + t27 * t67 + t54 * t64 - g(1) * (t182 + t190) - g(2) * (t198 + t225); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, 0, 0, 0, 0, 0, t99, -t98, 0, -qJD(4) * t183 + t157 * t16 + t161 * t17 - g(3), 0, 0, 0, 0, 0, 0, t46, -t45, t186 + t249, -t18 * t58 - t19 * t57 + t3 * t86 - t4 * t85 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t175 + t247, t181 + t248, 0, 0, t75, t59, t98, t76, t99, 0, t173 * t157 + (-t177 - t264) * t161 + t240, t173 * t161 + (t177 - t120) * t157 + t215, pkin(7) * t197 - t223 * t248 + t169, -t44 * pkin(3) + pkin(7) * t170 - g(1) * t201 - g(2) * t227 + t183 * t68 - t60 * t69, t10, t5, t45, t11, t46, 0, -t133 * t31 + t145 * t65 + t250 * t147 + t191 * t72 + t179, t133 * t30 - t145 * t66 - t251 * t147 + t191 * t74 + t180, -t250 * t74 - t251 * t72 + t30 * t65 - t31 * t66 + t184, -g(1) * t182 - g(2) * t198 - t27 * t133 + t250 * t18 + t251 * t19 + t191 * t54 + t3 * t66 + t4 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t224 * t144, t231, t209, t228, qJDD(4), -t253 + t137 + (t53 - t241) * qJD(4) + (t176 - t265) * t157, g(3) * t157 + (t52 + t243) * qJD(4) + t176 * t161 - t212, 0, 0, t252, t32, t22, -t252, -t187, t145, -t147 * t20 + (t145 * t160 - t147 * t219 - t233 * t72) * pkin(4) + t167, t147 * t21 + (-t145 * t156 - t147 * t218 - t233 * t74) * pkin(4) + t168, (t19 + t20) * t74 + (-t18 + t21) * t72 + (-t156 * t31 + t160 * t30 + (t156 * t74 - t160 * t72) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t253 + t156 * t3 + t160 * t4 + (-t156 * t18 + t160 * t19) * qJD(5) + (-t148 * t54 + t260) * t157) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t32, t22, -t252, -t187, t145, t147 * t19 + t167, t147 * t18 + t168, 0, 0;];
tau_reg = t1;
