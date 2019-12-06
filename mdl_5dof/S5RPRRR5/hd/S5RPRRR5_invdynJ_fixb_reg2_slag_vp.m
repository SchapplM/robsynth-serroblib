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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:16:36
% EndTime: 2019-12-05 18:16:39
% DurationCPUTime: 1.69s
% Computational Cost: add. (3465->281), mult. (6082->351), div. (0->0), fcn. (3792->16), ass. (0->185)
t150 = cos(pkin(9));
t125 = t150 * pkin(1) + pkin(2);
t149 = sin(pkin(9));
t248 = pkin(1) * t149;
t209 = qJD(1) * t248;
t253 = qJD(3) * t209 - t125 * qJDD(1);
t141 = qJDD(1) + qJDD(3);
t153 = sin(qJ(3));
t157 = cos(qJ(3));
t105 = t125 * qJD(1);
t230 = pkin(1) * qJDD(1);
t206 = t149 * t230;
t250 = qJD(3) * t105 + t206;
t192 = -t253 * t153 + t250 * t157;
t43 = t141 * pkin(7) + t192;
t252 = qJD(2) * qJD(4) + t43;
t144 = qJ(1) + pkin(9);
t134 = qJ(3) + t144;
t123 = sin(t134);
t124 = cos(t134);
t184 = g(2) * t124 + g(3) * t123;
t117 = g(2) * t123;
t220 = -g(3) * t124 + t117;
t155 = cos(qJ(5));
t156 = cos(qJ(4));
t213 = qJD(5) * t155;
t215 = qJD(4) * t156;
t251 = -t155 * t215 - t156 * t213;
t143 = qJD(1) + qJD(3);
t151 = sin(qJ(5));
t152 = sin(qJ(4));
t224 = t155 * t152;
t86 = t151 * t156 + t224;
t74 = t86 * t143;
t207 = t152 * qJDD(2) + t252 * t156;
t216 = qJD(4) * t152;
t69 = t153 * t105 + t157 * t209;
t61 = t143 * pkin(7) + t69;
t16 = -t61 * t216 + t207;
t133 = t156 * qJDD(2);
t212 = t152 * qJD(2);
t235 = t156 * t61;
t53 = t212 + t235;
t17 = -t53 * qJD(4) - t152 * t43 + t133;
t135 = t156 * qJD(2);
t237 = t152 * t61;
t52 = t135 - t237;
t165 = t16 * t156 + (-t152 * t53 - t156 * t52) * qJD(4) - t17 * t152;
t79 = t157 * t125 - t153 * t248;
t80 = t153 * t125 + t157 * t248;
t142 = qJD(4) + qJD(5);
t199 = pkin(8) * t143 + t61;
t51 = t156 * t199 + t212;
t13 = qJDD(4) * pkin(4) + t133 + (-pkin(8) * t141 - t43) * t152 - t51 * qJD(4);
t202 = t143 * t216;
t222 = t156 * t141;
t173 = t202 - t222;
t14 = -pkin(8) * t173 + t16;
t214 = qJD(5) * t151;
t50 = -t152 * t199 + t135;
t49 = qJD(4) * pkin(4) + t50;
t3 = (qJD(5) * t49 + t14) * t155 + t151 * t13 - t51 * t214;
t159 = -pkin(8) - pkin(7);
t78 = pkin(7) + t80;
t249 = -pkin(8) - t78;
t247 = g(1) * t156;
t246 = t141 * pkin(3);
t245 = t143 * pkin(3);
t244 = t156 * pkin(4);
t226 = t151 * t152;
t205 = t143 * t226;
t223 = t155 * t156;
t72 = -t143 * t223 + t205;
t243 = t74 * t72;
t109 = t159 * t152;
t138 = t156 * pkin(8);
t110 = t156 * pkin(7) + t138;
t65 = t155 * t109 - t151 * t110;
t68 = t157 * t105 - t153 * t209;
t85 = -t223 + t226;
t203 = qJD(4) * t159;
t91 = t152 * t203;
t92 = t156 * t203;
t242 = qJD(5) * t65 + t151 * t92 + t155 * t91 + t85 * t68;
t66 = t151 * t109 + t155 * t110;
t241 = -qJD(5) * t66 - t151 * t91 + t155 * t92 + t86 * t68;
t191 = t250 * t153 + t253 * t157;
t44 = t191 - t246;
t60 = -t68 - t245;
t240 = t44 * t152 + t60 * t215;
t225 = t152 * t141;
t181 = t151 * t225 - t155 * t222;
t58 = t142 * t86;
t31 = t143 * t58 + t181;
t179 = t142 * t226;
t57 = t179 + t251;
t239 = -t86 * t31 + t57 * t72;
t238 = t151 * t51;
t236 = t155 * t51;
t234 = t68 * t143;
t233 = t69 * t143;
t70 = t79 * qJD(3);
t232 = t70 * t143;
t71 = t80 * qJD(3);
t231 = t71 * t143;
t148 = qJ(4) + qJ(5);
t137 = cos(t148);
t229 = t123 * t137;
t228 = t124 * t137;
t227 = t143 * t152;
t145 = t152 ^ 2;
t146 = t156 ^ 2;
t219 = t145 - t146;
t218 = t145 + t146;
t210 = t184 * t156 + t60 * t216;
t208 = pkin(4) * t216;
t139 = t143 ^ 2;
t204 = t152 * t139 * t156;
t129 = pkin(3) + t244;
t198 = qJD(4) * t249;
t27 = pkin(4) * t173 + t44;
t54 = -t129 * t143 - t68;
t196 = g(2) * t228 + g(3) * t229 + t27 * t85 + t54 * t58;
t194 = t218 * t141;
t193 = t123 * t159 - t124 * t129;
t190 = -t141 * t224 + t251 * t143 - t151 * t222;
t188 = t156 * t202;
t187 = -t69 + t208;
t77 = -pkin(3) - t79;
t130 = sin(t144);
t154 = sin(qJ(1));
t186 = -t154 * pkin(1) - pkin(2) * t130;
t131 = cos(t144);
t158 = cos(qJ(1));
t185 = -t158 * pkin(1) - pkin(2) * t131;
t182 = g(2) * t158 + g(3) * t154;
t30 = t143 * t179 + t190;
t180 = -t85 * t30 + t74 * t58;
t18 = t155 * t49 - t238;
t19 = t151 * t49 + t236;
t4 = -qJD(5) * t19 + t155 * t13 - t151 * t14;
t178 = t18 * t57 - t19 * t58 - t3 * t85 - t4 * t86 + t220;
t140 = qJDD(4) + qJDD(5);
t45 = t86 * t140 - t57 * t142;
t62 = t249 * t152;
t63 = t156 * t78 + t138;
t33 = -t151 * t63 + t155 * t62;
t34 = t151 * t62 + t155 * t63;
t177 = t52 * t152 - t53 * t156;
t176 = -t123 * t129 - t124 * t159;
t175 = -t192 - t220;
t174 = -t191 + t184;
t160 = qJD(4) ^ 2;
t172 = -pkin(7) * t160 + t233 + t246;
t171 = -t143 * t60 - t220;
t170 = -t141 * t77 - t160 * t78 - t231;
t169 = -pkin(7) * qJDD(4) + (t68 - t245) * qJD(4);
t168 = -qJDD(4) * t78 + (t143 * t77 - t70) * qJD(4);
t136 = sin(t148);
t167 = -t136 * t184 + t27 * t86 - t54 * t57;
t164 = t220 + t165;
t163 = g(1) * t136 - g(2) * t229 + g(3) * t228 + t54 * t72 - t3;
t162 = -g(1) * t137 - t220 * t136 - t54 * t74 + t4;
t147 = qJDD(2) - g(1);
t114 = t124 * pkin(7);
t96 = qJDD(4) * t156 - t160 * t152;
t95 = qJDD(4) * t152 + t160 * t156;
t76 = t146 * t141 - 0.2e1 * t188;
t75 = t145 * t141 + 0.2e1 * t188;
t67 = t77 - t244;
t64 = t71 + t208;
t59 = -0.2e1 * t219 * t143 * qJD(4) + 0.2e1 * t152 * t222;
t46 = -t85 * t140 - t58 * t142;
t42 = -t152 * t70 + t156 * t198;
t41 = t152 * t198 + t156 * t70;
t32 = -t72 ^ 2 + t74 ^ 2;
t22 = -t190 + (-t205 + t72) * t142;
t21 = t155 * t50 - t238;
t20 = -t151 * t50 - t236;
t11 = t31 * t85 + t72 * t58;
t10 = -t30 * t86 - t74 * t57;
t7 = -qJD(5) * t34 - t151 * t41 + t155 * t42;
t6 = qJD(5) * t33 + t151 * t42 + t155 * t41;
t5 = -t180 + t239;
t1 = [0, 0, 0, 0, 0, qJDD(1), t182, -g(2) * t154 + g(3) * t158, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(2) * t131 + g(3) * t130 + 0.2e1 * t150 * t230, -g(2) * t130 + g(3) * t131 - 0.2e1 * t206, 0, (t182 + (t149 ^ 2 + t150 ^ 2) * t230) * pkin(1), 0, 0, 0, 0, 0, t141, t79 * t141 + t174 - t231, -t80 * t141 + t175 - t232, 0, -g(2) * t185 - g(3) * t186 - t191 * t79 + t192 * t80 - t68 * t71 + t69 * t70, t75, t59, t95, t76, t96, 0, t168 * t152 + (t170 - t44) * t156 + t210, t168 * t156 + (-t170 - t184) * t152 + t240, t78 * t194 + t218 * t232 + t164, t44 * t77 + t60 * t71 - g(2) * (-t124 * pkin(3) - t123 * pkin(7) + t185) - g(3) * (-t123 * pkin(3) + t114 + t186) - t177 * t70 + t165 * t78, t10, t5, t45, t11, t46, 0, t33 * t140 + t7 * t142 + t67 * t31 + t64 * t72 + t196, -t34 * t140 - t6 * t142 - t67 * t30 + t64 * t74 + t167, t33 * t30 - t34 * t31 - t6 * t72 - t7 * t74 + t178, t3 * t34 + t19 * t6 + t4 * t33 + t18 * t7 + t27 * t67 + t54 * t64 - g(2) * (t185 + t193) - g(3) * (t176 + t186); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, t96, -t95, 0, -qJD(4) * t177 + t16 * t152 + t17 * t156 - g(1), 0, 0, 0, 0, 0, 0, t46, -t45, t180 + t239, -t18 * t58 - t19 * t57 + t3 * t86 - t4 * t85 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t174 + t233, t175 + t234, 0, 0, t75, t59, t95, t76, t96, 0, t169 * t152 + (t172 - t44) * t156 + t210, t169 * t156 + (-t172 - t184) * t152 + t240, pkin(7) * t194 - t218 * t234 + t164, -g(3) * t114 - t60 * t69 + t177 * t68 + (t184 - t44) * pkin(3) + (t165 + t117) * pkin(7), t10, t5, t45, t11, t46, 0, -t129 * t31 + t65 * t140 + t241 * t142 + t187 * t72 + t196, t129 * t30 - t66 * t140 - t242 * t142 + t187 * t74 + t167, -t241 * t74 - t242 * t72 + t65 * t30 - t66 * t31 + t178, -g(2) * t193 - g(3) * t176 - t27 * t129 + t241 * t18 + t187 * t54 + t242 * t19 + t3 * t66 + t4 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t219 * t139, t225, t204, t222, qJDD(4), -t247 + t133 + (t53 - t235) * qJD(4) + (t171 - t252) * t152, g(1) * t152 + (t52 + t237) * qJD(4) + t171 * t156 - t207, 0, 0, t243, t32, t22, -t243, -t181, t140, -t20 * t142 + (t140 * t155 - t142 * t214 - t72 * t227) * pkin(4) + t162, t21 * t142 + (-t140 * t151 - t142 * t213 - t74 * t227) * pkin(4) + t163, (t19 + t20) * t74 + (-t18 + t21) * t72 + (-t151 * t31 + t155 * t30 + (t151 * t74 - t155 * t72) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t247 + t151 * t3 + t155 * t4 + (-t151 * t18 + t155 * t19) * qJD(5) + (-t143 * t54 - t220) * t152) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t32, t22, -t243, -t181, t140, t19 * t142 + t162, t18 * t142 + t163, 0, 0;];
tau_reg = t1;
