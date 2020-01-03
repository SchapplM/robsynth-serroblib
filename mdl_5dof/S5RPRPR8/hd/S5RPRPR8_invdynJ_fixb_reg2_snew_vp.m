% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:12
% EndTime: 2019-12-31 18:22:19
% DurationCPUTime: 2.31s
% Computational Cost: add. (8094->301), mult. (16903->431), div. (0->0), fcn. (11070->10), ass. (0->180)
t172 = qJD(1) ^ 2;
t161 = sin(pkin(9));
t163 = cos(pkin(9));
t167 = sin(qJ(3));
t195 = qJD(1) * t167;
t127 = -qJD(3) * t163 + t161 * t195;
t129 = qJD(3) * t161 + t163 * t195;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t103 = t127 * t169 + t129 * t166;
t170 = cos(qJ(3));
t192 = qJD(1) * qJD(3);
t186 = t170 * t192;
t191 = t167 * qJDD(1);
t135 = t186 + t191;
t111 = qJDD(3) * t161 + t135 * t163;
t181 = -t163 * qJDD(3) + t135 * t161;
t67 = -qJD(5) * t103 + t111 * t169 - t166 * t181;
t194 = t170 * qJD(1);
t146 = -qJD(5) + t194;
t91 = t103 * t146;
t220 = t91 + t67;
t150 = t167 * t192;
t190 = t170 * qJDD(1);
t136 = -t150 + t190;
t202 = t129 * t127;
t175 = -t136 - t202;
t219 = t161 * t175;
t218 = t163 * t175;
t130 = -qJDD(5) + t136;
t105 = -t127 * t166 + t129 * t169;
t203 = t105 * t103;
t174 = -t130 - t203;
t217 = t166 * t174;
t216 = t169 * t174;
t116 = t127 * t194;
t95 = -t111 + t116;
t117 = t129 * t194;
t93 = -t181 - t117;
t183 = t166 * t111 + t169 * t181;
t53 = (qJD(5) + t146) * t105 + t183;
t101 = t103 ^ 2;
t102 = t105 ^ 2;
t215 = t127 ^ 2;
t126 = t129 ^ 2;
t144 = t146 ^ 2;
t214 = qJD(3) ^ 2;
t177 = t135 + t186;
t168 = sin(qJ(1));
t171 = cos(qJ(1));
t185 = t168 * g(1) - g(2) * t171;
t131 = qJDD(1) * pkin(1) + t185;
t179 = g(1) * t171 + g(2) * t168;
t133 = -pkin(1) * t172 - t179;
t162 = sin(pkin(8));
t164 = cos(pkin(8));
t182 = t164 * t131 - t162 * t133;
t99 = -qJDD(1) * pkin(2) - t172 * pkin(6) - t182;
t75 = -t177 * qJ(4) + (-t136 + t150) * pkin(3) + t99;
t196 = t162 * t131 + t164 * t133;
t100 = -pkin(2) * t172 + qJDD(1) * pkin(6) + t196;
t178 = -pkin(3) * t170 - qJ(4) * t167;
t180 = t172 * t178 + t100;
t197 = -g(3) + qJDD(2);
t184 = t167 * t197;
t79 = -pkin(3) * t214 + qJDD(3) * qJ(4) + t170 * t180 + t184;
t36 = 0.2e1 * qJD(4) * t129 + t161 * t79 - t163 * t75;
t33 = t175 * pkin(4) + pkin(7) * t95 - t36;
t112 = -pkin(4) * t194 - pkin(7) * t129;
t37 = -0.2e1 * qJD(4) * t127 + t161 * t75 + t163 * t79;
t34 = -pkin(4) * t215 - pkin(7) * t181 + t112 * t194 + t37;
t13 = t166 * t34 - t169 * t33;
t14 = t166 * t33 + t169 * t34;
t7 = -t13 * t169 + t14 * t166;
t213 = t161 * t7;
t212 = t163 * t7;
t149 = t170 * t197;
t78 = -qJDD(3) * pkin(3) - t214 * qJ(4) + t167 * t180 + qJDD(4) - t149;
t211 = t161 * t78;
t96 = t136 - t202;
t210 = t161 * t96;
t209 = t163 * t78;
t208 = t163 * t96;
t40 = pkin(4) * t181 - pkin(7) * t215 + t129 * t112 + t78;
t207 = t166 * t40;
t70 = t130 - t203;
t206 = t166 * t70;
t205 = t169 * t40;
t204 = t169 * t70;
t201 = t146 * t166;
t200 = t146 * t169;
t145 = t170 * t172 * t167;
t140 = qJDD(3) + t145;
t199 = t167 * t140;
t141 = qJDD(3) - t145;
t198 = t170 * t141;
t189 = t170 * t203;
t188 = t170 * t202;
t187 = pkin(1) * t162 + pkin(6);
t8 = t13 * t166 + t14 * t169;
t21 = t161 * t36 + t163 * t37;
t85 = t100 * t167 - t149;
t86 = t170 * t100 + t184;
t60 = t167 * t85 + t170 * t86;
t176 = t161 * t37 - t163 * t36;
t173 = -pkin(1) * t164 - pkin(2) + t178;
t158 = t170 ^ 2;
t157 = t167 ^ 2;
t154 = t158 * t172;
t153 = t157 * t172;
t143 = -t154 - t214;
t142 = -t153 - t214;
t139 = t153 + t154;
t138 = (t157 + t158) * qJDD(1);
t137 = -0.2e1 * t150 + t190;
t134 = 0.2e1 * t186 + t191;
t124 = t170 * t136;
t115 = -t126 - t154;
t114 = -t126 + t154;
t113 = -t154 + t215;
t109 = -t142 * t167 - t198;
t108 = t143 * t170 - t199;
t106 = -t154 - t215;
t94 = t111 + t116;
t92 = -t117 + t181;
t89 = -t126 - t215;
t88 = -t102 + t144;
t87 = t101 - t144;
t84 = -t102 - t144;
t82 = -t115 * t161 + t208;
t81 = t115 * t163 + t210;
t80 = t102 - t101;
t76 = -t144 - t101;
t74 = t106 * t163 - t219;
t73 = t106 * t161 + t218;
t66 = -qJD(5) * t105 - t183;
t65 = -t161 * t95 + t163 * t93;
t63 = (t103 * t169 - t105 * t166) * t146;
t62 = (t103 * t166 + t105 * t169) * t146;
t61 = -t101 - t102;
t59 = t167 * t94 + t170 * t82;
t58 = t167 * t92 + t170 * t74;
t57 = -t91 + t67;
t52 = (qJD(5) - t146) * t105 + t183;
t51 = t169 * t87 + t206;
t50 = -t166 * t88 + t216;
t49 = t166 * t87 - t204;
t48 = t169 * t88 + t217;
t47 = t105 * t201 + t169 * t67;
t46 = -t105 * t200 + t166 * t67;
t45 = -t103 * t200 - t166 * t66;
t44 = -t103 * t201 + t169 * t66;
t42 = -t166 * t84 + t204;
t41 = t169 * t84 + t206;
t39 = t169 * t76 - t217;
t38 = t166 * t76 + t216;
t31 = t166 * t57 - t169 * t53;
t30 = -t166 * t220 - t169 * t52;
t29 = -t166 * t53 - t169 * t57;
t28 = -t166 * t52 + t169 * t220;
t27 = -t161 * t41 + t163 * t42;
t26 = t161 * t42 + t163 * t41;
t25 = -pkin(7) * t41 + t205;
t24 = -pkin(7) * t38 + t207;
t23 = -t161 * t38 + t163 * t39;
t22 = t161 * t39 + t163 * t38;
t19 = -pkin(4) * t220 + pkin(7) * t42 + t207;
t18 = t167 * t220 + t170 * t27;
t16 = -pkin(4) * t52 + pkin(7) * t39 - t205;
t15 = t167 * t52 + t170 * t23;
t11 = -t161 * t29 + t163 * t31;
t10 = t161 * t31 + t163 * t29;
t9 = t11 * t170 + t167 * t61;
t6 = -pkin(4) * t40 + pkin(7) * t8;
t5 = -pkin(7) * t29 - t7;
t4 = -pkin(4) * t61 + pkin(7) * t31 + t8;
t3 = t163 * t8 - t213;
t2 = t161 * t8 + t212;
t1 = t167 * t40 + t170 * t3;
t12 = [0, 0, 0, 0, 0, qJDD(1), t185, t179, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t164 - t162 * t172) + t182, pkin(1) * (-qJDD(1) * t162 - t164 * t172) - t196, 0, pkin(1) * (t162 * t196 + t164 * t182), t177 * t167, t134 * t170 + t137 * t167, t199 + t170 * (-t153 + t214), -t167 * t186 + t124, t167 * (t154 - t214) + t198, 0, -t170 * t99 + pkin(2) * t137 + pkin(6) * t108 + pkin(1) * (t108 * t162 + t137 * t164), t167 * t99 - pkin(2) * t134 + pkin(6) * t109 + pkin(1) * (t109 * t162 - t134 * t164), pkin(2) * t139 + pkin(6) * t138 + pkin(1) * (t138 * t162 + t139 * t164) + t60, -pkin(2) * t99 + pkin(6) * t60 + pkin(1) * (t162 * t60 - t164 * t99), t167 * (t111 * t163 + t117 * t161) - t188, t167 * (-t161 * t94 - t163 * t92) + t170 * (-t126 + t215), t167 * (-t114 * t161 + t218) + t170 * t95, t167 * (-t116 * t163 + t161 * t181) + t188, t167 * (t113 * t163 + t210) - t170 * t93, t124 + t167 * (t127 * t163 - t129 * t161) * t194, t167 * (-qJ(4) * t73 + t211) + t170 * (-pkin(3) * t73 + t36) - pkin(2) * t73 + pkin(6) * t58 + pkin(1) * (t162 * t58 - t164 * t73), t167 * (-qJ(4) * t81 + t209) + t170 * (-pkin(3) * t81 + t37) - pkin(2) * t81 + pkin(6) * t59 + pkin(1) * (t162 * t59 - t164 * t81), -t167 * t176 + t187 * (t167 * t89 + t170 * t65) + t173 * (t161 * t93 + t163 * t95), t187 * (t167 * t78 + t170 * t21) + t173 * t176, t167 * (-t161 * t46 + t163 * t47) - t189, t167 * (-t161 * t28 + t163 * t30) - t170 * t80, t167 * (-t161 * t48 + t163 * t50) - t170 * t57, t167 * (-t161 * t44 + t163 * t45) + t189, t167 * (-t161 * t49 + t163 * t51) + t170 * t53, t167 * (-t161 * t62 + t163 * t63) + t170 * t130, t167 * (-qJ(4) * t22 - t16 * t161 + t163 * t24) + t170 * (-pkin(3) * t22 - pkin(4) * t38 + t13) - pkin(2) * t22 + pkin(6) * t15 + pkin(1) * (t15 * t162 - t164 * t22), t167 * (-qJ(4) * t26 - t161 * t19 + t163 * t25) + t170 * (-pkin(3) * t26 - pkin(4) * t41 + t14) - pkin(2) * t26 + pkin(6) * t18 + pkin(1) * (t162 * t18 - t164 * t26), t167 * (-qJ(4) * t10 - t161 * t4 + t163 * t5) + t170 * (-pkin(3) * t10 - pkin(4) * t29) - pkin(2) * t10 + pkin(6) * t9 + pkin(1) * (-t10 * t164 + t162 * t9), t167 * (-pkin(7) * t212 - qJ(4) * t2 - t161 * t6) + t170 * (-pkin(3) * t2 - pkin(4) * t7) - pkin(2) * t2 + pkin(6) * t1 + pkin(1) * (t1 * t162 - t164 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, 0, 0, 0, 0, 0, 0, t140 * t170 + t143 * t167, -t141 * t167 + t142 * t170, 0, t167 * t86 - t170 * t85, 0, 0, 0, 0, 0, 0, t167 * t74 - t170 * t92, t167 * t82 - t170 * t94, t167 * t65 - t170 * t89, t167 * t21 - t170 * t78, 0, 0, 0, 0, 0, 0, t167 * t23 - t170 * t52, t167 * t27 - t170 * t220, t11 * t167 - t170 * t61, t167 * t3 - t170 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t153 - t154, t191, t145, t190, qJDD(3), -t85, -t86, 0, 0, t111 * t161 - t117 * t163, -t161 * t92 + t163 * t94, t114 * t163 + t219, -t116 * t161 - t163 * t181, t113 * t161 - t208, (t127 * t161 + t129 * t163) * t194, -pkin(3) * t92 + qJ(4) * t74 - t209, -pkin(3) * t94 + qJ(4) * t82 + t211, -pkin(3) * t89 + qJ(4) * t65 + t21, -pkin(3) * t78 + qJ(4) * t21, t161 * t47 + t163 * t46, t161 * t30 + t163 * t28, t161 * t50 + t163 * t48, t161 * t45 + t163 * t44, t161 * t51 + t163 * t49, t161 * t63 + t163 * t62, -pkin(3) * t52 + qJ(4) * t23 + t16 * t163 + t161 * t24, -pkin(3) * t220 + qJ(4) * t27 + t161 * t25 + t163 * t19, -pkin(3) * t61 + qJ(4) * t11 + t161 * t5 + t163 * t4, -pkin(3) * t40 - pkin(7) * t213 + qJ(4) * t3 + t163 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t94, t89, t78, 0, 0, 0, 0, 0, 0, t52, t220, t61, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, t80, t57, -t203, -t53, -t130, -t13, -t14, 0, 0;];
tauJ_reg = t12;
