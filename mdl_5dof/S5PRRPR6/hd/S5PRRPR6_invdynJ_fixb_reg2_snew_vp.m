% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:50
% EndTime: 2019-12-05 16:33:01
% DurationCPUTime: 2.60s
% Computational Cost: add. (8421->304), mult. (17586->451), div. (0->0), fcn. (12674->12), ass. (0->186)
t161 = sin(pkin(5));
t163 = cos(pkin(5));
t203 = sin(pkin(9));
t204 = cos(pkin(9));
t176 = t203 * g(1) - t204 * g(2);
t174 = t163 * t176;
t196 = -g(3) + qJDD(1);
t222 = t161 * t196 + t174;
t171 = qJD(2) ^ 2;
t160 = sin(pkin(10));
t162 = cos(pkin(10));
t166 = sin(qJ(3));
t195 = qJD(2) * t166;
t128 = -t162 * qJD(3) + t160 * t195;
t130 = t160 * qJD(3) + t162 * t195;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t105 = t168 * t128 + t165 * t130;
t169 = cos(qJ(3));
t192 = qJD(2) * qJD(3);
t186 = t169 * t192;
t191 = t166 * qJDD(2);
t134 = t186 + t191;
t114 = t160 * qJDD(3) + t162 * t134;
t182 = -t162 * qJDD(3) + t160 * t134;
t73 = -t105 * qJD(5) + t168 * t114 - t165 * t182;
t194 = t169 * qJD(2);
t147 = -qJD(5) + t194;
t93 = t105 * t147;
t221 = t93 + t73;
t150 = t166 * t192;
t190 = t169 * qJDD(2);
t135 = -t150 + t190;
t201 = t130 * t128;
t177 = -t135 - t201;
t220 = t160 * t177;
t219 = t162 * t177;
t131 = -qJDD(5) + t135;
t107 = -t165 * t128 + t168 * t130;
t202 = t107 * t105;
t175 = -t131 - t202;
t218 = t165 * t175;
t217 = t168 * t175;
t119 = t128 * t194;
t97 = -t114 + t119;
t120 = t130 * t194;
t95 = -t182 - t120;
t183 = t165 * t114 + t168 * t182;
t54 = (qJD(5) + t147) * t107 + t183;
t103 = t105 ^ 2;
t104 = t107 ^ 2;
t216 = t128 ^ 2;
t127 = t130 ^ 2;
t145 = t147 ^ 2;
t215 = qJD(3) ^ 2;
t173 = -t161 * t176 + t163 * t196;
t172 = t166 * t173;
t180 = -t169 * pkin(3) - t166 * qJ(4);
t139 = -t204 * g(1) - t203 * g(2);
t167 = sin(qJ(2));
t170 = cos(qJ(2));
t102 = t170 * t139 + t222 * t167;
t88 = -t171 * pkin(2) + qJDD(2) * pkin(7) + t102;
t184 = t171 * t180 + t88;
t64 = -t215 * pkin(3) + qJDD(3) * qJ(4) + t184 * t169 + t172;
t179 = t134 + t186;
t181 = t167 * t139 - t222 * t170;
t87 = -qJDD(2) * pkin(2) - t171 * pkin(7) + t181;
t71 = -t179 * qJ(4) + (-t135 + t150) * pkin(3) + t87;
t36 = 0.2e1 * qJD(4) * t130 + t160 * t64 - t162 * t71;
t27 = t177 * pkin(4) + t97 * pkin(8) - t36;
t115 = -pkin(4) * t194 - t130 * pkin(8);
t37 = -0.2e1 * qJD(4) * t128 + t160 * t71 + t162 * t64;
t30 = -t216 * pkin(4) - t182 * pkin(8) + t115 * t194 + t37;
t11 = t165 * t30 - t168 * t27;
t12 = t165 * t27 + t168 * t30;
t7 = -t168 * t11 + t165 * t12;
t214 = t160 * t7;
t213 = t162 * t7;
t112 = t169 * t173;
t63 = -qJDD(3) * pkin(3) - t215 * qJ(4) + t184 * t166 + qJDD(4) - t112;
t212 = t160 * t63;
t99 = t135 - t201;
t211 = t160 * t99;
t210 = t162 * t63;
t209 = t162 * t99;
t38 = t182 * pkin(4) - t216 * pkin(8) + t130 * t115 + t63;
t208 = t165 * t38;
t74 = t131 - t202;
t207 = t165 * t74;
t206 = t168 * t38;
t205 = t168 * t74;
t200 = t147 * t165;
t199 = t147 * t168;
t146 = t166 * t171 * t169;
t140 = qJDD(3) + t146;
t198 = t166 * t140;
t141 = qJDD(3) - t146;
t197 = t169 * t141;
t189 = t169 * t202;
t188 = t169 * t201;
t187 = t162 * t194;
t8 = t165 * t11 + t168 * t12;
t21 = t160 * t36 + t162 * t37;
t81 = t166 * t88 - t112;
t82 = t169 * t88 + t172;
t41 = t166 * t81 + t169 * t82;
t20 = t160 * t37 - t162 * t36;
t178 = -pkin(2) + t180;
t157 = t169 ^ 2;
t156 = t166 ^ 2;
t154 = t157 * t171;
t153 = t156 * t171;
t144 = -t154 - t215;
t143 = -t153 - t215;
t138 = t153 + t154;
t137 = (t156 + t157) * qJDD(2);
t136 = -0.2e1 * t150 + t190;
t133 = 0.2e1 * t186 + t191;
t125 = t169 * t135;
t118 = -t127 - t154;
t117 = -t127 + t154;
t116 = -t154 + t216;
t111 = -t166 * t143 - t197;
t110 = t169 * t144 - t198;
t108 = -t154 - t216;
t96 = t114 + t119;
t94 = -t120 + t182;
t91 = -t127 - t216;
t90 = -t104 + t145;
t89 = t103 - t145;
t86 = -t104 - t145;
t85 = -t160 * t118 + t209;
t84 = t162 * t118 + t211;
t83 = t104 - t103;
t79 = -t145 - t103;
t78 = t162 * t108 - t220;
t77 = t160 * t108 + t219;
t72 = -t107 * qJD(5) - t183;
t70 = -t160 * t97 + t162 * t95;
t69 = t160 * t95 + t162 * t97;
t68 = (t105 * t168 - t107 * t165) * t147;
t67 = (t105 * t165 + t107 * t168) * t147;
t61 = -t103 - t104;
t60 = t166 * t96 + t169 * t85;
t59 = t166 * t94 + t169 * t78;
t58 = -t93 + t73;
t53 = (qJD(5) - t147) * t107 + t183;
t52 = t168 * t89 + t207;
t51 = -t165 * t90 + t217;
t50 = t165 * t89 - t205;
t49 = t168 * t90 + t218;
t48 = t107 * t200 + t168 * t73;
t47 = -t107 * t199 + t165 * t73;
t46 = -t105 * t199 - t165 * t72;
t45 = -t105 * t200 + t168 * t72;
t44 = t166 * t91 + t169 * t70;
t43 = -t165 * t86 + t205;
t42 = t168 * t86 + t207;
t40 = t168 * t79 - t218;
t39 = t165 * t79 + t217;
t34 = t165 * t58 - t168 * t54;
t33 = -t165 * t221 - t168 * t53;
t32 = -t165 * t54 - t168 * t58;
t31 = -t165 * t53 + t168 * t221;
t29 = -t160 * t42 + t162 * t43;
t28 = t160 * t43 + t162 * t42;
t25 = -pkin(8) * t42 + t206;
t24 = -t160 * t39 + t162 * t40;
t23 = t160 * t40 + t162 * t39;
t22 = -pkin(8) * t39 + t208;
t19 = t166 * t221 + t169 * t29;
t18 = -pkin(4) * t221 + pkin(8) * t43 + t208;
t17 = -pkin(4) * t53 + pkin(8) * t40 - t206;
t16 = t166 * t53 + t169 * t24;
t15 = t166 * t63 + t169 * t21;
t14 = -t160 * t32 + t162 * t34;
t13 = t160 * t34 + t162 * t32;
t9 = t169 * t14 + t166 * t61;
t6 = -pkin(4) * t38 + pkin(8) * t8;
t5 = -pkin(8) * t32 - t7;
t4 = -pkin(4) * t61 + pkin(8) * t34 + t8;
t3 = t162 * t8 - t214;
t2 = t160 * t8 + t213;
t1 = t166 * t38 + t169 * t3;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, 0, 0, 0, 0, 0, (qJDD(2) * t170 - t167 * t171) * t161, (-qJDD(2) * t167 - t170 * t171) * t161, 0, t163 ^ 2 * t196 + (t167 * t102 - t170 * t181 - t174) * t161, 0, 0, 0, 0, 0, 0, t163 * (t169 * t140 + t166 * t144) + (t167 * t110 + t170 * t136) * t161, t163 * (-t166 * t141 + t169 * t143) + (t167 * t111 - t170 * t133) * t161, (t137 * t167 + t138 * t170) * t161, t163 * (t166 * t82 - t169 * t81) + (t167 * t41 - t170 * t87) * t161, 0, 0, 0, 0, 0, 0, t163 * (t166 * t78 - t169 * t94) + (t167 * t59 - t170 * t77) * t161, t163 * (t166 * t85 - t169 * t96) + (t167 * t60 - t170 * t84) * t161, t163 * (t166 * t70 - t169 * t91) + (t167 * t44 - t170 * t69) * t161, t163 * (t166 * t21 - t169 * t63) + (t167 * t15 - t170 * t20) * t161, 0, 0, 0, 0, 0, 0, t163 * (t166 * t24 - t169 * t53) + (t167 * t16 - t170 * t23) * t161, t163 * (t166 * t29 - t169 * t221) + (t167 * t19 - t170 * t28) * t161, t163 * (t166 * t14 - t169 * t61) + (-t170 * t13 + t167 * t9) * t161, t163 * (t166 * t3 - t169 * t38) + (t167 * t1 - t170 * t2) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t181, -t102, 0, 0, t179 * t166, t169 * t133 + t166 * t136, t198 + t169 * (-t153 + t215), -t166 * t186 + t125, t166 * (t154 - t215) + t197, 0, pkin(2) * t136 + pkin(7) * t110 - t169 * t87, -pkin(2) * t133 + pkin(7) * t111 + t166 * t87, pkin(2) * t138 + pkin(7) * t137 + t41, -pkin(2) * t87 + pkin(7) * t41, t166 * (t162 * t114 + t160 * t120) - t188, t166 * (-t160 * t96 - t162 * t94) + t169 * (-t127 + t216), t166 * (-t160 * t117 + t219) + t169 * t97, t166 * (-t128 * t187 + t160 * t182) + t188, t166 * (t162 * t116 + t211) - t169 * t95, t125 + t166 * (t128 * t162 - t130 * t160) * t194, t166 * (-qJ(4) * t77 + t212) + t169 * (-pkin(3) * t77 + t36) - pkin(2) * t77 + pkin(7) * t59, t166 * (-qJ(4) * t84 + t210) + t169 * (-pkin(3) * t84 + t37) - pkin(2) * t84 + pkin(7) * t60, pkin(7) * t44 - t166 * t20 + t178 * t69, pkin(7) * t15 + t178 * t20, t166 * (-t160 * t47 + t162 * t48) - t189, t166 * (-t160 * t31 + t162 * t33) - t169 * t83, t166 * (-t160 * t49 + t162 * t51) - t169 * t58, t166 * (-t160 * t45 + t162 * t46) + t189, t166 * (-t160 * t50 + t162 * t52) + t169 * t54, t166 * (-t160 * t67 + t162 * t68) + t169 * t131, t166 * (-qJ(4) * t23 - t160 * t17 + t162 * t22) + t169 * (-pkin(3) * t23 - pkin(4) * t39 + t11) - pkin(2) * t23 + pkin(7) * t16, t166 * (-qJ(4) * t28 - t160 * t18 + t162 * t25) + t169 * (-pkin(3) * t28 - pkin(4) * t42 + t12) - pkin(2) * t28 + pkin(7) * t19, t166 * (-qJ(4) * t13 - t160 * t4 + t162 * t5) + t169 * (-pkin(3) * t13 - pkin(4) * t32) - pkin(2) * t13 + pkin(7) * t9, t166 * (-pkin(8) * t213 - qJ(4) * t2 - t160 * t6) + t169 * (-pkin(3) * t2 - pkin(4) * t7) - pkin(2) * t2 + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t153 - t154, t191, t146, t190, qJDD(3), -t81, -t82, 0, 0, t160 * t114 - t130 * t187, -t160 * t94 + t162 * t96, t162 * t117 + t220, -t160 * t119 - t162 * t182, t160 * t116 - t209, (t128 * t160 + t130 * t162) * t194, -pkin(3) * t94 + qJ(4) * t78 - t210, -pkin(3) * t96 + qJ(4) * t85 + t212, -pkin(3) * t91 + qJ(4) * t70 + t21, -pkin(3) * t63 + qJ(4) * t21, t160 * t48 + t162 * t47, t160 * t33 + t162 * t31, t160 * t51 + t162 * t49, t160 * t46 + t162 * t45, t160 * t52 + t162 * t50, t160 * t68 + t162 * t67, -pkin(3) * t53 + qJ(4) * t24 + t160 * t22 + t162 * t17, -pkin(3) * t221 + qJ(4) * t29 + t160 * t25 + t162 * t18, -pkin(3) * t61 + qJ(4) * t14 + t160 * t5 + t162 * t4, -pkin(3) * t38 - pkin(8) * t214 + qJ(4) * t3 + t162 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t96, t91, t63, 0, 0, 0, 0, 0, 0, t53, t221, t61, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t83, t58, -t202, -t54, -t131, -t11, -t12, 0, 0;];
tauJ_reg = t10;
