% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:18
% EndTime: 2020-01-03 11:28:34
% DurationCPUTime: 2.50s
% Computational Cost: add. (7501->276), mult. (17010->410), div. (0->0), fcn. (12044->10), ass. (0->162)
t195 = sin(qJ(1));
t196 = cos(qJ(1));
t164 = t195 * g(2) - t196 * g(3);
t197 = qJD(1) ^ 2;
t128 = -t197 * pkin(1) - t164;
t148 = sin(pkin(8));
t150 = cos(pkin(8));
t165 = -t196 * g(2) - t195 * g(3);
t161 = qJDD(1) * pkin(1) + t165;
t179 = t150 * t128 + t148 * t161;
t208 = -t197 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t179;
t152 = sin(qJ(4));
t147 = sin(pkin(9));
t149 = cos(pkin(9));
t154 = cos(qJ(4));
t120 = (-t147 * t152 + t149 * t154) * qJD(1);
t166 = t147 * t154 + t149 * t152;
t122 = t166 * qJD(1);
t182 = t122 * t120;
t205 = qJDD(4) + t182;
t210 = t152 * t205;
t209 = t154 * t205;
t156 = t147 ^ 2;
t158 = t149 ^ 2;
t130 = (t156 + t158) * t197;
t151 = sin(qJ(5));
t142 = qJDD(4) + qJDD(5);
t153 = cos(qJ(5));
t94 = -t153 * t120 + t122 * t151;
t96 = t120 * t151 + t122 * t153;
t70 = t96 * t94;
t202 = -t70 + t142;
t207 = t151 * t202;
t206 = t153 * t202;
t168 = -t148 * t128 + t150 * t161;
t91 = -qJDD(1) * pkin(2) - t197 * qJ(3) + qJDD(3) - t168;
t204 = (pkin(1) * t148 + qJ(3)) * t130 + t91 - (pkin(1) * t150 + pkin(2)) * qJDD(1);
t175 = t122 * qJD(4);
t172 = t149 * qJDD(1);
t173 = t147 * qJDD(1);
t83 = -t152 * t173 + t154 * t172;
t105 = t83 - t175;
t119 = t166 * qJDD(1);
t176 = t120 * qJD(4);
t107 = t119 + t176;
t62 = -qJD(5) * t94 + t105 * t151 + t107 * t153;
t143 = qJD(4) + qJD(5);
t89 = t143 * t94;
t203 = t62 - t89;
t145 = -g(1) + qJDD(2);
t135 = t149 * t145;
t200 = t135 + (pkin(3) * t197 * t149 - pkin(6) * qJDD(1) - t208) * t147;
t177 = t197 * t158;
t82 = t147 * t145 + t208 * t149;
t76 = -pkin(3) * t177 + pkin(6) * t172 + t82;
t50 = t152 * t76 - t154 * t200;
t198 = -t50 + (-t107 + t176) * pkin(7);
t92 = t94 ^ 2;
t93 = t96 ^ 2;
t117 = t120 ^ 2;
t118 = t122 ^ 2;
t141 = t143 ^ 2;
t160 = pkin(4) * t205 + t198;
t109 = qJD(4) * pkin(4) - pkin(7) * t122;
t51 = t200 * t152 + t154 * t76;
t33 = -t117 * pkin(4) + t105 * pkin(7) - qJD(4) * t109 + t51;
t17 = t151 * t33 - t153 * t160;
t188 = t153 * t33;
t18 = t151 * t160 + t188;
t8 = t151 * t18 - t153 * t17;
t194 = t152 * t8;
t193 = t154 * t8;
t27 = t152 * t51 - t154 * t50;
t192 = t147 * t27;
t80 = -pkin(3) * t172 + t91 + (-t156 * t197 - t177) * pkin(6);
t49 = -t105 * pkin(4) - t117 * pkin(7) + t109 * t122 + t80;
t191 = t151 * t49;
t67 = t70 + t142;
t190 = t151 * t67;
t189 = t152 * t80;
t187 = t153 * t49;
t186 = t153 * t67;
t185 = t154 * t80;
t102 = qJDD(4) - t182;
t184 = t102 * t152;
t183 = t102 * t154;
t181 = t143 * t151;
t180 = t143 * t153;
t81 = t208 * t147 - t135;
t59 = t147 * t81 + t149 * t82;
t9 = t151 * t17 + t153 * t18;
t28 = t152 * t50 + t154 * t51;
t169 = -t153 * t105 + t107 * t151;
t163 = (-qJD(5) + t143) * t96 - t169;
t155 = qJD(4) ^ 2;
t140 = t158 * qJDD(1);
t139 = t156 * qJDD(1);
t129 = t140 + t139;
t112 = -t118 - t155;
t111 = -t118 + t155;
t110 = t117 - t155;
t106 = t119 + 0.2e1 * t176;
t104 = -t83 + 0.2e1 * t175;
t100 = -t155 - t117;
t87 = -t93 + t141;
t86 = t92 - t141;
t85 = -t93 - t141;
t84 = -t117 - t118;
t78 = -t112 * t152 - t183;
t77 = t112 * t154 - t184;
t74 = t119 * t152 + t154 * t83;
t73 = -t119 * t154 + t152 * t83;
t72 = t100 * t154 - t210;
t71 = t100 * t152 + t209;
t69 = t93 - t92;
t65 = -t141 - t92;
t64 = (t151 * t96 - t153 * t94) * t143;
t63 = (-t151 * t94 - t153 * t96) * t143;
t61 = -qJD(5) * t96 - t169;
t60 = -t92 - t93;
t58 = -t147 * t77 + t149 * t78;
t57 = t153 * t86 - t190;
t56 = -t151 * t87 + t206;
t55 = t151 * t86 + t186;
t54 = t153 * t87 + t207;
t53 = -t151 * t85 - t186;
t52 = t153 * t85 - t190;
t47 = -t147 * t73 + t149 * t74;
t46 = t62 + t89;
t41 = (qJD(5) + t143) * t96 + t169;
t40 = -t147 * t71 + t149 * t72;
t39 = t153 * t62 - t96 * t181;
t38 = t151 * t62 + t96 * t180;
t37 = -t151 * t61 + t94 * t180;
t36 = t153 * t61 + t94 * t181;
t35 = t153 * t65 - t207;
t34 = t151 * t65 + t206;
t31 = -t152 * t52 + t154 * t53;
t30 = t152 * t53 + t154 * t52;
t29 = -pkin(7) * t52 + t187;
t26 = t151 * t46 + t153 * t163;
t25 = -t151 * t203 - t153 * t41;
t24 = t151 * t163 - t153 * t46;
t23 = -t151 * t41 + t153 * t203;
t22 = -pkin(7) * t34 + t191;
t21 = -t152 * t34 + t154 * t35;
t20 = t152 * t35 + t154 * t34;
t19 = -pkin(4) * t203 + pkin(7) * t53 + t191;
t15 = -pkin(4) * t41 + pkin(7) * t35 - t187;
t14 = -t147 * t30 + t149 * t31;
t13 = t149 * t28 - t192;
t12 = -t152 * t24 + t154 * t26;
t11 = t152 * t26 + t154 * t24;
t10 = -t147 * t20 + t149 * t21;
t7 = -pkin(4) * t49 + pkin(7) * t9;
t6 = -pkin(7) * t24 - t8;
t5 = -t147 * t11 + t149 * t12;
t4 = -pkin(4) * t60 + pkin(7) * t26 + t9;
t3 = t154 * t9 - t194;
t2 = t152 * t9 + t193;
t1 = -t147 * t2 + t149 * t3;
t16 = [0, 0, 0, 0, 0, qJDD(1), t165, t164, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t150 - t197 * t148) + t168, pkin(1) * (-qJDD(1) * t148 - t197 * t150) - t179, 0, pkin(1) * (t148 * t179 + t150 * t168), t139, 0.2e1 * t147 * t172, 0, t140, 0, 0, -t204 * t149, t204 * t147, pkin(2) * t130 + qJ(3) * t129 + pkin(1) * (t129 * t148 + t130 * t150) + t59, -pkin(2) * t91 + qJ(3) * t59 + pkin(1) * (t148 * t59 - t150 * t91), t147 * (t107 * t154 - t152 * t175) + t149 * (t107 * t152 + t154 * t175), t147 * (-t104 * t154 - t106 * t152) + t149 * (-t104 * t152 + t106 * t154), t147 * (-t111 * t152 + t209) + t149 * (t111 * t154 + t210), t147 * (-t105 * t152 - t154 * t176) + t149 * (t105 * t154 - t152 * t176), t147 * (t110 * t154 - t184) + t149 * (t110 * t152 + t183), (t147 * (t120 * t154 + t122 * t152) + t149 * (t120 * t152 - t122 * t154)) * qJD(4), t147 * (-pkin(6) * t71 + t189) + t149 * (-pkin(3) * t104 + pkin(6) * t72 - t185) - pkin(2) * t104 + qJ(3) * t40 + pkin(1) * (-t104 * t150 + t148 * t40), t147 * (-pkin(6) * t77 + t185) + t149 * (-pkin(3) * t106 + pkin(6) * t78 + t189) - pkin(2) * t106 + qJ(3) * t58 + pkin(1) * (-t106 * t150 + t148 * t58), t147 * (-pkin(6) * t73 - t27) + t149 * (-pkin(3) * t84 + pkin(6) * t74 + t28) - pkin(2) * t84 + qJ(3) * t47 + pkin(1) * (t148 * t47 - t150 * t84), -pkin(6) * t192 + t149 * (-pkin(3) * t80 + pkin(6) * t28) - pkin(2) * t80 + qJ(3) * t13 + pkin(1) * (t13 * t148 - t150 * t80), t147 * (-t152 * t38 + t154 * t39) + t149 * (t152 * t39 + t154 * t38), t147 * (-t152 * t23 + t154 * t25) + t149 * (t152 * t25 + t154 * t23), t147 * (-t152 * t54 + t154 * t56) + t149 * (t152 * t56 + t154 * t54), t147 * (-t152 * t36 + t154 * t37) + t149 * (t152 * t37 + t154 * t36), t147 * (-t152 * t55 + t154 * t57) + t149 * (t152 * t57 + t154 * t55), t147 * (-t152 * t63 + t154 * t64) + t149 * (t152 * t64 + t154 * t63), t147 * (-pkin(6) * t20 - t15 * t152 + t154 * t22) + t149 * (-pkin(3) * t41 + pkin(6) * t21 + t15 * t154 + t152 * t22) - pkin(2) * t41 + qJ(3) * t10 + pkin(1) * (t10 * t148 - t150 * t41), t147 * (-pkin(6) * t30 - t152 * t19 + t154 * t29) + t149 * (-pkin(3) * t203 + pkin(6) * t31 + t152 * t29 + t154 * t19) - pkin(2) * t203 + qJ(3) * t14 + pkin(1) * (t14 * t148 - t150 * t203), t147 * (-pkin(6) * t11 - t152 * t4 + t154 * t6) + t149 * (-pkin(3) * t60 + pkin(6) * t12 + t152 * t6 + t154 * t4) - pkin(2) * t60 + qJ(3) * t5 + pkin(1) * (t148 * t5 - t150 * t60), t147 * (-pkin(6) * t2 - pkin(7) * t193 - t152 * t7) + t149 * (-pkin(3) * t49 + pkin(6) * t3 - pkin(7) * t194 + t154 * t7) - pkin(2) * t49 + qJ(3) * t1 + pkin(1) * (t1 * t148 - t150 * t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t82 - t149 * t81, 0, 0, 0, 0, 0, 0, t147 * t72 + t149 * t71, t147 * t78 + t149 * t77, t147 * t74 + t149 * t73, t147 * t28 + t149 * t27, 0, 0, 0, 0, 0, 0, t147 * t21 + t149 * t20, t147 * t31 + t149 * t30, t11 * t149 + t12 * t147, t147 * t3 + t149 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t173, -t130, t91, 0, 0, 0, 0, 0, 0, t104, t106, t84, t80, 0, 0, 0, 0, 0, 0, t41, t203, t60, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t118 - t117, t119, t182, t83, qJDD(4), -t50, -t51, 0, 0, t70, t69, t46, -t70, t163, t142, pkin(4) * t34 - t17, -t188 - t151 * t198 + (-t151 * t205 + t52) * pkin(4), pkin(4) * t24, pkin(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, t46, -t70, t163, t142, -t17, -t18, 0, 0;];
tauJ_reg = t16;
