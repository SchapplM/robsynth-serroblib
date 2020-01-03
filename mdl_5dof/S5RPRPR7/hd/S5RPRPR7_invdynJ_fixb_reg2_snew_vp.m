% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:52
% EndTime: 2019-12-31 18:19:57
% DurationCPUTime: 2.13s
% Computational Cost: add. (7065->306), mult. (15361->436), div. (0->0), fcn. (10162->10), ass. (0->178)
t149 = sin(pkin(9));
t151 = cos(pkin(9));
t159 = cos(qJ(3));
t192 = qJD(1) * t159;
t156 = sin(qJ(3));
t193 = qJD(1) * t156;
t113 = t149 * t193 - t151 * t192;
t115 = t149 * t192 + t151 * t193;
t94 = t115 * t113;
t212 = qJDD(3) - t94;
t217 = t149 * t212;
t216 = t151 * t212;
t155 = sin(qJ(5));
t158 = cos(qJ(5));
t101 = qJD(3) * t155 + t115 * t158;
t99 = -t158 * qJD(3) + t115 * t155;
t77 = t101 * t99;
t184 = qJD(1) * qJD(3);
t176 = t159 * t184;
t183 = t156 * qJDD(1);
t123 = t176 + t183;
t140 = t159 * qJDD(1);
t177 = t156 * t184;
t124 = t140 - t177;
t171 = t123 * t149 - t151 * t124;
t93 = qJDD(5) + t171;
t213 = -t77 + t93;
t215 = t155 * t213;
t214 = t158 * t213;
t190 = qJD(3) * t115;
t78 = t171 + t190;
t146 = -g(3) + qJDD(2);
t162 = qJD(1) ^ 2;
t157 = sin(qJ(1));
t160 = cos(qJ(1));
t170 = g(1) * t160 + g(2) * t157;
t121 = -pkin(1) * t162 - t170;
t150 = sin(pkin(8));
t152 = cos(pkin(8));
t169 = g(1) * t157 - g(2) * t160;
t165 = qJDD(1) * pkin(1) + t169;
t194 = t152 * t121 + t150 * t165;
t166 = -pkin(2) * t162 + qJDD(1) * pkin(6) + t194;
t164 = t156 * t166;
t196 = t156 * t162;
t163 = -t164 - t123 * qJ(4) + qJDD(3) * pkin(3) + (pkin(3) * t196 + qJ(4) * t184 + t146) * t159;
t128 = qJD(3) * pkin(3) - qJ(4) * t193;
t145 = t159 ^ 2;
t143 = t145 * t162;
t83 = t156 * t146 + t159 * t166;
t63 = -pkin(3) * t143 + t124 * qJ(4) - qJD(3) * t128 + t83;
t34 = -0.2e1 * qJD(4) * t113 + t149 * t163 + t151 * t63;
t110 = qJD(5) + t113;
t96 = t123 * t151 + t124 * t149;
t174 = -t158 * qJDD(3) + t155 * t96;
t47 = (qJD(5) - t110) * t101 + t174;
t97 = t99 ^ 2;
t98 = t101 ^ 2;
t109 = t110 ^ 2;
t111 = t113 ^ 2;
t112 = t115 ^ 2;
t211 = 0.2e1 * qJD(4);
t210 = pkin(4) * t149;
t172 = -t150 * t121 + t152 * t165;
t87 = -qJDD(1) * pkin(2) - t162 * pkin(6) - t172;
t67 = -t124 * pkin(3) - qJ(4) * t143 + t128 * t193 + qJDD(4) + t87;
t208 = t149 * t67;
t91 = qJDD(3) + t94;
t207 = t149 * t91;
t206 = t151 * t67;
t205 = t151 * t91;
t161 = qJD(3) ^ 2;
t175 = t149 * t63 - t151 * t163;
t88 = pkin(4) * t113 - pkin(7) * t115;
t26 = -qJDD(3) * pkin(4) - t161 * pkin(7) + (t211 + t88) * t115 + t175;
t204 = t155 * t26;
t61 = t77 + t93;
t203 = t155 * t61;
t33 = t115 * t211 + t175;
t18 = t149 * t34 - t151 * t33;
t202 = t156 * t18;
t201 = t158 * t26;
t200 = t158 * t61;
t199 = t110 * t155;
t198 = t110 * t158;
t136 = t159 * t196;
t129 = qJDD(3) + t136;
t197 = t156 * t129;
t130 = qJDD(3) - t136;
t195 = t159 * t130;
t191 = qJD(3) * t113;
t189 = qJD(3) * t149;
t188 = qJD(3) * t151;
t185 = qJD(5) + t110;
t182 = t149 * t77;
t181 = t151 * t77;
t180 = -pkin(1) * t152 - pkin(2);
t179 = pkin(1) * t150 + pkin(6);
t178 = -pkin(4) * t151 - pkin(3);
t27 = -pkin(4) * t161 + qJDD(3) * pkin(7) - t113 * t88 + t34;
t173 = -t96 + t191;
t36 = t78 * pkin(4) + t173 * pkin(7) + t67;
t14 = t155 * t27 - t158 * t36;
t15 = t155 * t36 + t158 * t27;
t6 = t14 * t155 + t158 * t15;
t19 = t149 * t33 + t151 * t34;
t82 = -t159 * t146 + t164;
t55 = t156 * t82 + t159 * t83;
t5 = -t14 * t158 + t15 * t155;
t167 = -t155 * qJDD(3) - t158 * t96;
t79 = -t171 + t190;
t70 = -qJD(5) * t99 - t167;
t144 = t156 ^ 2;
t141 = t144 * t162;
t135 = -t143 - t161;
t134 = -t141 - t161;
t127 = t141 + t143;
t126 = (t144 + t145) * qJDD(1);
t125 = t140 - 0.2e1 * t177;
t122 = 0.2e1 * t176 + t183;
t106 = -t112 - t161;
t105 = -t112 + t161;
t104 = t111 - t161;
t103 = -t134 * t156 - t195;
t102 = t135 * t159 - t197;
t89 = -t161 - t111;
t86 = t110 * t99;
t85 = -t98 + t109;
t84 = t97 - t109;
t81 = t96 + t191;
t76 = -t111 - t112;
t75 = t98 - t97;
t73 = -t98 - t109;
t72 = -t106 * t149 - t205;
t71 = t106 * t151 - t207;
t69 = -qJD(5) * t101 - t174;
t68 = -t109 - t97;
t66 = t97 + t98;
t65 = t151 * t89 - t217;
t64 = t149 * t89 + t216;
t56 = (t101 * t155 - t158 * t99) * t110;
t54 = t149 * t81 + t151 * t79;
t53 = t149 * t79 - t151 * t81;
t52 = t185 * t99 + t167;
t51 = t70 + t86;
t50 = t70 - t86;
t48 = -t185 * t101 - t174;
t46 = -t101 * t199 + t158 * t70;
t45 = -t155 * t69 + t99 * t198;
t44 = -t156 * t71 + t159 * t72;
t43 = t158 * t84 - t203;
t42 = -t155 * t85 + t214;
t41 = -t155 * t73 - t200;
t40 = t158 * t73 - t203;
t39 = t158 * t68 - t215;
t38 = t155 * t68 + t214;
t37 = -t156 * t64 + t159 * t65;
t31 = -t156 * t53 + t159 * t54;
t30 = t155 * t51 - t158 * t47;
t29 = -t155 * t50 + t158 * t48;
t28 = -t155 * t47 - t158 * t51;
t25 = -t149 * t52 + t151 * t41;
t24 = t149 * t41 + t151 * t52;
t23 = -t149 * t48 + t151 * t39;
t22 = t149 * t39 + t151 * t48;
t21 = -t149 * t66 + t151 * t30;
t20 = t149 * t30 + t151 * t66;
t17 = -pkin(7) * t40 + t201;
t16 = -pkin(7) * t38 + t204;
t12 = -t156 * t24 + t159 * t25;
t11 = -t156 * t22 + t159 * t23;
t10 = -pkin(4) * t40 + t15;
t9 = -pkin(4) * t38 + t14;
t7 = t159 * t19 - t202;
t4 = -pkin(7) * t28 - t5;
t3 = t149 * t26 + t151 * t6;
t2 = t149 * t6 - t151 * t26;
t1 = [0, 0, 0, 0, 0, qJDD(1), t169, t170, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t152 - t150 * t162) + t172, pkin(1) * (-qJDD(1) * t150 - t152 * t162) - t194, 0, pkin(1) * (t150 * t194 + t152 * t172), (t123 + t176) * t156, t122 * t159 + t125 * t156, t197 + t159 * (-t141 + t161), (t124 - t177) * t159, t156 * (t143 - t161) + t195, 0, -t159 * t87 + pkin(2) * t125 + pkin(6) * t102 + pkin(1) * (t102 * t150 + t125 * t152), t156 * t87 - pkin(2) * t122 + pkin(6) * t103 + pkin(1) * (t103 * t150 - t122 * t152), pkin(2) * t127 + pkin(6) * t126 + pkin(1) * (t126 * t150 + t127 * t152) + t55, -pkin(2) * t87 + pkin(6) * t55 + pkin(1) * (t150 * t55 - t152 * t87), t156 * (-t115 * t189 + t151 * t96) + t159 * (t115 * t188 + t149 * t96), t156 * (t149 * t173 - t151 * t78) + t159 * (-t149 * t78 - t151 * t173), t156 * (-t105 * t149 + t216) + t159 * (t105 * t151 + t217), t156 * (t113 * t188 + t149 * t171) + t159 * (t113 * t189 - t151 * t171), t156 * (t104 * t151 - t207) + t159 * (t104 * t149 + t205), (t156 * (-t113 * t151 + t115 * t149) + t159 * (-t113 * t149 - t115 * t151)) * qJD(3), t156 * (-qJ(4) * t64 + t208) + t159 * (-pkin(3) * t78 + qJ(4) * t65 - t206) - pkin(2) * t78 + pkin(6) * t37 + pkin(1) * (t150 * t37 - t152 * t78), t156 * (-qJ(4) * t71 + t206) + t159 * (pkin(3) * t173 + qJ(4) * t72 + t208) + pkin(2) * t173 + pkin(6) * t44 + pkin(1) * (t150 * t44 + t152 * t173), t156 * (-qJ(4) * t53 - t18) + t159 * (-pkin(3) * t76 + qJ(4) * t54 + t19) - pkin(2) * t76 + pkin(6) * t31 + pkin(1) * (t150 * t31 - t152 * t76), -qJ(4) * t202 + t159 * (-pkin(3) * t67 + qJ(4) * t19) - pkin(2) * t67 + pkin(6) * t7 + pkin(1) * (t150 * t7 - t152 * t67), t156 * (t151 * t46 + t182) + t159 * (t149 * t46 - t181), t156 * (t149 * t75 + t151 * t29) + t159 * (t149 * t29 - t151 * t75), t156 * (t149 * t51 + t151 * t42) + t159 * (t149 * t42 - t151 * t51), t156 * (t151 * t45 - t182) + t159 * (t149 * t45 + t181), t156 * (-t149 * t47 + t151 * t43) + t159 * (t149 * t43 + t151 * t47), t156 * (t149 * t93 + t151 * t56) + t159 * (t149 * t56 - t151 * t93), t156 * (-qJ(4) * t22 - t149 * t9 + t151 * t16) + t159 * (-pkin(3) * t38 + qJ(4) * t23 + t149 * t16 + t151 * t9) - pkin(2) * t38 + pkin(6) * t11 + pkin(1) * (t11 * t150 - t152 * t38), t156 * (-qJ(4) * t24 - t10 * t149 + t151 * t17) + t159 * (-pkin(3) * t40 + qJ(4) * t25 + t10 * t151 + t149 * t17) - pkin(2) * t40 + pkin(6) * t12 + pkin(1) * (t12 * t150 - t152 * t40), t156 * (-qJ(4) * t20 + t151 * t4) + t159 * (qJ(4) * t21 + t149 * t4) + t179 * (-t156 * t20 + t159 * t21) + (t156 * t210 + t159 * t178 + t180) * t28, (t156 * (-pkin(7) * t151 + t210) + t159 * (-pkin(7) * t149 + t178) + t180) * t5 + (t179 + qJ(4)) * (-t156 * t2 + t159 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, 0, 0, 0, 0, 0, 0, t129 * t159 + t135 * t156, -t130 * t156 + t134 * t159, 0, t156 * t83 - t159 * t82, 0, 0, 0, 0, 0, 0, t156 * t65 + t159 * t64, t156 * t72 + t159 * t71, t156 * t54 + t159 * t53, t156 * t19 + t159 * t18, 0, 0, 0, 0, 0, 0, t156 * t23 + t159 * t22, t156 * t25 + t159 * t24, t156 * t21 + t159 * t20, t156 * t3 + t159 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t141 - t143, t183, t136, t140, qJDD(3), -t82, -t83, 0, 0, t94, t112 - t111, t81, -t94, t79, qJDD(3), pkin(3) * t64 - t33, pkin(3) * t71 - t34, pkin(3) * t53, pkin(3) * t18, t101 * t198 + t155 * t70, t155 * t48 + t158 * t50, t158 * t85 + t215, t158 * t69 + t99 * t199, t155 * t84 + t200, (-t101 * t158 - t155 * t99) * t110, pkin(3) * t22 + pkin(4) * t48 + pkin(7) * t39 - t201, pkin(3) * t24 + pkin(4) * t52 + pkin(7) * t41 + t204, pkin(3) * t20 + pkin(4) * t66 + pkin(7) * t30 + t6, pkin(3) * t2 - pkin(4) * t26 + pkin(7) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t173, t76, t67, 0, 0, 0, 0, 0, 0, t38, t40, t28, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t75, t51, -t77, -t47, t93, -t14, -t15, 0, 0;];
tauJ_reg = t1;
