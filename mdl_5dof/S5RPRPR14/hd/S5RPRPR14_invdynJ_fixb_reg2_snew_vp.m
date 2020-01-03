% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR14_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:17
% EndTime: 2019-12-31 18:35:24
% DurationCPUTime: 2.50s
% Computational Cost: add. (6166->277), mult. (13563->379), div. (0->0), fcn. (8839->8), ass. (0->177)
t136 = sin(pkin(8));
t137 = cos(pkin(8));
t140 = sin(qJ(3));
t143 = cos(qJ(3));
t184 = t143 * t136;
t115 = (t140 * t137 + t184) * qJD(1);
t180 = qJD(1) * t143;
t117 = -t136 * t140 * qJD(1) + t137 * t180;
t189 = t117 * t115;
t209 = qJDD(3) - t189;
t211 = t136 * t209;
t210 = t137 * t209;
t179 = qJD(3) * t115;
t172 = qJD(1) * qJD(3);
t163 = t143 * t172;
t170 = t140 * qJDD(1);
t121 = -t163 - t170;
t129 = t143 * qJDD(1);
t164 = t140 * t172;
t122 = t129 - t164;
t92 = t136 * t121 + t137 * t122;
t79 = t92 - t179;
t139 = sin(qJ(5));
t142 = cos(qJ(5));
t97 = -t142 * qJD(3) + t139 * t117;
t99 = t139 * qJD(3) + t142 * t117;
t76 = t99 * t97;
t160 = -t137 * t121 + t136 * t122;
t89 = qJDD(5) + t160;
t206 = -t76 + t89;
t208 = t139 * t206;
t207 = t142 * t206;
t178 = qJD(3) * t117;
t77 = t160 + t178;
t141 = sin(qJ(1));
t144 = cos(qJ(1));
t157 = t141 * g(1) - t144 * g(2);
t153 = qJDD(2) - t157;
t146 = qJD(1) ^ 2;
t182 = t146 * qJ(2);
t151 = t153 - t182;
t203 = pkin(6) + pkin(1);
t149 = -t203 * qJDD(1) + t151;
t148 = t143 * t149;
t183 = t143 * t146;
t147 = t148 - t122 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t183 - qJ(4) * t172 + g(3)) * t140;
t152 = qJD(3) * pkin(3) - qJ(4) * t180;
t133 = t140 ^ 2;
t188 = t133 * t146;
t94 = t143 * g(3) - t140 * t149;
t72 = -pkin(3) * t188 + t121 * qJ(4) - qJD(3) * t152 - t94;
t41 = -0.2e1 * qJD(4) * t115 + t136 * t147 + t137 * t72;
t111 = qJD(5) + t115;
t161 = -t142 * qJDD(3) + t139 * t92;
t47 = (qJD(5) - t111) * t99 + t161;
t132 = qJDD(1) * qJ(2);
t158 = t144 * g(1) + t141 * g(2);
t154 = -t132 + t158;
t205 = -t121 * pkin(3) - (qJ(4) * t133 + t203) * t146 + t152 * t180 + qJDD(4) - t154;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t110 = t111 ^ 2;
t113 = t115 ^ 2;
t114 = t117 ^ 2;
t204 = 0.2e1 * qJD(4);
t171 = qJD(2) * qJD(1);
t167 = -0.2e1 * t171;
t73 = t167 - t205;
t201 = t136 * t73;
t87 = qJDD(3) + t189;
t200 = t136 * t87;
t199 = t137 * t73;
t198 = t137 * t87;
t145 = qJD(3) ^ 2;
t162 = t136 * t72 - t137 * t147;
t84 = t115 * pkin(4) - t117 * pkin(7);
t30 = -qJDD(3) * pkin(4) - t145 * pkin(7) + (t204 + t84) * t117 + t162;
t197 = t139 * t30;
t57 = t76 + t89;
t196 = t139 * t57;
t195 = t142 * t30;
t194 = t142 * t57;
t40 = t117 * t204 + t162;
t20 = t136 * t41 - t137 * t40;
t193 = t143 * t20;
t192 = qJDD(1) * pkin(1);
t191 = t111 * t139;
t190 = t111 * t142;
t134 = t143 ^ 2;
t187 = t134 * t146;
t166 = t140 * t183;
t186 = t140 * (qJDD(3) + t166);
t185 = t143 * (qJDD(3) - t166);
t181 = t133 + t134;
t177 = qJD(3) * t136;
t176 = qJD(3) * t137;
t173 = qJD(5) + t111;
t169 = t136 * t76;
t168 = t137 * t76;
t165 = -pkin(4) * t137 - pkin(3);
t31 = -t145 * pkin(4) + qJDD(3) * pkin(7) - t115 * t84 + t41;
t131 = 0.2e1 * t171;
t34 = t77 * pkin(4) - t79 * pkin(7) + t131 + t205;
t14 = t139 * t31 - t142 * t34;
t15 = t139 * t34 + t142 * t31;
t6 = t139 * t14 + t142 * t15;
t21 = t136 * t40 + t137 * t41;
t3 = t136 * t6 - t137 * t30;
t1 = t140 * (t136 * t30 + t137 * t6) + t143 * t3;
t5 = t139 * t15 - t142 * t14;
t93 = t140 * g(3) + t148;
t67 = -t140 * t94 + t143 * t93;
t155 = -t139 * qJDD(3) - t142 * t92;
t78 = -t160 + t178;
t64 = -t97 * qJD(5) - t155;
t124 = t181 * qJDD(1);
t123 = t129 - 0.2e1 * t164;
t120 = 0.2e1 * t163 + t170;
t112 = -t151 + t192;
t106 = -t114 - t145;
t105 = -t114 + t145;
t104 = t113 - t145;
t103 = t203 * t146 + t154 + t167;
t101 = -t186 + t143 * (-t145 - t187);
t100 = t140 * (-t145 - t188) + t185;
t85 = -t145 - t113;
t83 = t111 * t97;
t82 = -t96 + t110;
t81 = t95 - t110;
t80 = t179 + t92;
t75 = -t113 - t114;
t74 = t96 - t95;
t71 = -t96 - t110;
t70 = -t136 * t106 - t198;
t69 = t137 * t106 - t200;
t63 = -t99 * qJD(5) - t161;
t62 = -t110 - t95;
t61 = t95 + t96;
t60 = t137 * t85 - t211;
t59 = t136 * t85 + t210;
t55 = (t139 * t99 - t142 * t97) * t111;
t54 = t136 * t80 + t137 * t78;
t53 = t136 * t78 - t137 * t80;
t52 = t173 * t97 + t155;
t51 = t64 + t83;
t50 = t64 - t83;
t48 = -t173 * t99 - t161;
t46 = t142 * t64 - t99 * t191;
t45 = -t139 * t63 + t97 * t190;
t44 = t140 * t70 + t143 * t69;
t43 = t142 * t81 - t196;
t42 = -t139 * t82 + t207;
t38 = -t139 * t71 - t194;
t37 = t142 * t71 - t196;
t36 = t142 * t62 - t208;
t35 = t139 * t62 + t207;
t32 = t140 * t60 + t143 * t59;
t29 = t140 * t54 + t143 * t53;
t28 = t139 * t51 - t142 * t47;
t27 = -t139 * t50 + t142 * t48;
t26 = -t139 * t47 - t142 * t51;
t25 = -t136 * t52 + t137 * t38;
t24 = t136 * t38 + t137 * t52;
t23 = -t136 * t48 + t137 * t36;
t22 = t136 * t36 + t137 * t48;
t19 = -t136 * t61 + t137 * t28;
t18 = t136 * t28 + t137 * t61;
t17 = -pkin(7) * t37 + t195;
t16 = -pkin(7) * t35 + t197;
t12 = -pkin(4) * t37 + t15;
t11 = -pkin(4) * t35 + t14;
t10 = t140 * t25 + t143 * t24;
t9 = t140 * t23 + t143 * t22;
t8 = t140 * t21 + t193;
t7 = t140 * t19 + t143 * t18;
t2 = -pkin(7) * t26 - t5;
t4 = [0, 0, 0, 0, 0, qJDD(1), t157, t158, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t153 - 0.2e1 * t192, t131 + 0.2e1 * t132 - t158, pkin(1) * t112 + qJ(2) * (-t146 * pkin(1) + t131 - t154), (t122 - t164) * t143, -t143 * t120 - t140 * t123, t185 - t140 * (t145 - t187), (-t121 + t163) * t140, t143 * (-t145 + t188) - t186, 0, qJ(2) * t120 - t203 * t100 - t140 * t103, qJ(2) * t123 - t203 * t101 - t143 * t103, t203 * t124 - t181 * t182 - t67, -qJ(2) * t103 - t203 * t67, t143 * (-t117 * t177 + t137 * t92) - t140 * (t117 * t176 + t136 * t92), t143 * (-t136 * t79 - t137 * t77) - t140 * (-t136 * t77 + t137 * t79), t143 * (-t136 * t105 + t210) - t140 * (t137 * t105 + t211), t143 * (t115 * t176 + t136 * t160) - t140 * (t115 * t177 - t137 * t160), t143 * (t137 * t104 - t200) - t140 * (t136 * t104 + t198), (t143 * (-t115 * t137 + t117 * t136) - t140 * (-t115 * t136 - t117 * t137)) * qJD(3), t143 * (-qJ(4) * t59 - t201) - t140 * (-pkin(3) * t77 + qJ(4) * t60 + t199) + qJ(2) * t77 - t203 * t32, t143 * (-qJ(4) * t69 - t199) - t140 * (-pkin(3) * t79 + qJ(4) * t70 - t201) + qJ(2) * t79 - t203 * t44, t143 * (-qJ(4) * t53 - t20) - t140 * (-pkin(3) * t75 + qJ(4) * t54 + t21) + qJ(2) * t75 - t203 * t29, -qJ(4) * t193 - t140 * (pkin(3) * t73 + qJ(4) * t21) - qJ(2) * t73 - t203 * t8, t143 * (t137 * t46 + t169) - t140 * (t136 * t46 - t168), t143 * (t136 * t74 + t137 * t27) - t140 * (t136 * t27 - t137 * t74), t143 * (t136 * t51 + t137 * t42) - t140 * (t136 * t42 - t137 * t51), t143 * (t137 * t45 - t169) - t140 * (t136 * t45 + t168), t143 * (-t136 * t47 + t137 * t43) - t140 * (t136 * t43 + t137 * t47), t143 * (t136 * t89 + t137 * t55) - t140 * (t136 * t55 - t137 * t89), t143 * (-qJ(4) * t22 - t136 * t11 + t137 * t16) - t140 * (-pkin(3) * t35 + qJ(4) * t23 + t137 * t11 + t136 * t16) + qJ(2) * t35 - t203 * t9, t143 * (-qJ(4) * t24 - t136 * t12 + t137 * t17) - t140 * (-pkin(3) * t37 + qJ(4) * t25 + t137 * t12 + t136 * t17) + qJ(2) * t37 - t203 * t10, t143 * (-qJ(4) * t18 + t137 * t2) - t140 * (qJ(4) * t19 + t136 * t2) - t203 * t7 + (pkin(4) * t184 - t140 * t165 + qJ(2)) * t26, (t143 * (pkin(4) * t136 - pkin(7) * t137) - t140 * (-pkin(7) * t136 + t165) + qJ(2)) * t5 + (-t203 - qJ(4)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t146, -t112, 0, 0, 0, 0, 0, 0, t100, t101, -t124, t67, 0, 0, 0, 0, 0, 0, t32, t44, t29, t8, 0, 0, 0, 0, 0, 0, t9, t10, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, (-t133 + t134) * t146, t129, -t166, -t170, qJDD(3), t93, t94, 0, 0, t189, t114 - t113, t80, -t189, t78, qJDD(3), pkin(3) * t59 - t40, pkin(3) * t69 - t41, pkin(3) * t53, pkin(3) * t20, t139 * t64 + t190 * t99, t139 * t48 + t142 * t50, t142 * t82 + t208, t142 * t63 + t191 * t97, t139 * t81 + t194, (-t139 * t97 - t142 * t99) * t111, pkin(3) * t22 + pkin(4) * t48 + pkin(7) * t36 - t195, pkin(3) * t24 + pkin(4) * t52 + pkin(7) * t38 + t197, pkin(3) * t18 + pkin(4) * t61 + pkin(7) * t28 + t6, pkin(3) * t3 - pkin(4) * t30 + pkin(7) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t79, t75, -t73, 0, 0, 0, 0, 0, 0, t35, t37, t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t74, t51, -t76, -t47, t89, -t14, -t15, 0, 0;];
tauJ_reg = t4;
