% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:36
% EndTime: 2019-12-05 16:13:45
% DurationCPUTime: 2.10s
% Computational Cost: add. (1684->343), mult. (3735->451), div. (0->0), fcn. (2577->8), ass. (0->163)
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t119 = cos(qJ(3));
t177 = qJD(3) * t119;
t154 = t118 * t177;
t120 = cos(qJ(2));
t179 = qJD(2) * t120;
t220 = t117 * t179 + t154;
t114 = sin(pkin(7));
t116 = cos(pkin(7));
t145 = g(1) * t116 + g(2) * t114;
t108 = t119 * qJDD(2);
t115 = cos(pkin(8));
t178 = qJD(3) * t117;
t155 = t118 * t178;
t113 = sin(pkin(8));
t185 = t119 * t120;
t67 = t113 * t118 + t115 * t185;
t31 = qJD(2) * t67 - t115 * t155;
t170 = qJD(2) * qJD(3);
t150 = t119 * t170;
t168 = t117 * qJDD(2);
t131 = t150 + t168;
t38 = t113 * qJDD(3) + t115 * t131;
t186 = t118 * t119;
t63 = -t113 * t120 + t115 * t186;
t176 = t113 * qJD(3);
t181 = qJD(2) * t117;
t77 = t115 * t181 + t176;
t219 = qJD(2) * (-t119 * t31 + (qJD(3) * t63 - t120 * t77) * t117) - (t117 * t38 + t177 * t77) * t118 - t63 * t108;
t144 = pkin(3) * t117 - qJ(4) * t119;
t58 = qJD(3) * t144 - t117 * qJD(4);
t218 = -t67 * qJD(1) + t113 * t58;
t73 = t77 ^ 2;
t175 = t115 * qJD(3);
t75 = t113 * t181 - t175;
t217 = -t75 ^ 2 - t73;
t188 = t117 * t120;
t60 = t114 * t188 + t116 * t119;
t64 = -t114 * t119 + t116 * t188;
t147 = g(1) * t64 + g(2) * t60;
t189 = t117 * t118;
t165 = g(3) * t189;
t130 = t147 + t165;
t171 = qJD(1) * qJD(2);
t71 = qJDD(2) * pkin(6) + t118 * qJDD(1) + t120 * t171;
t55 = t117 * t71;
t174 = t118 * qJD(1);
t91 = qJD(2) * pkin(6) + t174;
t22 = -qJDD(3) * pkin(3) + t177 * t91 + qJDD(4) + t55;
t125 = -t130 + t22;
t161 = t113 * t185;
t216 = qJD(1) * t161 + (-t174 + t58) * t115;
t215 = t145 * t118;
t151 = t117 * t170;
t214 = t151 - t108;
t106 = t115 * qJDD(3);
t37 = t113 * t131 - t106;
t213 = t37 * pkin(4) - t38 * qJ(5);
t81 = t117 * t91;
t19 = qJDD(3) * qJ(4) + t119 * t71 + (qJD(4) - t81) * qJD(3);
t137 = pkin(3) * t119 + qJ(4) * t117 + pkin(2);
t104 = t118 * t171;
t148 = -t120 * qJDD(1) + t104;
t9 = qJD(2) * t58 - qJDD(2) * t137 + t148;
t5 = t113 * t9 + t115 * t19;
t209 = pkin(6) * t119;
t208 = pkin(6) * t120;
t205 = g(3) * t120;
t173 = t119 * qJD(5);
t204 = -t173 + (-pkin(6) * t115 + qJ(5)) * t178 + t218;
t160 = pkin(6) * t113 + pkin(4);
t203 = -t160 * t178 - t216;
t164 = pkin(6) * t178;
t202 = t113 * t164 + t216;
t201 = -t115 * t164 + t218;
t172 = t120 * qJD(1);
t48 = -qJD(2) * t137 - t172;
t82 = t119 * t91;
t69 = qJD(3) * qJ(4) + t82;
t12 = t113 * t48 + t115 * t69;
t200 = t145 * t189;
t44 = -t113 * t137 + t115 * t209;
t199 = qJD(2) * pkin(2);
t80 = t144 * qJD(2);
t197 = t115 * t80;
t196 = t115 * t137;
t193 = qJ(4) * t115;
t192 = qJD(5) * t77;
t190 = t115 * t117;
t187 = t118 * t115;
t184 = qJDD(1) - g(3);
t111 = t117 ^ 2;
t112 = t119 ^ 2;
t183 = t111 - t112;
t121 = qJD(3) ^ 2;
t122 = qJD(2) ^ 2;
t182 = t121 + t122;
t180 = qJD(2) * t119;
t169 = qJDD(3) * t117;
t167 = t120 * qJDD(2);
t152 = qJ(4) * t108;
t157 = t113 * t180;
t166 = qJD(4) * t157 + t113 * t152 + t115 * t165;
t163 = t75 * t172;
t162 = t77 * t172;
t159 = qJ(4) * t178;
t158 = t117 * t172;
t4 = -t113 * t19 + t115 * t9;
t149 = pkin(2) * t120 + pkin(3) * t185 + pkin(6) * t118 + qJ(4) * t188;
t61 = t114 * t185 - t116 * t117;
t65 = t114 * t117 + t116 * t185;
t146 = g(1) * t65 + g(2) * t61;
t59 = -qJD(3) * pkin(3) + qJD(4) + t81;
t143 = pkin(4) * t113 - qJ(5) * t115;
t92 = -t172 - t199;
t142 = g(3) * t118 - qJD(2) * t92;
t11 = -t113 * t69 + t115 * t48;
t141 = qJ(4) * t38 + qJD(4) * t77;
t139 = g(3) * (-pkin(3) * t189 + qJ(4) * t186);
t138 = pkin(4) * t108 + qJDD(5) - t4;
t136 = pkin(6) + t143;
t62 = t113 * t186 + t115 * t120;
t30 = -qJD(2) * t187 - t113 * t155 + t120 * t157;
t134 = t30 * t77 - t31 * t75 - t37 * t63 + t38 * t62;
t39 = t62 * t114;
t41 = t62 * t116;
t66 = t161 - t187;
t133 = -g(1) * t41 - g(2) * t39 + g(3) * t66;
t40 = t63 * t114;
t42 = t63 * t116;
t132 = g(1) * t42 + g(2) * t40 - g(3) * t67;
t129 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t121 - t148 - t205;
t128 = -qJD(4) * t115 * t75 - g(3) * t186 - t193 * t37 - t146;
t127 = -pkin(6) * qJDD(3) + (t172 + t92 - t199) * qJD(3);
t126 = t113 * t168 - t106 + (-t77 + t176) * t180;
t124 = t62 * t108 + t37 * t189 + (t119 * t30 - t178 * t62) * qJD(2) + t220 * t75;
t123 = t137 * t215;
t98 = t116 * t208;
t95 = t114 * t208;
t83 = -pkin(4) * t115 - qJ(5) * t113 - pkin(3);
t68 = t113 * t80;
t57 = t64 * pkin(3);
t56 = t60 * pkin(3);
t54 = t75 * t180;
t49 = t136 * t117;
t43 = -t113 * t209 - t196;
t35 = t119 * t160 + t196;
t34 = -qJ(5) * t119 + t44;
t27 = t143 * t180 + t82;
t26 = -t190 * t91 + t68;
t25 = t113 * t81 + t197;
t24 = -qJD(5) * t190 + t136 * t177;
t21 = -t197 + (-pkin(4) * qJD(2) - t113 * t91) * t117;
t20 = t68 + (qJ(5) * qJD(2) - t115 * t91) * t117;
t16 = t54 + t38;
t10 = pkin(4) * t75 - qJ(5) * t77 + t59;
t8 = -qJ(5) * t180 + t12;
t7 = pkin(4) * t180 + qJD(5) - t11;
t3 = -t192 + t22 + t213;
t2 = -pkin(4) * t151 + t138;
t1 = qJ(5) * t214 - qJD(2) * t173 + t5;
t6 = [t184, 0, -t118 * t122 + t167, -qJDD(2) * t118 - t120 * t122, 0, 0, 0, 0, 0, (-0.2e1 * t151 + t108) * t120 + (-t119 * t182 - t169) * t118, (-qJDD(3) * t118 - 0.2e1 * t120 * t170) * t119 + (t118 * t182 - t167) * t117, t124, -t219, t134, t59 * t154 - t11 * t30 + t12 * t31 - t4 * t62 + t5 * t63 - g(3) + (t118 * t22 + t179 * t59) * t117, t124, t134, t219, t1 * t63 + t10 * t220 + t3 * t189 + t2 * t62 + t7 * t30 + t8 * t31 - g(3); 0, qJDD(2), t120 * t184 + t215, -t118 * t184 + t120 * t145, qJDD(2) * t111 + 0.2e1 * t117 * t150, 0.2e1 * t108 * t117 - 0.2e1 * t170 * t183, t119 * t121 + t169, qJDD(3) * t119 - t117 * t121, 0, t127 * t117 + ((t145 + t171) * t118 + t129) * t119, t127 * t119 + (-t129 - t104) * t117 - t200, (-t163 + pkin(6) * t37 + t22 * t113 + (qJD(2) * t43 + t11) * qJD(3)) * t117 + (-qJDD(2) * t43 - t4 + (pkin(6) * t75 + t113 * t59) * qJD(3) - t202 * qJD(2)) * t119 + t132, (-t162 + pkin(6) * t38 + t22 * t115 + (-qJD(2) * t44 - t12) * qJD(3)) * t117 + (qJDD(2) * t44 + t5 + (pkin(6) * t77 + t115 * t59) * qJD(3) + t201 * qJD(2)) * t119 + t133, -t44 * t37 - t43 * t38 - t202 * t77 - t201 * t75 + (-t11 * t115 - t113 * t12) * t177 + (-t113 * t5 - t115 * t4 - t205) * t117 + t200, t5 * t44 + t4 * t43 - t59 * t158 - g(1) * t98 - g(2) * t95 - g(3) * t149 + t201 * t12 + t202 * t11 + (t117 * t22 + t177 * t59) * pkin(6) + t123, t24 * t75 + t49 * t37 + (-t163 + t3 * t113 + (-qJD(2) * t35 - t7) * qJD(3)) * t117 + (qJD(2) * t203 + qJDD(2) * t35 + t10 * t176 + t2) * t119 + t132, -t34 * t37 + t35 * t38 + t203 * t77 - t204 * t75 + (-t113 * t8 + t115 * t7) * t177 + (-t1 * t113 + t115 * t2 - t205) * t117 + t200, -t24 * t77 - t49 * t38 + (t162 - t3 * t115 + (qJD(2) * t34 + t8) * qJD(3)) * t117 + (-qJD(2) * t204 - qJDD(2) * t34 - t10 * t175 - t1) * t119 - t133, t1 * t34 + t3 * t49 + t2 * t35 - g(1) * (-t42 * pkin(4) - t41 * qJ(5) + t98) - g(2) * (-t40 * pkin(4) - t39 * qJ(5) + t95) - g(3) * (pkin(4) * t67 + qJ(5) * t66 + t149) + t204 * t8 + t203 * t7 + (t24 - t158) * t10 + t123; 0, 0, 0, 0, -t117 * t122 * t119, t183 * t122, t168, t108, qJDD(3), t117 * t142 + t147 - t55, (t142 - t71) * t119 + t146, -t75 * t82 - pkin(3) * t37 + (t147 - t22) * t115 + (-t11 * t117 + t119 * t25 + (-t119 * t59 - t159) * t113) * qJD(2) + t166, -pkin(3) * t38 + (qJDD(2) * t193 - t77 * t91) * t119 + t125 * t113 + (t117 * t12 - t119 * t26 + (-t159 + (qJD(4) - t59) * t119) * t115) * qJD(2), t25 * t77 + t26 * t75 + (t11 * t180 + t5) * t115 + (t12 * t180 + t141 - t4) * t113 + t128, -t22 * pkin(3) - t12 * t26 - t11 * t25 - t59 * t82 + g(1) * t57 + g(2) * t56 - t139 + (-t11 * t113 + t115 * t12) * qJD(4) + (-t4 * t113 + t5 * t115 - t146) * qJ(4), t83 * t37 + (-qJD(5) * t113 - t27) * t75 + (t147 - t3) * t115 + (t117 * t7 - t119 * t21 + (-t10 * t119 - t159) * t113) * qJD(2) + t166, t20 * t75 - t21 * t77 + (-t180 * t7 + t1) * t115 + (t180 * t8 + t141 + t2) * t113 + t128, -t115 * t152 + t27 * t77 - t83 * t38 + (t130 - t3 + t192) * t113 + (-t117 * t8 + t119 * t20 + (t159 + (-qJD(4) + t10) * t119) * t115) * qJD(2), t3 * t83 - t8 * t20 - t10 * t27 - t7 * t21 - g(1) * (t65 * qJ(4) - t57) - g(2) * (qJ(4) * t61 - t56) - t139 + (pkin(4) * t130 + qJ(4) * t1 + qJD(4) * t8) * t115 + (t2 * qJ(4) + qJ(5) * t130 + t7 * qJD(4) - t10 * qJD(5)) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t16, t217, t11 * t77 + t12 * t75 + t125, t126, t217, -t16, t8 * t75 + (-qJD(5) - t7) * t77 + t125 + t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t77 - t214, -t54 + t38, -t112 * t122 - t73, t10 * t77 - g(1) * (t113 * t65 - t116 * t187) - g(2) * (t113 * t61 - t114 * t187) - g(3) * t62 + (-pkin(4) * t178 + t119 * t8) * qJD(2) + t138;];
tau_reg = t6;
