% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:23:32
% EndTime: 2019-05-05 17:23:37
% DurationCPUTime: 1.86s
% Computational Cost: add. (3780->278), mult. (7819->287), div. (0->0), fcn. (3795->6), ass. (0->179)
t207 = -pkin(7) - pkin(1);
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t129 = qJD(1) ^ 2;
t219 = t123 * t129;
t102 = t126 * t219;
t93 = -t102 + qJDD(3);
t193 = t126 * t93;
t128 = qJD(3) ^ 2;
t119 = t123 ^ 2;
t184 = t119 * t129;
t97 = t128 + t184;
t53 = t123 * t97 - t193;
t227 = t207 * t53;
t92 = t102 + qJDD(3);
t197 = t123 * t92;
t120 = t126 ^ 2;
t183 = t120 * t129;
t99 = t128 + t183;
t52 = t126 * t99 + t197;
t172 = qJD(1) * qJD(3);
t110 = t123 * t172;
t169 = t126 * qJDD(1);
t86 = -0.2e1 * t110 + t169;
t160 = -qJ(2) * t86 + t207 * t52;
t179 = t119 + t120;
t89 = t179 * qJDD(1);
t90 = t179 * t129;
t158 = -qJ(2) * t90 - t207 * t89;
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t178 = qJD(1) * t123;
t77 = qJD(3) * t125 + t122 * t178;
t79 = -qJD(3) * t122 + t125 * t178;
t48 = t79 * t77;
t85 = -t110 + t169;
t75 = qJDD(6) + t85;
t223 = -t48 + t75;
t226 = t122 * t223;
t225 = t125 * t223;
t224 = 0.2e1 * qJD(1);
t209 = pkin(3) + pkin(4);
t162 = t126 * t172;
t170 = t123 * qJDD(1);
t84 = t162 + t170;
t41 = -qJD(6) * t77 - qJDD(3) * t122 + t125 * t84;
t177 = qJD(1) * t126;
t104 = qJD(6) + t177;
t64 = t104 * t77;
t222 = -t64 + t41;
t94 = -qJD(3) * pkin(4) - qJ(5) * t177;
t221 = qJ(5) * t84 + qJD(3) * t94;
t220 = pkin(3) * t99 + qJ(4) * t92;
t154 = t85 + t110;
t218 = t154 * qJ(5);
t217 = (t128 - t97) * qJ(4);
t163 = t207 * t129;
t216 = qJ(5) + t207;
t215 = -pkin(3) * t162 - qJ(4) * t110;
t171 = qJD(2) * qJD(1);
t115 = 0.2e1 * t171;
t165 = -0.2e1 * t177;
t214 = qJD(4) * t165 + t115;
t210 = 0.2e1 * qJD(4);
t124 = sin(qJ(1));
t127 = cos(qJ(1));
t161 = g(1) * t124 - g(2) * t127;
t155 = qJDD(2) - t161;
t144 = -qJ(2) * t129 + t155;
t63 = qJDD(1) * t207 + t144;
t59 = t123 * t63;
t213 = qJDD(3) * qJ(4) + qJD(3) * t210 + t59;
t212 = -pkin(4) * t84 - qJ(5) * t184 + qJDD(5);
t73 = t77 ^ 2;
t74 = t79 ^ 2;
t100 = t104 ^ 2;
t211 = -0.2e1 * qJD(2);
t208 = pkin(3) + pkin(8);
t206 = pkin(3) * t84;
t205 = pkin(1) * t129;
t204 = g(3) * t126;
t203 = t123 * g(3);
t202 = -qJ(4) - pkin(5);
t201 = qJ(4) * t85;
t148 = t213 + t221;
t157 = pkin(5) * t126 - pkin(8) * t123;
t167 = pkin(4) * t184;
t182 = t126 * qJ(4);
t153 = pkin(3) * t123 - t182;
t80 = t153 * qJD(1);
t16 = t167 + t204 - qJDD(3) * pkin(5) - t157 * t219 + t208 * t128 + (-(2 * qJD(5)) + t80) * t178 - t148;
t200 = t122 * t16;
t39 = t48 + t75;
t199 = t122 * t39;
t198 = t123 * t86;
t196 = t125 * t16;
t195 = t125 * t39;
t194 = t126 * t63;
t191 = t128 - t90;
t189 = qJ(4) * t123;
t188 = qJ(5) * t126;
t187 = qJDD(1) * pkin(1);
t186 = t104 * t122;
t185 = t104 * t125;
t181 = t128 * qJ(4);
t180 = t210 + t94;
t176 = qJD(3) * t123;
t175 = qJDD(3) + t93;
t174 = -t177 * t80 - qJDD(4);
t173 = pkin(4) + t208;
t168 = -t74 - t100;
t166 = t126 * t48;
t164 = qJD(5) * t224;
t118 = qJDD(1) * qJ(2);
t156 = t127 * g(1) + g(2) * t124;
t149 = -t118 + t156;
t139 = -t163 + t149;
t136 = t139 + t215;
t10 = -t202 * t85 - t208 * t84 + (-pkin(5) * t176 + t211 + (-pkin(8) * qJD(3) + t180) * t126) * qJD(1) + t136 + t212;
t45 = t194 + t203;
t147 = t45 + t174;
t96 = pkin(4) * t102;
t141 = qJD(5) * t165 - t147 + t96;
t17 = -t85 * qJ(5) + t202 * t128 - t173 * qJDD(3) + (-qJ(5) * t176 - t157 * t177) * qJD(1) + t141;
t5 = -t10 * t125 + t122 * t17;
t83 = 0.2e1 * t162 + t170;
t159 = qJ(2) * t83 - t227;
t6 = t10 * t122 + t125 * t17;
t2 = t122 * t6 - t125 * t5;
t3 = t122 * t5 + t125 * t6;
t1 = -t123 * t16 - t126 * t3;
t142 = -pkin(3) * t128 - t178 * t80 + t213;
t34 = t142 - t204;
t133 = t123 * t164 + t221 + t34;
t19 = t133 - t167;
t134 = -t96 + (t63 + t164) * t126 + t209 * qJDD(3) + t174;
t20 = t134 + t181 + t203 + t218;
t7 = t123 * t19 + t126 * t20;
t46 = -t59 + t204;
t26 = -t123 * t46 + t126 * t45;
t152 = t126 * t83 + t198;
t151 = -t126 * (-t128 + t184) + t197;
t150 = qJDD(3) * t125 + t122 * t84;
t146 = (t126 * t180 + t211) * qJD(1);
t145 = -t123 * t209 - qJ(2) + t182;
t140 = -t123 * t173 - t126 * t202 - qJ(2);
t138 = (-qJD(6) + t104) * t79 - t150;
t137 = -qJDD(3) * pkin(3) - t174 - t181 - t194;
t135 = -t149 - t201 + t206 - t215;
t132 = -t135 + t212;
t131 = (qJD(4) * t126 - qJD(2)) * t224 + t136;
t130 = pkin(7) * t129 + t132 + t205;
t91 = (-t119 + t120) * t129;
t65 = -t144 + t187;
t62 = -t74 + t100;
t61 = t73 - t100;
t60 = t139 - 0.2e1 * t171;
t58 = (t85 - t110) * t126;
t56 = t193 - t123 * (t128 - t183);
t51 = (t84 + t162) * t123;
t47 = t74 - t73;
t42 = -t100 - t73;
t40 = -qJD(6) * t79 - t150;
t37 = -t73 - t74;
t35 = -t137 + t203;
t33 = t64 + t41;
t28 = (qJD(6) + t104) * t79 + t150;
t24 = -t122 * t168 - t195;
t23 = t125 * t168 - t199;
t22 = t125 * t42 - t226;
t21 = t122 * t42 + t225;
t18 = -t163 + t146 + t132;
t15 = t123 * t34 + t126 * t35;
t14 = t122 * t33 + t125 * t138;
t13 = t122 * t138 - t125 * t33;
t12 = t123 * t222 - t126 * t24;
t11 = t123 * t28 - t126 * t22;
t8 = t123 * t37 - t126 * t14;
t4 = [0, 0, 0, 0, 0, qJDD(1), t161, t156, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t155 - 0.2e1 * t187, t115 + 0.2e1 * t118 - t156, pkin(1) * t65 + qJ(2) * (t115 - t149 - t205), t58, -t152, t56, t51, -t151, 0, -t123 * t60 + t159, -t126 * t60 - t160, t158 - t26, -qJ(2) * t60 + t207 * t26, t58, t56, t152, 0, t151, t51, -t83 * t182 - t123 * (t201 + (-t83 - t84) * pkin(3) + t131) + t159, t126 * (qJ(4) * t90 + t137) - t123 * (pkin(3) * t90 + t142) + t158, t126 * (-t206 + (t85 + t86) * qJ(4) + t131) - pkin(3) * t198 + t160, t207 * t15 + (qJ(2) + t153) * (t135 + t214 + t163), t51, -t152, -t151, t58, t56, 0, -t123 * (-qJ(5) * t92 + t209 * t86) + t160 + (qJ(4) * t86 + qJ(5) * t99 + t130 + t146) * t126, -t93 * t188 - t123 * (-qJ(5) * t97 - t177 * t94 - t130 + t214) + t145 * t83 + t227, ((qJ(5) * qJDD(1) - qJD(1) * t80 + t164) * t123 + (t90 - t184) * pkin(4) - t191 * pkin(3) + t148) * t123 + ((t154 + t169) * qJ(5) + t191 * qJ(4) + t134) * t126 - t158, t145 * t18 + t216 * t7, t166 - t123 * (-t125 * t41 + t186 * t79), t126 * t47 - t123 * (t122 * t222 + t125 * t28), t126 * t33 - t123 * (t122 * t62 - t225), -t166 - t123 * (t122 * t40 - t185 * t77), t126 * t138 - t123 * (-t125 * t61 + t199), t126 * t75 - t123 * (-t122 * t79 + t125 * t77) * t104, t126 * (-qJ(5) * t22 - t5) - t123 * (-qJ(5) * t28 + t200) + t207 * t11 + t140 * t21, t126 * (-qJ(5) * t24 - t6) - t123 * (-qJ(5) * t222 + t196) + t207 * t12 + t140 * t23, -t14 * t188 - t123 * (-qJ(5) * t37 + t2) + t207 * t8 + t140 * t13, t1 * t216 + t140 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t129, -t65, 0, 0, 0, 0, 0, 0, -t53, -t52, -t89, t26, 0, 0, 0, 0, 0, 0, -t53, -t89, t52, t15, 0, 0, 0, 0, 0, 0, t52, t53, t89, t7, 0, 0, 0, 0, 0, 0, t11, t12, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t91, t169, -t102, -t170, qJDD(3), t45, t46, 0, 0, t102, t169, -t91, qJDD(3), t170, -t102, pkin(3) * t175 + t147 + t217, (-pkin(3) * t126 - t189) * qJDD(1), t34 + t220, pkin(3) * t35 + qJ(4) * t34, -t102, t91, -t170, t102, t169, qJDD(3), (t99 - t184) * pkin(4) + t133 + t220, -t175 * t209 + t141 - t217 - t218, (t126 * t209 + t189) * qJDD(1), qJ(4) * t19 + t20 * t209, -t122 * t41 - t185 * t79, t122 * t28 - t125 * t222, -t125 * t62 - t226, -t125 * t40 - t186 * t77, -t122 * t61 - t195, (t122 * t77 + t125 * t79) * t104, -t173 * t22 - t202 * t28 - t196, -t173 * t24 - t202 * t222 + t200, -t14 * t173 - t202 * t37 - t3, t16 * t202 - t173 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t169, -t99, -t35, 0, 0, 0, 0, 0, 0, -t99, t93, -t169, -t20, 0, 0, 0, 0, 0, 0, t22, t24, t14, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t83, -t90, t18, 0, 0, 0, 0, 0, 0, t21, t23, t13, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, t33, -t48, t138, t75, -t5, -t6, 0, 0;];
tauJ_reg  = t4;
