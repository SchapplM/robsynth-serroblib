% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:43
% EndTime: 2019-12-31 20:53:49
% DurationCPUTime: 1.85s
% Computational Cost: add. (1655->309), mult. (2330->334), div. (0->0), fcn. (1163->8), ass. (0->180)
t112 = qJDD(1) + qJDD(2);
t120 = sin(qJ(3));
t103 = t120 * qJ(4);
t123 = cos(qJ(3));
t188 = t123 * pkin(3) + t103;
t225 = -pkin(2) - t188;
t227 = t112 * t225;
t113 = qJD(1) + qJD(2);
t226 = t113 * t225;
t224 = pkin(2) + t103;
t191 = t113 * t123;
t217 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t208 = pkin(3) + qJ(5);
t223 = t208 * qJD(3);
t222 = t208 * qJDD(3);
t179 = qJD(3) * t123;
t166 = t113 * t179;
t83 = t120 * t112;
t221 = (-t166 - t83) * pkin(4);
t121 = sin(qJ(2));
t176 = qJDD(1) * t121;
t124 = cos(qJ(2));
t183 = qJD(2) * t124;
t36 = t112 * pkin(7) + (qJD(1) * t183 + t176) * pkin(1);
t30 = t120 * t36;
t201 = pkin(1) * qJD(1);
t172 = t121 * t201;
t55 = t113 * pkin(7) + t172;
t41 = t55 * t179;
t169 = qJDD(4) + t30 + t41;
t197 = qJDD(3) * pkin(3);
t10 = t169 - t197;
t180 = qJD(3) * t120;
t31 = t123 * t36;
t9 = t55 * t180 - t217 - t31;
t220 = t10 * t120 - t9 * t123;
t118 = qJ(1) + qJ(2);
t101 = sin(t118);
t102 = cos(t118);
t203 = g(1) * t102 + g(2) * t101;
t182 = qJD(3) * qJ(4);
t157 = -qJD(5) - t182;
t46 = t123 * t55;
t28 = pkin(4) * t191 + t46;
t19 = -t157 + t28;
t45 = t120 * t55;
t219 = -qJD(4) - t45;
t218 = 0.2e1 * t217;
t194 = t102 * t120;
t196 = t101 * t120;
t216 = g(1) * t194 + g(2) * t196 - g(3) * t123;
t126 = qJD(3) ^ 2;
t215 = pkin(7) * t126;
t92 = g(1) * t101;
t122 = sin(qJ(1));
t214 = g(1) * t122;
t213 = g(3) * t120;
t212 = t112 * pkin(2);
t211 = t113 * pkin(2);
t210 = t124 * pkin(1);
t207 = g(1) * t196 - g(2) * t194;
t193 = t102 * t123;
t195 = t101 * t123;
t206 = g(1) * t195 - g(2) * t193;
t84 = t123 * t112;
t205 = qJ(4) * t84 + qJD(4) * t191;
t88 = t102 * pkin(7);
t204 = t102 * pkin(4) + t88;
t202 = -qJD(2) * t172 + qJDD(1) * t210;
t95 = t121 * pkin(1) + pkin(7);
t200 = t126 * t95;
t33 = -qJD(3) * pkin(3) - t219;
t199 = t33 * t123;
t198 = pkin(7) * qJDD(3);
t192 = t113 * t120;
t190 = t123 * qJ(5);
t27 = -pkin(4) * t192 - t45;
t189 = qJD(4) - t27;
t116 = t120 ^ 2;
t117 = t123 ^ 2;
t187 = t116 - t117;
t186 = t116 + t117;
t185 = qJD(1) * t124;
t184 = qJD(2) * t121;
t181 = qJD(3) * t113;
t178 = qJDD(3) * t95;
t177 = t120 * qJD(4);
t158 = t113 * t172;
t53 = t120 * t158;
t175 = t53 + t207;
t98 = pkin(3) * t180;
t174 = t113 * t98 - t202;
t173 = pkin(1) * t183;
t171 = pkin(1) * t185;
t111 = t113 ^ 2;
t170 = t120 * t111 * t123;
t168 = pkin(4) * t84 + qJDD(5) + t31;
t96 = -pkin(2) - t210;
t167 = t113 * t184;
t165 = -pkin(4) * t113 - t55;
t21 = -t171 + t226;
t164 = -t21 - t226;
t35 = -t202 - t212;
t56 = -t171 - t211;
t163 = t35 * t120 + t56 * t179 - t207;
t162 = pkin(3) * t193 + t101 * pkin(7) + t224 * t102;
t161 = -t123 * t158 - t171 * t180 - t206;
t160 = t186 * t112;
t159 = -t30 + t216;
t156 = -t212 + t215;
t154 = t188 + t190;
t153 = t165 * qJD(3);
t152 = -qJDD(4) + t159;
t38 = -t46 - t182;
t150 = t38 * t120 + t199;
t149 = t120 * t33 - t123 * t38;
t148 = -g(2) * t102 + t202 + t92;
t147 = t101 * pkin(4) + t102 * t190 + t162;
t146 = -qJ(4) * t123 + qJ(5) * t120;
t44 = -pkin(2) - t154;
t145 = t186 * t171;
t138 = -t208 * t123 - t224;
t1 = (t146 * qJD(3) - t123 * qJD(5) - t177) * t113 + t138 * t112 + t174;
t100 = pkin(1) * t184;
t20 = qJ(5) * t180 + t157 * t123 - t177 + t98;
t16 = t100 + t20;
t34 = t44 - t210;
t143 = -t112 * t34 - t113 * t16 - t1;
t142 = -t112 * t44 - t113 * t20 - t1;
t141 = -t41 + t152;
t140 = -qJ(4) * t179 - t177;
t139 = t33 * t179 + t38 * t180 - t203 + t220;
t39 = t140 + t98;
t137 = t113 * t39 + t215 + t227;
t29 = t100 + t39;
t47 = t96 - t188;
t136 = t112 * t47 + t113 * t29 + t200;
t17 = t189 - t223;
t134 = t169 - t221 - t222;
t4 = -qJD(3) * qJD(5) + t134;
t6 = t120 * t153 + t168 + t217;
t135 = t4 * t120 + t6 * t123 + t17 * t179 - t19 * t180 - t203;
t133 = pkin(1) * t167 + t112 * t96 + t200;
t132 = t138 * t92;
t131 = -g(1) * t88 - t225 * t92;
t14 = t138 * t113 - t171;
t130 = -t123 * t203 + t14 * t191 + t168;
t129 = -t178 + (t113 * t96 - t173) * qJD(3);
t128 = t178 + (-t113 * t47 + t173 - t21) * qJD(3);
t127 = t150 * qJD(3) + t220;
t125 = cos(qJ(1));
t109 = t125 * pkin(1);
t107 = t123 * pkin(4);
t106 = t120 * pkin(4);
t99 = pkin(4) * t179;
t82 = pkin(3) * t192;
t69 = qJ(4) * t193;
t66 = qJ(4) * t195;
t63 = t123 * pkin(7) + t107;
t62 = t120 * pkin(7) + t106;
t60 = -t116 * t111 - t126;
t59 = qJDD(3) * t123 - t126 * t120;
t58 = qJDD(3) * t120 + t126 * t123;
t52 = qJDD(3) + t170;
t51 = pkin(7) * t179 + t99;
t50 = (-pkin(4) - pkin(7)) * t180;
t49 = t123 * t95 + t107;
t48 = t120 * t95 + t106;
t42 = t56 * t180;
t40 = -qJ(4) * t191 + t82;
t37 = t116 * t112 + 0.2e1 * t120 * t166;
t24 = t146 * t113 + t82;
t23 = t120 * t173 + t95 * t179 + t99;
t22 = t123 * t173 + (-pkin(4) - t95) * t180;
t18 = 0.2e1 * t120 * t84 - 0.2e1 * t187 * t181;
t15 = t21 * t192;
t12 = t14 * t180;
t7 = t140 * t113 + t174 + t227;
t5 = t7 * t123;
t2 = [qJDD(1), -g(2) * t125 + t214, g(1) * t125 + g(2) * t122, t112, (t112 * t124 - t167) * pkin(1) + t148, ((-qJDD(1) - t112) * t121 + (-qJD(1) - t113) * t183) * pkin(1) + t203, t37, t18, t58, t59, 0, t42 + t129 * t120 + (-t133 - t35) * t123 + t206, t133 * t120 + t129 * t123 + t163, t186 * t113 * t173 + t95 * t160 + t139, t128 * t120 + t136 * t123 - t206 + t5, t128 * t123 + (-t136 - t7) * t120 + t207, t7 * t47 + t21 * t29 - g(2) * (t109 + t162) + (t149 * t183 + t214) * pkin(1) + t127 * t95 + t131, (t120 * t48 + t123 * t49) * t112 + (t120 * t23 + t123 * t22 + (-t120 * t49 + t123 * t48) * qJD(3)) * t113 + t135, t49 * qJDD(3) + t143 * t120 + (t22 + (-t113 * t34 - t14) * t123) * qJD(3) + t207, -t48 * qJDD(3) + t12 + (t34 * t192 - t23) * qJD(3) + t143 * t123 + t206, t1 * t34 + t14 * t16 + t4 * t48 + t17 * t23 + t6 * t49 + t19 * t22 - g(1) * (-t122 * pkin(1) + t204) - g(2) * (t109 + t147) - t132; 0, 0, 0, t112, t148 + t158, (-t176 + (-qJD(2) + t113) * t185) * pkin(1) + t203, t37, t18, t58, t59, 0, t42 + (-pkin(2) * t181 - t198) * t120 + (-t156 - t35) * t123 - t161, -t53 + t156 * t120 + (-t198 + (t171 - t211) * qJD(3)) * t123 + t163, pkin(7) * t160 - t113 * t145 + t139, t5 + t137 * t123 + (t164 * qJD(3) + t198) * t120 + t161, (t198 + (t164 - t171) * qJD(3)) * t123 + (-t137 - t7) * t120 + t175, t7 * t225 + t21 * t39 - g(2) * t162 + (-t121 * t21 - t149 * t124) * t201 + t127 * pkin(7) + t131, (t120 * t62 + t123 * t63) * t112 + (t120 * t51 + t123 * t50 + (-t120 * t63 + t123 * t62) * qJD(3) - t145) * t113 + t135, t63 * qJDD(3) + t142 * t120 + (t50 + (-t113 * t44 - t14 - t171) * t123) * qJD(3) + t175, -t62 * qJDD(3) + t12 + (t44 * t192 - t51) * qJD(3) + t142 * t123 - t161, t1 * t44 + t14 * t20 + t4 * t62 + t17 * t51 + t6 * t63 + t19 * t50 - g(1) * t204 - g(2) * t147 - t132 + (-t121 * t14 + (-t120 * t17 - t123 * t19) * t124) * t201; 0, 0, 0, 0, 0, 0, -t170, t187 * t111, t83, t84, qJDD(3), -t56 * t192 + t159, t213 - t31 + (-t113 * t56 + t203) * t123, -pkin(3) * t83 + (-qJD(3) * t188 - t150) * t113 + t205, -t40 * t191 + t15 - t152 - 0.2e1 * t197, t31 + (t113 * t40 - g(3)) * t120 + (t113 * t21 - t203) * t123 + t218, -t9 * qJ(4) - t10 * pkin(3) - t21 * t40 - t55 * t199 - g(1) * (-pkin(3) * t194 + t69) - g(2) * (-pkin(3) * t196 + t66) - g(3) * t188 + t219 * t38, -t208 * t83 + (-t17 - t27 - t223) * t191 + t205, -t27 * qJD(3) + (t113 * t24 - g(3) + t153) * t120 + t130 + t218, (-t120 * t14 + t123 * t24) * t113 + (0.2e1 * qJD(5) + t28) * qJD(3) + 0.2e1 * t222 + t141 + t221, t6 * qJ(4) - t14 * t24 - g(1) * t69 - g(2) * t66 - g(3) * t154 + t189 * t19 + (-qJD(5) - t28) * t17 + (t203 * t120 - t4) * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t52, t60, qJD(3) * t38 - t141 + t15 - t197, t83, t60, -t52, t14 * t192 + (-qJD(5) - t19) * qJD(3) + t134 - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, qJDD(3) - t170, -t111 * t117 - t126, -t213 + (t120 * t165 + t17) * qJD(3) + t130 + t217;];
tau_reg = t2;
