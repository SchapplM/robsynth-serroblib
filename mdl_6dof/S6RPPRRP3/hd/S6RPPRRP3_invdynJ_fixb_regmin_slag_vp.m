% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:52
% EndTime: 2019-03-09 02:03:58
% DurationCPUTime: 2.32s
% Computational Cost: add. (2799->357), mult. (4939->432), div. (0->0), fcn. (2936->10), ass. (0->179)
t105 = sin(pkin(9));
t84 = pkin(1) * t105 + qJ(3);
t200 = qJDD(1) * t84;
t111 = cos(qJ(4));
t100 = qJ(1) + pkin(9);
t91 = sin(t100);
t224 = g(1) * t91;
t92 = cos(t100);
t86 = g(2) * t92;
t161 = t86 - t224;
t108 = sin(qJ(4));
t220 = g(3) * t108;
t236 = -t161 * t111 - t220;
t106 = cos(pkin(9));
t87 = -pkin(1) * t106 - pkin(2);
t78 = -pkin(7) + t87;
t60 = qJD(1) * t78 + qJD(3);
t38 = -t108 * qJD(2) + t111 * t60;
t235 = t38 * qJD(4);
t79 = qJD(1) * t108 + qJD(5);
t74 = qJD(1) * t84;
t234 = -t74 * qJD(1) + t161;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t187 = qJD(4) * t108;
t160 = t107 * t187;
t174 = t111 * qJDD(1);
t188 = qJD(4) * t107;
t189 = qJD(1) * t111;
t67 = t110 * t189 + t188;
t197 = qJD(5) * t67;
t27 = -qJD(1) * t160 - t110 * qJDD(4) + t107 * t174 + t197;
t144 = pkin(4) * t108 - t111 * pkin(8);
t233 = qJDD(1) * t87;
t35 = -qJD(4) * pkin(4) - t38;
t181 = t110 * qJD(4);
t65 = t107 * t189 - t181;
t11 = pkin(5) * t65 - qJ(6) * t67 + t35;
t180 = qJD(1) * qJD(4);
t156 = t111 * t180;
t175 = t108 * qJDD(1);
t62 = qJDD(5) + t156 + t175;
t225 = pkin(8) * t62;
t232 = t11 * t79 - t225;
t39 = t111 * qJD(2) + t108 * t60;
t36 = qJD(4) * pkin(8) + t39;
t56 = t144 + t84;
t40 = t56 * qJD(1);
t12 = -t107 * t36 + t110 * t40;
t193 = qJD(6) - t12;
t6 = -pkin(5) * t79 + t193;
t13 = t107 * t40 + t110 * t36;
t7 = qJ(6) * t79 + t13;
t141 = t107 * t6 + t110 * t7;
t179 = qJD(2) * qJD(4);
t170 = -t108 * qJDD(2) - t111 * t179 - t60 * t187;
t57 = qJDD(1) * t78 + qJDD(3);
t16 = -qJDD(4) * pkin(4) - t111 * t57 - t170;
t183 = qJD(5) * t111;
t159 = t107 * t183;
t127 = t181 * t108 + t159;
t26 = qJD(1) * t127 - qJD(5) * t181 - t107 * qJDD(4) - t110 * t174;
t3 = pkin(5) * t27 + qJ(6) * t26 - qJD(6) * t67 + t16;
t231 = qJD(4) * t141 - t3;
t207 = pkin(8) * qJD(5);
t167 = t79 * t207;
t230 = t167 + t236;
t229 = 0.2e1 * qJD(4) * t74 + qJDD(4) * t78;
t185 = qJD(5) * t107;
t205 = t110 * t62;
t131 = -t79 * t185 + t205;
t206 = t107 * t79;
t228 = t108 * (qJD(4) * t67 - t131) - t111 * (t181 * t79 - t26) + qJD(1) * t206;
t227 = t67 ^ 2;
t226 = pkin(5) * t62;
t223 = t7 * t79;
t219 = g(3) * t111;
t217 = t13 * t79;
t216 = t65 * t79;
t215 = t67 * t65;
t214 = t67 * t79;
t22 = t108 * t27;
t186 = qJD(4) * t111;
t51 = t65 * t186;
t213 = t22 + t51;
t145 = pkin(4) * t111 + pkin(8) * t108;
t69 = t145 * qJD(1);
t212 = t107 * t69 + t110 * t38;
t23 = t108 * t26;
t49 = t67 * t186;
t211 = t49 - t23;
t138 = pkin(5) * t107 - qJ(6) * t110;
t210 = -t107 * qJD(6) + t79 * t138 - t39;
t195 = t108 * t110;
t209 = t107 * t56 + t78 * t195;
t194 = t110 * t111;
t208 = g(3) * t195 + t194 * t86;
t63 = qJD(4) * t145 + qJD(3);
t204 = t110 * t63;
t202 = t111 * t67;
t201 = t62 * qJ(6);
t199 = qJD(4) * t65;
t198 = qJD(5) * t65;
t196 = t107 * t108;
t104 = qJDD(2) - g(3);
t103 = t111 ^ 2;
t192 = t108 ^ 2 - t103;
t113 = qJD(4) ^ 2;
t114 = qJD(1) ^ 2;
t191 = -t113 - t114;
t184 = qJD(5) * t110;
t178 = qJD(3) * qJD(1);
t177 = qJDD(4) * t108;
t176 = qJDD(4) * t111;
t173 = t111 * t224;
t15 = qJDD(4) * pkin(8) + t111 * qJDD(2) + t108 * t57 + t235;
t25 = qJD(1) * t63 + qJDD(1) * t144 + t200;
t172 = -t107 * t25 - t110 * t15 - t40 * t184;
t171 = t111 * t78 * t181 + t107 * t63 + t56 * t184;
t112 = cos(qJ(1));
t168 = t112 * pkin(1) + t92 * pkin(2) + t91 * qJ(3);
t166 = t79 * t195;
t50 = t65 * t187;
t165 = t67 * t187;
t164 = t79 * t184;
t163 = t67 * t183;
t109 = sin(qJ(1));
t158 = -t109 * pkin(1) + t92 * qJ(3);
t157 = t107 * t78 - pkin(5);
t155 = t107 * t15 - t110 * t25 + t36 * t184 + t40 * t185;
t153 = -t26 + t198;
t152 = -t27 + t197;
t151 = t79 * t159;
t41 = -t92 * t110 + t196 * t91;
t43 = t110 * t91 + t196 * t92;
t150 = g(1) * t43 + g(2) * t41;
t42 = t107 * t92 + t195 * t91;
t44 = -t107 * t91 + t195 * t92;
t149 = -g(1) * t44 - g(2) * t42;
t148 = g(1) * t92 + g(2) * t91;
t143 = g(1) * t109 - g(2) * t112;
t142 = t107 * t7 - t110 * t6;
t139 = pkin(5) * t110 + qJ(6) * t107;
t136 = pkin(4) + t139;
t134 = t138 - t78;
t133 = -t111 * t62 + t187 * t79;
t132 = t107 * t62 + t164;
t130 = -t36 * t185 - t172;
t129 = qJDD(3) + t233;
t128 = -t148 + t200;
t126 = t57 + t234;
t125 = t35 * t79 - t225;
t124 = g(1) * t41 - g(2) * t43 + t107 * t219 - t155;
t122 = t108 * t161 - t219;
t1 = qJD(6) * t79 + t130 + t201;
t2 = qJDD(6) + t155 - t226;
t121 = -qJD(5) * t142 + t1 * t110 + t2 * t107;
t120 = t11 * t67 + qJDD(6) - t124;
t61 = t178 + t200;
t119 = -t113 * t78 + t128 + t178 + t61;
t118 = -g(1) * t42 + g(2) * t44 - g(3) * t194 + t130;
t117 = qJD(4) * t11 + t121;
t116 = -t62 * t196 - t111 * t27 + t50 + (-t107 * t186 + (-qJD(5) * t108 - qJD(1)) * t110) * t79;
t73 = -t108 * t113 + t176;
t72 = -t111 * t113 - t177;
t59 = t79 * t160;
t46 = t62 * t194;
t37 = t134 * t111;
t31 = pkin(5) * t67 + qJ(6) * t65;
t29 = t108 * t157 - t110 * t56;
t28 = qJ(6) * t108 + t209;
t19 = t27 * t194;
t18 = -pkin(5) * t189 + t107 * t38 - t110 * t69;
t17 = qJ(6) * t189 + t212;
t14 = (qJD(5) * t139 - qJD(6) * t110) * t111 - t134 * t187;
t8 = t216 - t26;
t5 = qJD(5) * t209 + t157 * t186 - t204;
t4 = qJ(6) * t186 + (-t78 * t185 + qJD(6)) * t108 + t171;
t9 = [qJDD(1), t143, g(1) * t112 + g(2) * t109 (t143 + (t105 ^ 2 + t106 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + t161 + 0.2e1 * t233, t128 + 0.2e1 * t178 + t200, t61 * t84 + t74 * qJD(3) + t129 * t87 - g(1) * (-t91 * pkin(2) + t158) - g(2) * t168, qJDD(1) * t103 - 0.2e1 * t108 * t156, -0.2e1 * t108 * t174 + 0.2e1 * t180 * t192, t73, t72, 0, t119 * t108 + t229 * t111, -t229 * t108 + t119 * t111, -t67 * t159 + (-t111 * t26 - t165) * t110, -t19 + (-t163 + t50) * t110 + (t165 + (t26 + t198) * t111) * t107, -t127 * t79 + t211 + t46, -t22 + t59 + (-t132 - t199) * t111, t108 * t62 + t186 * t79 (-t185 * t56 + t204) * t79 + t56 * t205 + (-t35 * t188 + (-t132 + t199) * t78 - t155) * t108 + (t35 * t184 + t16 * t107 - t78 * t27 + (-t206 * t78 + t12) * qJD(4)) * t111 + t149, -t171 * t79 - t209 * t62 + ((t78 * t79 + t36) * t185 + (-t110 * t35 + t67 * t78) * qJD(4) + t172) * t108 + (-qJD(4) * t13 + t16 * t110 - t185 * t35 + t26 * t78) * t111 + t150, t14 * t65 + t37 * t27 - t29 * t62 - t5 * t79 + (-t11 * t188 - t2) * t108 + (-qJD(4) * t6 + t3 * t107 + t11 * t184) * t111 + t149, -t29 * t26 - t28 * t27 - t4 * t65 + t5 * t67 + t142 * t187 + (-qJD(5) * t141 - t1 * t107 + t110 * t2 + t148) * t111, -t14 * t67 + t37 * t26 + t28 * t62 + t4 * t79 + (t11 * t181 + t1) * t108 + (qJD(4) * t7 + t11 * t185 - t3 * t110) * t111 - t150, t1 * t28 + t7 * t4 + t3 * t37 + t11 * t14 + t2 * t29 + t6 * t5 - g(1) * (t44 * pkin(5) + t43 * qJ(6) + t144 * t92 + t158) - g(2) * (t42 * pkin(5) + t92 * pkin(7) + t41 * qJ(6) + t168) + (-g(1) * (-pkin(2) - pkin(7)) - g(2) * t144) * t91; 0, 0, 0, t104, 0, 0, t104, 0, 0, 0, 0, 0, t72, -t73, 0, 0, 0, 0, 0, -t111 * t132 + t213 + t59, t110 * t133 + t151 + t211, t107 * t133 - t111 * t164 + t213, -t19 + (t163 + t50) * t110 + (t111 * t153 - t165) * t107, -t151 + t23 + t46 + (-t166 - t202) * qJD(4), -t231 * t108 + t117 * t111 - g(3); 0, 0, 0, 0, qJDD(1), -t114, t129 + t234, 0, 0, 0, 0, 0, t108 * t191 + t176, t111 * t191 - t177, 0, 0, 0, 0, 0, t116, t228, t116 (qJD(1) * t67 + t108 * t152 - t51) * t110 + (qJD(1) * t65 + t108 * t153 + t49) * t107, -t228, -t142 * qJD(1) + t117 * t108 + t231 * t111 + t161; 0, 0, 0, 0, 0, 0, 0, t111 * t114 * t108, -t192 * t114, t174, -t175, qJDD(4), t39 * qJD(4) + t111 * t126 + t170 + t220, t235 + (-qJD(4) * t60 - t104) * t111 + (-t126 + t179) * t108, -t26 * t107 + t110 * t214 (-t26 - t216) * t110 + (-t27 - t214) * t107 (t166 - t202) * qJD(1) + t132 (t111 * t65 - t196 * t79) * qJD(1) + t131, -t79 * t189, -t12 * t189 - pkin(4) * t27 - t39 * t65 + (-t173 - t16 + (-t69 - t207) * t79) * t110 + (t38 * t79 + t125) * t107 + t208, pkin(4) * t26 + t212 * t79 + t13 * t189 - t39 * t67 + t125 * t110 + (t16 + t230) * t107, t6 * t189 + t18 * t79 - t136 * t27 + t210 * t65 + (-t167 - t3 - t173) * t110 + t232 * t107 + t208, t17 * t65 - t18 * t67 + (pkin(8) * t152 + t6 * t79 + t1) * t110 + (pkin(8) * t153 + t2 - t223) * t107 + t122, -t7 * t189 - t17 * t79 - t136 * t26 - t210 * t67 - t232 * t110 + (-t230 - t3) * t107, -t7 * t17 - t6 * t18 + t210 * t11 + (t121 + t122) * pkin(8) + (-t236 - t3) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t65 ^ 2 + t227, t8, t214 - t27, t62, -t35 * t67 + t124 + t217, t12 * t79 + t35 * t65 - t118, -t31 * t65 - t120 + t217 + 0.2e1 * t226, pkin(5) * t26 - t27 * qJ(6) + (-t13 + t7) * t67 + (t6 - t193) * t65, 0.2e1 * t201 - t11 * t65 + t31 * t67 + (0.2e1 * qJD(6) - t12) * t79 + t118, t1 * qJ(6) - t2 * pkin(5) - t11 * t31 - t6 * t13 - g(1) * (-pkin(5) * t41 + qJ(6) * t42) - g(2) * (pkin(5) * t43 - qJ(6) * t44) + t193 * t7 + t138 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 + t215, t8, -t79 ^ 2 - t227, t120 - t223 - t226;];
tau_reg  = t9;
