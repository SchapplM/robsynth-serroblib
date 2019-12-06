% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:36
% EndTime: 2019-12-05 17:31:43
% DurationCPUTime: 1.86s
% Computational Cost: add. (1409->286), mult. (3592->437), div. (0->0), fcn. (2888->10), ass. (0->160)
t115 = sin(pkin(9));
t118 = cos(pkin(9));
t120 = cos(pkin(7));
t117 = sin(pkin(7));
t119 = cos(pkin(8));
t190 = t117 * t119;
t71 = t115 * t190 + t118 * t120;
t63 = t71 * qJD(1);
t58 = qJD(5) + t63;
t121 = sin(qJ(5));
t123 = cos(qJ(5));
t124 = cos(qJ(1));
t188 = t117 * t124;
t116 = sin(pkin(8));
t122 = sin(qJ(1));
t185 = t120 * t124;
t78 = t116 * t122 + t119 * t185;
t45 = t115 * t188 + t118 * t78;
t183 = t122 * t119;
t77 = t116 * t185 - t183;
t207 = t121 * t45 - t123 * t77;
t206 = t121 * t77 + t123 * t45;
t205 = -t58 + qJD(5);
t172 = qJ(2) * qJDD(1);
t204 = g(2) * t122;
t203 = g(2) * t124;
t202 = g(3) * t124;
t170 = qJD(1) * qJD(4);
t166 = t120 * qJDD(1);
t152 = t119 * t166;
t171 = qJD(1) * qJD(2);
t157 = t120 * t171;
t137 = pkin(2) * t120 + qJ(3) * t117 + pkin(1);
t174 = qJD(3) * t117;
t52 = -qJD(1) * t174 - t137 * qJDD(1) + qJDD(2);
t25 = qJ(2) * t152 + t116 * t52 + t119 * t157;
t18 = (-qJ(4) * qJDD(1) - t170) * t120 + t25;
t141 = pkin(3) * t116 - qJ(4) * t119;
t167 = t117 * qJDD(1);
t79 = qJ(2) * t167 + t117 * t171 + qJDD(3);
t29 = (t141 * qJDD(1) - t119 * t170) * t117 + t79;
t4 = t115 * t29 + t118 * t18;
t176 = qJD(1) * t120;
t163 = qJ(2) * t176;
t70 = -t137 * qJD(1) + qJD(2);
t40 = t116 * t70 + t119 * t163;
t32 = -qJ(4) * t176 + t40;
t178 = qJD(1) * t117;
t93 = qJ(2) * t178 + qJD(3);
t50 = t141 * t178 + t93;
t11 = t115 * t50 + t118 * t32;
t187 = t119 * t120;
t55 = qJ(2) * t187 - t116 * t137;
t47 = -qJ(4) * t120 + t55;
t59 = (qJ(2) + t141) * t117;
t20 = t115 * t59 + t118 * t47;
t61 = t71 * qJDD(1);
t57 = qJDD(5) + t61;
t201 = t121 * t57;
t199 = t123 * t57;
t197 = qJD(5) * t58;
t196 = qJDD(1) * pkin(1);
t112 = t117 ^ 2;
t125 = qJD(1) ^ 2;
t195 = t112 * t125;
t194 = t115 * t116;
t193 = t116 * t117;
t192 = t116 * t121;
t191 = t116 * t123;
t189 = t117 * t122;
t186 = t120 * t122;
t184 = t120 * t125;
t182 = g(1) * t120 + g(2) * t189;
t111 = t116 ^ 2;
t181 = -t119 ^ 2 - t111;
t180 = t120 ^ 2 + t112;
t179 = qJD(1) * t116;
t177 = qJD(1) * t119;
t175 = qJD(2) * t120;
t173 = qJD(5) * t121;
t169 = qJDD(1) * t116;
t168 = qJDD(1) * t119;
t165 = t111 * t195;
t155 = t116 * t167;
t2 = pkin(6) * t155 + t4;
t153 = t116 * t166;
t24 = -qJ(2) * t153 - t116 * t157 + t119 * t52;
t21 = pkin(3) * t166 + qJDD(4) - t24;
t154 = t118 * t168;
t94 = t115 * t166;
t62 = t117 * t154 - t94;
t6 = pkin(4) * t61 - pkin(6) * t62 + t21;
t164 = -t121 * t2 + t123 * t6;
t162 = t116 * t178;
t161 = t121 * t179;
t160 = t123 * t179;
t159 = t115 * t176;
t158 = t118 * t178;
t156 = t111 * t167;
t39 = -t116 * t163 + t119 * t70;
t54 = -t116 * t120 * qJ(2) - t119 * t137;
t151 = t121 * t62 - t123 * t155;
t150 = t58 ^ 2;
t149 = t180 * t125;
t148 = 0.2e1 * t180;
t147 = t117 * t160;
t75 = t116 * t186 + t119 * t124;
t146 = -g(2) * t77 - g(3) * t75;
t49 = t120 * pkin(3) - t54;
t145 = g(3) * t122 + t203;
t74 = -t116 * t174 + t119 * t175;
t144 = t121 * t6 + t123 * t2;
t31 = pkin(3) * t176 + qJD(4) - t39;
t66 = t119 * t158 - t159;
t7 = pkin(4) * t63 - pkin(6) * t66 + t31;
t9 = pkin(6) * t162 + t11;
t143 = t121 * t9 - t123 * t7;
t142 = -t121 * t7 - t123 * t9;
t3 = -t115 * t18 + t118 * t29;
t10 = -t115 * t32 + t118 * t50;
t19 = -t115 * t47 + t118 * t59;
t72 = -t115 * t120 + t118 * t190;
t14 = pkin(4) * t71 - pkin(6) * t72 + t49;
t16 = pkin(6) * t193 + t20;
t140 = -t121 * t16 + t123 * t14;
t139 = t121 * t14 + t123 * t16;
t138 = (-qJD(4) * t119 + qJD(2)) * t117;
t136 = t117 * t191 - t121 * t72;
t42 = t117 * t192 + t123 * t72;
t135 = t118 * t191 - t119 * t121;
t134 = t118 * t192 + t119 * t123;
t73 = t116 * t175 + t119 * t174;
t133 = qJD(1) * t73 - qJDD(1) * t54 - t24;
t132 = qJD(1) * t74 + qJDD(1) * t55 + t25;
t12 = qJD(5) * t147 + t121 * t155 + t123 * t62 - t66 * t173;
t35 = t117 * t161 + t123 * t66;
t131 = -t145 - t196;
t130 = t135 * t58;
t107 = qJDD(2) - t196;
t129 = -t107 - t131;
t128 = t148 * t171 + t204;
t127 = t117 * t79 + (t171 + t172) * t112;
t109 = t124 * qJ(2);
t76 = t116 * t124 - t120 * t183;
t68 = (t115 * t117 + t118 * t187) * qJD(1);
t65 = t119 * t159 - t158;
t56 = -qJD(4) * t120 + t74;
t44 = -t115 * t189 + t118 * t76;
t37 = t42 * qJD(5);
t36 = t136 * qJD(5);
t33 = t121 * t66 - t147;
t28 = t115 * t138 + t118 * t56;
t27 = t115 * t56 - t118 * t138;
t23 = -t121 * t75 + t123 * t44;
t22 = -t121 * t44 - t123 * t75;
t15 = -pkin(4) * t193 - t19;
t13 = t35 * qJD(5) + t151;
t8 = -pkin(4) * t162 - t10;
t1 = -pkin(4) * t155 - t3;
t5 = [qJDD(1), t145, t202 - t204, t129 * t120, -t129 * t117, t148 * t172 + t128 - t202, -g(3) * t109 + (-t107 + t145) * pkin(1) + (t180 * t172 + t128) * qJ(2), g(2) * t78 - g(3) * t76 + t127 * t116 + t133 * t120, t127 * t119 + t132 * t120 + t146, (-t132 * t116 + t133 * t119 + t145) * t117, t25 * t55 + t40 * t74 + t24 * t54 - t39 * t73 - g(2) * (-pkin(1) * t124 - pkin(2) * t185 - t122 * qJ(2)) - g(3) * (-pkin(1) * t122 - pkin(2) * t186 + t109) + (t79 * qJ(2) + t145 * qJ(3) + t93 * qJD(2)) * t117, g(2) * t45 - g(3) * t44 + t21 * t71 + t49 * t61 + t63 * t73 + (-qJD(1) * t27 + qJDD(1) * t19 + t3) * t193, t73 * t66 + t49 * t62 + t21 * t72 - g(2) * (t115 * t78 - t118 * t188) - g(3) * (-t115 * t76 - t118 * t189) + (-qJD(1) * t28 - qJDD(1) * t20 - t4) * t193, -t19 * t62 - t20 * t61 + t27 * t66 - t28 * t63 - t3 * t72 - t4 * t71 - t146, t4 * t20 + t11 * t28 + t3 * t19 - t10 * t27 + t21 * t49 + t31 * t73 - g(2) * (-pkin(3) * t78 - qJ(4) * t77) - g(3) * (pkin(3) * t76 - qJ(4) * t75 + t109) + t137 * t203 + (g(2) * qJ(2) + g(3) * t137) * t122, t12 * t42 + t35 * t36, t12 * t136 - t13 * t42 - t33 * t36 - t35 * t37, t12 * t71 + t36 * t58 + t42 * t57, -t13 * t71 + t136 * t57 - t37 * t58, t57 * t71, (-t121 * t28 + t123 * t73) * t58 + t140 * t57 + t164 * t71 + t27 * t33 + t15 * t13 - t1 * t136 + t8 * t37 + g(2) * t206 - g(3) * t23 + (-t139 * t58 + t142 * t71) * qJD(5), -(t121 * t73 + t123 * t28) * t58 - t139 * t57 - t144 * t71 + t27 * t35 + t15 * t12 + t1 * t42 + t8 * t36 - g(2) * t207 - g(3) * t22 + (-t140 * t58 + t143 * t71) * qJD(5); 0, 0, 0, -t166, t167, -t149, -qJ(2) * t149 + qJDD(2) + t131, -t116 * t149 - t152, -t119 * t149 + t153, t181 * t167, t116 * t25 + t119 * t24 + (-t117 * t93 + (t116 * t39 - t119 * t40) * t120) * qJD(1) - t145, -t115 * t156 - t119 * t61 + (t117 * t65 - t120 * t63) * t179, -t118 * t156 - t119 * t62 + (t117 * t68 - t120 * t66) * t179, t68 * t63 - t65 * t66 + (t115 * t62 - t118 * t61) * t116, t10 * t65 - t11 * t68 - t119 * t21 + (-t115 * t3 + t118 * t4 - t31 * t176) * t116 - t145, 0, 0, 0, 0, 0, -t134 * t57 + t13 * t194 - (t120 * t160 - t121 * t68) * t58 - t65 * t33 - qJD(5) * t130, -t135 * t57 + t12 * t194 + (t120 * t161 + t123 * t68) * t58 - t65 * t35 + t134 * t197; 0, 0, 0, 0, 0, 0, 0, (-t119 * t184 + t169) * t117, (t116 * t184 + t168) * t117, t181 * t195, (-t202 + (t116 * t40 + t119 * t39) * qJD(1)) * t117 + t79 + t182, -t115 * t165 + (t118 * t169 - t63 * t177) * t117, -t118 * t165 + (-t115 * t169 - t66 * t177) * t117, -t115 * t61 - t118 * t62 + (t115 * t66 - t118 * t63) * t162, t115 * t4 + t118 * t3 + (-t202 + (-t119 * t31 + (-t10 * t115 + t11 * t118) * t116) * qJD(1)) * t117 + t182, 0, 0, 0, 0, 0, -t118 * t13 + (-t123 * t197 - t201) * t115 + (-t134 * t58 + t33 * t194) * t178, -t118 * t12 + (t58 * t173 - t199) * t115 + (t35 * t194 - t130) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t162 + t61, -t94 + (-t63 * t179 + t154) * t117, -t63 ^ 2 - t66 ^ 2, -g(1) * t193 + g(2) * t75 - g(3) * t77 + t10 * t66 + t11 * t63 + t21, 0, 0, 0, 0, 0, -t121 * t150 - t66 * t33 + t199, -t123 * t150 - t66 * t35 - t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t33 * t58 + t12, -t205 * t35 - t151, t57, -g(1) * t136 - g(2) * t22 + g(3) * t207 + t205 * t142 - t8 * t35 + t164, g(1) * t42 + g(2) * t23 + g(3) * t206 + t205 * t143 + t8 * t33 - t144;];
tau_reg = t5;
