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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:23:01
% EndTime: 2020-01-03 11:23:08
% DurationCPUTime: 1.85s
% Computational Cost: add. (1409->284), mult. (3592->435), div. (0->0), fcn. (2888->10), ass. (0->162)
t123 = sin(pkin(9));
t126 = cos(pkin(9));
t128 = cos(pkin(7));
t125 = sin(pkin(7));
t127 = cos(pkin(8));
t200 = t125 * t127;
t71 = t123 * t200 + t126 * t128;
t63 = t71 * qJD(1);
t58 = qJD(5) + t63;
t206 = qJDD(1) * pkin(1);
t111 = qJDD(2) - t206;
t130 = sin(qJ(1));
t132 = cos(qJ(1));
t215 = -g(2) * t132 - g(3) * t130;
t218 = t215 - t111;
t129 = sin(qJ(5));
t131 = cos(qJ(5));
t198 = t125 * t132;
t124 = sin(pkin(8));
t195 = t128 * t132;
t78 = t124 * t130 + t127 * t195;
t45 = t123 * t198 + t126 * t78;
t193 = t130 * t127;
t77 = t124 * t195 - t193;
t217 = t129 * t45 - t131 * t77;
t216 = t129 * t77 + t131 * t45;
t214 = -t58 + qJD(5);
t180 = qJ(2) * qJDD(1);
t213 = g(2) * t130;
t212 = g(3) * t132;
t178 = qJD(1) * qJD(4);
t174 = t128 * qJDD(1);
t160 = t127 * t174;
t179 = qJD(1) * qJD(2);
t165 = t128 * t179;
t182 = qJD(3) * t125;
t82 = -pkin(2) * t128 - qJ(3) * t125 - pkin(1);
t52 = -qJD(1) * t182 + t82 * qJDD(1) + qJDD(2);
t25 = qJ(2) * t160 + t124 * t52 + t127 * t165;
t18 = (-qJ(4) * qJDD(1) - t178) * t128 + t25;
t148 = pkin(3) * t124 - qJ(4) * t127;
t175 = t125 * qJDD(1);
t79 = qJ(2) * t175 + t125 * t179 + qJDD(3);
t29 = (t148 * qJDD(1) - t127 * t178) * t125 + t79;
t4 = t123 * t29 + t126 * t18;
t184 = qJD(1) * t128;
t171 = qJ(2) * t184;
t70 = t82 * qJD(1) + qJD(2);
t40 = t124 * t70 + t127 * t171;
t32 = -qJ(4) * t184 + t40;
t186 = qJD(1) * t125;
t93 = qJ(2) * t186 + qJD(3);
t50 = t148 * t186 + t93;
t11 = t123 * t50 + t126 * t32;
t197 = t127 * t128;
t55 = qJ(2) * t197 + t124 * t82;
t47 = -qJ(4) * t128 + t55;
t59 = (qJ(2) + t148) * t125;
t20 = t123 * t59 + t126 * t47;
t61 = t71 * qJDD(1);
t57 = qJDD(5) + t61;
t211 = t129 * t57;
t209 = t131 * t57;
t207 = qJD(5) * t58;
t120 = t125 ^ 2;
t133 = qJD(1) ^ 2;
t205 = t120 * t133;
t204 = t123 * t124;
t203 = t124 * t125;
t202 = t124 * t129;
t201 = t124 * t131;
t199 = t125 * t130;
t196 = t128 * t130;
t194 = t128 * t133;
t192 = g(1) * t128 + g(3) * t198;
t191 = t132 * pkin(1) + t130 * qJ(2);
t119 = t124 ^ 2;
t189 = -t127 ^ 2 - t119;
t188 = t128 ^ 2 + t120;
t187 = qJD(1) * t124;
t185 = qJD(1) * t127;
t183 = qJD(2) * t128;
t181 = qJD(5) * t129;
t177 = qJDD(1) * t124;
t176 = qJDD(1) * t127;
t173 = t119 * t205;
t163 = t124 * t175;
t2 = pkin(6) * t163 + t4;
t161 = t124 * t174;
t24 = -qJ(2) * t161 - t124 * t165 + t127 * t52;
t21 = pkin(3) * t174 + qJDD(4) - t24;
t162 = t126 * t176;
t94 = t123 * t174;
t62 = t125 * t162 - t94;
t6 = pkin(4) * t61 - pkin(6) * t62 + t21;
t172 = -t129 * t2 + t131 * t6;
t170 = t124 * t186;
t169 = t129 * t187;
t168 = t131 * t187;
t167 = t123 * t184;
t166 = t126 * t186;
t164 = t119 * t175;
t39 = -t124 * t171 + t127 * t70;
t54 = -t124 * t128 * qJ(2) + t127 * t82;
t159 = t129 * t62 - t131 * t163;
t158 = t58 ^ 2;
t157 = t188 * t133;
t156 = pkin(2) * t195 + qJ(3) * t198 + t191;
t155 = 0.2e1 * t188;
t154 = t125 * t168;
t75 = t124 * t196 + t127 * t132;
t153 = g(2) * t77 + g(3) * t75;
t49 = t128 * pkin(3) - t54;
t74 = -t124 * t182 + t127 * t183;
t151 = t129 * t6 + t131 * t2;
t31 = pkin(3) * t184 + qJD(4) - t39;
t66 = t127 * t166 - t167;
t7 = pkin(4) * t63 - pkin(6) * t66 + t31;
t9 = pkin(6) * t170 + t11;
t150 = t129 * t9 - t131 * t7;
t149 = -t129 * t7 - t131 * t9;
t3 = -t123 * t18 + t126 * t29;
t10 = -t123 * t32 + t126 * t50;
t19 = -t123 * t47 + t126 * t59;
t72 = -t123 * t128 + t126 * t200;
t14 = pkin(4) * t71 - pkin(6) * t72 + t49;
t16 = pkin(6) * t203 + t20;
t147 = -t129 * t16 + t131 * t14;
t146 = t129 * t14 + t131 * t16;
t145 = (-qJD(4) * t127 + qJD(2)) * t125;
t115 = t130 * pkin(1);
t144 = pkin(2) * t196 - qJ(2) * t132 + qJ(3) * t199 + t115;
t143 = t125 * t201 - t129 * t72;
t42 = t125 * t202 + t131 * t72;
t142 = t126 * t201 - t127 * t129;
t141 = t126 * t202 + t127 * t131;
t73 = t124 * t183 + t127 * t182;
t140 = qJD(1) * t73 - qJDD(1) * t54 - t24;
t139 = qJD(1) * t74 + qJDD(1) * t55 + t25;
t12 = qJD(5) * t154 + t129 * t163 + t131 * t62 - t66 * t181;
t35 = t125 * t169 + t131 * t66;
t138 = t142 * t58;
t137 = t206 + t218;
t136 = t155 * t179 + t212;
t135 = t125 * t79 + (t179 + t180) * t120;
t76 = -t124 * t132 + t128 * t193;
t68 = (t123 * t125 + t126 * t197) * qJD(1);
t65 = t127 * t167 - t166;
t56 = -qJD(4) * t128 + t74;
t44 = t123 * t199 + t126 * t76;
t37 = t42 * qJD(5);
t36 = t143 * qJD(5);
t33 = t129 * t66 - t154;
t28 = t123 * t145 + t126 * t56;
t27 = t123 * t56 - t126 * t145;
t23 = t129 * t75 + t131 * t44;
t22 = -t129 * t44 + t131 * t75;
t15 = -pkin(4) * t203 - t19;
t13 = t35 * qJD(5) + t159;
t8 = -pkin(4) * t170 - t10;
t1 = -pkin(4) * t163 - t3;
t5 = [qJDD(1), t215, -t212 + t213, t137 * t128, -t137 * t125, t155 * t180 + t136 - t213, -t111 * pkin(1) - g(2) * t191 - g(3) * t115 + (t188 * t180 + t136) * qJ(2), -g(2) * t78 - g(3) * t76 + t135 * t124 + t140 * t128, t135 * t127 + t139 * t128 + t153, (-t139 * t124 + t140 * t127 + t215) * t125, t25 * t55 + t40 * t74 + t24 * t54 - t39 * t73 - g(2) * t156 - g(3) * t144 + (qJ(2) * t79 + qJD(2) * t93) * t125, -g(2) * t45 - g(3) * t44 + t21 * t71 + t49 * t61 + t63 * t73 + (-qJD(1) * t27 + qJDD(1) * t19 + t3) * t203, t73 * t66 + t49 * t62 + t21 * t72 - g(2) * (-t123 * t78 + t126 * t198) - g(3) * (-t123 * t76 + t126 * t199) + (-qJD(1) * t28 - qJDD(1) * t20 - t4) * t203, -t19 * t62 - t20 * t61 + t27 * t66 - t28 * t63 - t3 * t72 - t4 * t71 - t153, t4 * t20 + t11 * t28 + t3 * t19 - t10 * t27 + t21 * t49 + t31 * t73 - g(2) * (pkin(3) * t78 + qJ(4) * t77 + t156) - g(3) * (pkin(3) * t76 + qJ(4) * t75 + t144), t12 * t42 + t35 * t36, t12 * t143 - t13 * t42 - t33 * t36 - t35 * t37, t12 * t71 + t36 * t58 + t42 * t57, -t13 * t71 + t143 * t57 - t37 * t58, t57 * t71, (-t129 * t28 + t131 * t73) * t58 + t147 * t57 + t172 * t71 + t27 * t33 + t15 * t13 - t1 * t143 + t8 * t37 - g(2) * t216 - g(3) * t23 + (-t146 * t58 + t149 * t71) * qJD(5), -(t129 * t73 + t131 * t28) * t58 - t146 * t57 - t151 * t71 + t27 * t35 + t15 * t12 + t1 * t42 + t8 * t36 + g(2) * t217 - g(3) * t22 + (-t147 * t58 + t150 * t71) * qJD(5); 0, 0, 0, -t174, t175, -t157, -qJ(2) * t157 - t218, -t124 * t157 - t160, -t127 * t157 + t161, t189 * t175, t25 * t124 + t24 * t127 + (-t125 * t93 + (t124 * t39 - t127 * t40) * t128) * qJD(1) - t215, -t123 * t164 - t127 * t61 + (t125 * t65 - t128 * t63) * t187, -t126 * t164 - t127 * t62 + (t125 * t68 - t128 * t66) * t187, t68 * t63 - t65 * t66 + (t123 * t62 - t126 * t61) * t124, t10 * t65 - t11 * t68 - t127 * t21 + (-t123 * t3 + t126 * t4 - t31 * t184) * t124 - t215, 0, 0, 0, 0, 0, -t141 * t57 + t13 * t204 - (t128 * t168 - t129 * t68) * t58 - t65 * t33 - qJD(5) * t138, -t142 * t57 + t12 * t204 + (t128 * t169 + t131 * t68) * t58 - t65 * t35 + t141 * t207; 0, 0, 0, 0, 0, 0, 0, (-t127 * t194 + t177) * t125, (t124 * t194 + t176) * t125, t189 * t205, (-t213 + (t124 * t40 + t127 * t39) * qJD(1)) * t125 + t79 + t192, -t123 * t173 + (t126 * t177 - t63 * t185) * t125, -t126 * t173 + (-t123 * t177 - t66 * t185) * t125, -t123 * t61 - t126 * t62 + (t123 * t66 - t126 * t63) * t170, t123 * t4 + t126 * t3 + (-t213 + (-t127 * t31 + (-t10 * t123 + t11 * t126) * t124) * qJD(1)) * t125 + t192, 0, 0, 0, 0, 0, -t126 * t13 + (-t131 * t207 - t211) * t123 + (-t141 * t58 + t33 * t204) * t186, -t126 * t12 + (t181 * t58 - t209) * t123 + (t35 * t204 - t138) * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t170 + t61, -t94 + (-t63 * t187 + t162) * t125, -t63 ^ 2 - t66 ^ 2, -g(1) * t203 - g(2) * t75 + g(3) * t77 + t10 * t66 + t11 * t63 + t21, 0, 0, 0, 0, 0, -t129 * t158 - t66 * t33 + t209, -t131 * t158 - t66 * t35 - t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t33 * t58 + t12, -t214 * t35 - t159, t57, -g(1) * t143 - g(2) * t22 - g(3) * t217 + t214 * t149 - t8 * t35 + t172, g(1) * t42 + g(2) * t23 - g(3) * t216 + t214 * t150 + t8 * t33 - t151;];
tau_reg = t5;
