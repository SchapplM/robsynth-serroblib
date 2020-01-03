% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:29
% EndTime: 2020-01-03 11:28:36
% DurationCPUTime: 1.86s
% Computational Cost: add. (3153->267), mult. (6878->334), div. (0->0), fcn. (5142->16), ass. (0->148)
t144 = sin(qJ(5));
t202 = cos(qJ(5));
t147 = cos(qJ(4));
t141 = cos(pkin(9));
t184 = qJD(1) * t141;
t172 = t147 * t184;
t139 = sin(pkin(9));
t145 = sin(qJ(4));
t188 = t145 * t139;
t173 = qJD(1) * t188;
t85 = -t172 + t173;
t96 = t147 * t139 + t145 * t141;
t87 = t96 * qJD(1);
t157 = t144 * t85 - t202 * t87;
t177 = t141 * qJDD(1);
t178 = t139 * qJDD(1);
t174 = qJD(4) * t172 + t145 * t177 + t147 * t178;
t48 = qJD(4) * t173 - t174;
t162 = t145 * t178 - t147 * t177;
t90 = t96 * qJD(4);
t49 = qJD(1) * t90 + t162;
t153 = qJD(5) * t157 + t144 * t48 - t202 * t49;
t137 = qJD(4) + qJD(5);
t191 = t157 * t137;
t213 = t153 - t191;
t171 = qJD(5) * t202;
t181 = qJD(5) * t144;
t156 = -t144 * t49 - t85 * t171 - t181 * t87 - t202 * t48;
t40 = -t144 * t87 - t202 * t85;
t190 = t40 * t137;
t212 = t156 - t190;
t199 = t157 ^ 2;
t200 = t40 ^ 2;
t211 = t199 - t200;
t198 = t40 * t157;
t142 = cos(pkin(8));
t119 = -t142 * pkin(1) - pkin(2);
t179 = qJDD(1) * t119;
t100 = qJDD(3) + t179;
t138 = qJ(1) + pkin(8);
t125 = sin(t138);
t127 = cos(t138);
t207 = g(2) * t127 + g(3) * t125;
t210 = t100 + t207;
t206 = g(2) * t125 - g(3) * t127;
t140 = sin(pkin(8));
t111 = t140 * pkin(1) + qJ(3);
t102 = t111 * qJD(1);
t70 = t139 * qJD(2) + t141 * t102;
t64 = pkin(6) * t184 + t70;
t192 = t145 * t64;
t123 = t141 * qJD(2);
t63 = t123 + (-pkin(6) * qJD(1) - t102) * t139;
t29 = t147 * t63 - t192;
t23 = -t87 * pkin(7) + t29;
t22 = qJD(4) * pkin(4) + t23;
t30 = t145 * t63 + t147 * t64;
t24 = -t85 * pkin(7) + t30;
t121 = t141 * qJDD(2);
t93 = qJD(1) * qJD(3) + qJDD(1) * t111;
t59 = t121 + (-pkin(6) * qJDD(1) - t93) * t139;
t66 = t139 * qJDD(2) + t141 * t93;
t60 = pkin(6) * t177 + t66;
t168 = -t145 * t60 + t147 * t59;
t14 = -t30 * qJD(4) + t168;
t6 = qJDD(4) * pkin(4) + t48 * pkin(7) + t14;
t182 = qJD(4) * t147;
t176 = -t145 * t59 - t147 * t60 - t63 * t182;
t183 = qJD(4) * t145;
t13 = -t183 * t64 - t176;
t7 = -t49 * pkin(7) + t13;
t1 = t144 * t6 + t22 * t171 - t24 * t181 + t202 * t7;
t136 = pkin(9) + qJ(4);
t128 = qJ(5) + t136;
t116 = sin(t128);
t117 = cos(t128);
t129 = t141 * pkin(3);
t205 = t119 - t129;
t82 = qJD(1) * t205 + qJD(3);
t47 = t85 * pkin(4) + t82;
t209 = g(1) * t116 + t206 * t117 - t47 * t40 - t1;
t175 = t202 * t24;
t9 = t144 * t22 + t175;
t2 = -qJD(5) * t9 - t144 * t7 + t202 * t6;
t208 = -g(1) * t117 + t206 * t116 + t47 * t157 + t2;
t189 = pkin(1) * qJDD(1);
t204 = t87 ^ 2;
t203 = pkin(4) * t90;
t197 = t87 * t85;
t143 = -pkin(6) - qJ(3);
t196 = pkin(6) + t111;
t118 = t129 + pkin(2);
t89 = t139 * t183 - t141 * t182;
t187 = t147 * t141;
t95 = -t187 + t188;
t25 = t144 * t90 + t95 * t171 + t96 * t181 + t202 * t89;
t53 = -t144 * t95 + t202 * t96;
t195 = t153 * t53 - t25 * t40;
t194 = -t96 * t49 + t89 * t85;
t91 = t196 * t139;
t92 = t196 * t141;
t43 = -t145 * t91 + t147 * t92;
t193 = t144 * t24;
t134 = t139 ^ 2;
t135 = t141 ^ 2;
t185 = t134 + t135;
t42 = -t145 * t92 - t147 * t91;
t146 = sin(qJ(1));
t148 = cos(qJ(1));
t165 = -g(2) * t148 - g(3) * t146;
t26 = t53 * qJD(5) - t144 * t89 + t202 * t90;
t52 = t144 * t96 + t202 * t95;
t164 = t156 * t52 - t157 * t26;
t163 = -t95 * t48 + t87 * t90;
t132 = qJDD(4) + qJDD(5);
t161 = t53 * t132 - t25 * t137;
t65 = -t139 * t93 + t121;
t160 = -t65 * t139 + t66 * t141;
t159 = (-t139 * t102 + t123) * t139 - t70 * t141;
t34 = -t96 * pkin(7) + t42;
t35 = -t95 * pkin(7) + t43;
t19 = -t144 * t35 + t202 * t34;
t20 = t144 * t34 + t202 * t35;
t155 = t179 + t210;
t80 = qJDD(1) * t205 + qJDD(3);
t124 = sin(t136);
t126 = cos(t136);
t154 = -g(1) * t126 + t124 * t206;
t31 = qJD(3) * t187 - t91 * t182 + (-qJD(3) * t139 - qJD(4) * t92) * t145;
t33 = t49 * pkin(4) + t80;
t32 = -qJD(3) * t96 - qJD(4) * t43;
t133 = -pkin(7) + t143;
t131 = t148 * pkin(1);
t130 = t146 * pkin(1);
t97 = pkin(4) * t126 + t118;
t83 = t85 ^ 2;
t62 = t95 * pkin(4) + t205;
t51 = -t90 * qJD(4) - t95 * qJDD(4);
t50 = -t89 * qJD(4) + t96 * qJDD(4);
t28 = t89 * pkin(7) + t32;
t27 = -t90 * pkin(7) + t31;
t15 = -t52 * t132 - t26 * t137;
t11 = t202 * t23 - t193;
t10 = -t144 * t23 - t175;
t8 = t202 * t22 - t193;
t4 = -t20 * qJD(5) - t144 * t27 + t202 * t28;
t3 = t19 * qJD(5) + t144 * t28 + t202 * t27;
t5 = [0, 0, 0, 0, 0, qJDD(1), t165, g(2) * t146 - g(3) * t148, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t142 * t189 - t207, -0.2e1 * t140 * t189 + t206, 0, (t165 + (t140 ^ 2 + t142 ^ 2) * t189) * pkin(1), t134 * qJDD(1), 0.2e1 * t139 * t177, 0, t135 * qJDD(1), 0, 0, -t155 * t141, t155 * t139, t93 * t185 + t160 - t206, t100 * t119 - g(2) * (t127 * pkin(2) + t125 * qJ(3) + t131) - g(3) * (t125 * pkin(2) - t127 * qJ(3) + t130) + t160 * t111 - t159 * qJD(3), -t48 * t96 - t87 * t89, -t163 + t194, t50, t49 * t95 + t85 * t90, t51, 0, t32 * qJD(4) + t42 * qJDD(4) - t126 * t207 + t205 * t49 + t80 * t95 + t82 * t90, -t31 * qJD(4) - t43 * qJDD(4) + t124 * t207 - t205 * t48 + t80 * t96 - t82 * t89, -t13 * t95 - t14 * t96 + t29 * t89 - t30 * t90 - t31 * t85 - t32 * t87 + t42 * t48 - t43 * t49 - t206, t13 * t43 + t30 * t31 + t14 * t42 + t29 * t32 + t80 * t205 - g(2) * (t127 * t118 - t125 * t143 + t131) - g(3) * (t125 * t118 + t127 * t143 + t130), t156 * t53 + t157 * t25, -t164 + t195, t161, -t153 * t52 - t26 * t40, t15, 0, -t117 * t207 + t19 * t132 + t4 * t137 - t153 * t62 - t203 * t40 + t47 * t26 + t33 * t52, t116 * t207 - t20 * t132 - t3 * t137 + t156 * t62 - t157 * t203 - t47 * t25 + t33 * t53, -t1 * t52 + t153 * t20 - t156 * t19 + t157 * t4 - t2 * t53 + t8 * t25 - t9 * t26 + t3 * t40 - t206, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t33 * t62 + t47 * t203 - g(2) * (-t125 * t133 + t127 * t97 + t131) - g(3) * (t125 * t97 + t127 * t133 + t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t139 + t65 * t141 - g(1), 0, 0, 0, 0, 0, 0, t51, -t50, t163 + t194, t13 * t96 - t14 * t95 - t29 * t90 - t30 * t89 - g(1), 0, 0, 0, 0, 0, 0, t15, -t161, t164 + t195, t1 * t53 - t2 * t52 - t9 * t25 - t8 * t26 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, t178, -t185 * qJD(1) ^ 2, qJD(1) * t159 + t210, 0, 0, 0, 0, 0, 0, 0.2e1 * t87 * qJD(4) + t162, (-t85 - t173) * qJD(4) + t174, -t83 - t204, t29 * t87 + t30 * t85 + t207 + t80, 0, 0, 0, 0, 0, 0, -t153 - t191, t156 + t190, -t199 - t200, -t157 * t8 - t9 * t40 + t207 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t83 + t204, (t85 - t173) * qJD(4) + t174, -t197, -t162, qJDD(4), -t82 * t87 + t154 + t168, g(1) * t124 + t82 * t85 + t206 * t126 + (t29 + t192) * qJD(4) + t176, 0, 0, t198, t211, t212, -t198, t213, t132, -t10 * t137 + (t202 * t132 - t137 * t181 + t40 * t87) * pkin(4) + t208, t11 * t137 + (-t132 * t144 - t137 * t171 + t157 * t87) * pkin(4) + t209, -t10 * t157 - t11 * t40 - t9 * t157 + t8 * t40 + (-t202 * t156 + t144 * t153 + (-t144 * t157 + t202 * t40) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t202 * t2 + t1 * t144 - t47 * t87 + (-t144 * t8 + t202 * t9) * qJD(5) + t154) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t211, t212, -t198, t213, t132, t9 * t137 + t208, t8 * t137 + t209, 0, 0;];
tau_reg = t5;
