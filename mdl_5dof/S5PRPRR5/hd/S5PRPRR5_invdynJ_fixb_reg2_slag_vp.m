% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:48
% EndTime: 2019-12-05 15:54:54
% DurationCPUTime: 2.80s
% Computational Cost: add. (3033->304), mult. (7080->402), div. (0->0), fcn. (5514->14), ass. (0->156)
t138 = cos(qJ(2));
t184 = t138 * qJD(1);
t166 = qJD(3) - t184;
t128 = sin(pkin(9));
t205 = pkin(6) + qJ(3);
t102 = t205 * t128;
t130 = cos(pkin(9));
t103 = t205 * t130;
t134 = sin(qJ(4));
t137 = cos(qJ(4));
t193 = t137 * t130;
t194 = t134 * t128;
t97 = -t193 + t194;
t152 = t97 * t138;
t188 = qJD(4) * t137;
t204 = qJD(3) * t193 - t102 * t188 + (-qJD(3) * t128 - qJD(4) * t103) * t134 + qJD(1) * t152;
t60 = -t102 * t134 + t103 * t137;
t98 = t128 * t137 + t130 * t134;
t203 = -t60 * qJD(4) - t166 * t98;
t133 = sin(qJ(5));
t136 = cos(qJ(5));
t176 = qJD(2) * t193;
t177 = qJD(2) * t194;
t86 = -t176 + t177;
t88 = t98 * qJD(2);
t161 = t133 * t86 - t136 * t88;
t181 = t130 * qJDD(2);
t182 = t128 * qJDD(2);
t178 = qJD(4) * t176 + t134 * t181 + t137 * t182;
t49 = qJD(4) * t177 - t178;
t162 = t134 * t182 - t137 * t181;
t91 = t98 * qJD(4);
t50 = qJD(2) * t91 + t162;
t146 = qJD(5) * t161 + t133 * t49 - t136 * t50;
t127 = qJD(4) + qJD(5);
t199 = t161 * t127;
t226 = t146 - t199;
t186 = qJD(5) * t136;
t187 = qJD(5) * t133;
t153 = -t133 * t50 - t136 * t49 - t186 * t86 - t187 * t88;
t44 = -t133 * t88 - t136 * t86;
t198 = t44 * t127;
t225 = t153 - t198;
t208 = t161 ^ 2;
t209 = t44 ^ 2;
t224 = t208 - t209;
t223 = -t91 * pkin(7) + t204;
t189 = qJD(4) * t134;
t90 = t128 * t189 - t130 * t188;
t222 = t90 * pkin(7) + t203;
t207 = t44 * t161;
t121 = g(3) * t138;
t135 = sin(qJ(2));
t129 = sin(pkin(8));
t131 = cos(pkin(8));
t165 = g(1) * t131 + g(2) * t129;
t150 = t165 * t135 - t121;
t185 = t135 * qJD(1);
t106 = qJD(2) * qJ(3) + t185;
t170 = pkin(6) * qJD(2) + t106;
t78 = t170 * t130;
t201 = t134 * t78;
t77 = t170 * t128;
t33 = -t137 * t77 - t201;
t26 = -pkin(7) * t88 + t33;
t23 = qJD(4) * pkin(4) + t26;
t34 = -t134 * t77 + t137 * t78;
t27 = -pkin(7) * t86 + t34;
t79 = qJDD(2) * qJ(3) + t135 * qJDD(1) + (qJD(3) + t184) * qJD(2);
t169 = pkin(6) * qJDD(2) + t79;
t57 = t169 * t128;
t58 = t169 * t130;
t172 = -t134 * t58 - t137 * t57;
t17 = -t34 * qJD(4) + t172;
t6 = qJDD(4) * pkin(4) + t49 * pkin(7) + t17;
t179 = t134 * t57 - t137 * t58 + t188 * t77;
t16 = -t189 * t78 - t179;
t7 = -pkin(7) * t50 + t16;
t1 = (qJD(5) * t23 + t7) * t136 + t133 * t6 - t27 * t187;
t126 = pkin(9) + qJ(4);
t120 = qJ(5) + t126;
t114 = sin(t120);
t115 = cos(t120);
t195 = t131 * t138;
t196 = t129 * t138;
t210 = g(3) * t135;
t116 = pkin(3) * t130 + pkin(2);
t92 = -qJD(2) * t116 + t166;
t53 = t86 * pkin(4) + t92;
t221 = -t53 * t44 - g(1) * (-t114 * t129 - t115 * t195) - g(2) * (t114 * t131 - t115 * t196) + t115 * t210 - t1;
t200 = t136 * t27;
t11 = t133 * t23 + t200;
t2 = -qJD(5) * t11 - t133 * t7 + t136 * t6;
t220 = t53 * t161 - g(1) * (-t114 * t195 + t115 * t129) - g(2) * (-t114 * t196 - t115 * t131) + t114 * t210 + t2;
t81 = t97 * t135;
t218 = t165 * t138;
t183 = qJD(1) * qJD(2);
t197 = qJDD(2) * pkin(2);
t180 = t138 * qJDD(1);
t158 = t135 * t183 + qJDD(3) - t180;
t85 = t158 - t197;
t217 = t135 * (t165 + t183) + t197 - t85 - t121;
t216 = t88 ^ 2;
t215 = t50 * pkin(4);
t59 = -t102 * t137 - t103 * t134;
t37 = -pkin(7) * t98 + t59;
t38 = -pkin(7) * t97 + t60;
t18 = -t133 * t38 + t136 * t37;
t214 = qJD(5) * t18 + t133 * t222 + t136 * t223;
t19 = t133 * t37 + t136 * t38;
t213 = -qJD(5) * t19 - t133 * t223 + t136 * t222;
t206 = t88 * t86;
t202 = t133 * t27;
t192 = qJDD(1) - g(3);
t124 = t128 ^ 2;
t125 = t130 ^ 2;
t191 = t124 + t125;
t190 = qJD(2) * t135;
t171 = t191 * t79;
t168 = t191 * t138;
t167 = t191 * qJDD(2);
t164 = t116 * qJDD(2);
t163 = pkin(4) * t91 - t185;
t80 = t98 * t135;
t35 = t133 * t81 - t136 * t80;
t36 = -t133 * t80 - t136 * t81;
t52 = -t133 * t97 + t136 * t98;
t139 = qJD(2) ^ 2;
t157 = qJDD(2) * t138 - t135 * t139;
t148 = -t210 - t218;
t147 = t166 * t191;
t66 = -t164 + t158;
t145 = t158 - t150;
t144 = t171 + t148;
t141 = -t164 + t145;
t118 = sin(t126);
t119 = cos(t126);
t140 = t118 * t210 - g(1) * (-t118 * t195 + t119 * t129) - g(2) * (-t118 * t196 - t119 * t131);
t123 = -pkin(7) - t205;
t122 = qJDD(4) + qJDD(5);
t104 = -qJD(2) * pkin(2) + t166;
t99 = pkin(4) * t119 + t116;
t83 = t86 ^ 2;
t70 = pkin(4) * t97 - t116;
t51 = t133 * t98 + t136 * t97;
t40 = qJD(4) * t81 - t138 * t88;
t39 = -qJD(2) * t152 - qJD(4) * t80;
t28 = t66 + t215;
t21 = qJD(5) * t52 - t133 * t90 + t136 * t91;
t20 = t133 * t91 + t136 * t90 + t186 * t97 + t187 * t98;
t13 = t136 * t26 - t202;
t12 = -t133 * t26 - t200;
t10 = t136 * t23 - t202;
t9 = -qJD(5) * t36 - t133 * t39 + t136 * t40;
t8 = qJD(5) * t35 + t133 * t40 + t136 * t39;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, 0, t157, -qJDD(2) * t135 - t138 * t139, 0, -g(3) + (t135 ^ 2 + t138 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t157 * t130, -t157 * t128, t135 * t167 + t139 * t168, -t85 * t138 - g(3) + t135 * t171 + (t104 * t135 + t106 * t168) * qJD(2), 0, 0, 0, 0, 0, 0, qJD(4) * t40 - qJDD(4) * t80 - t138 * t50 + t190 * t86, -qJD(4) * t39 + qJDD(4) * t81 + t138 * t49 + t190 * t88, -t39 * t86 - t40 * t88 - t49 * t80 + t50 * t81, -t138 * t66 - t16 * t81 - t17 * t80 + t190 * t92 + t33 * t40 + t34 * t39 - g(3), 0, 0, 0, 0, 0, 0, t122 * t35 + t127 * t9 + t138 * t146 - t190 * t44, -t122 * t36 - t127 * t8 - t138 * t153 - t161 * t190, t146 * t36 - t153 * t35 + t161 * t9 + t44 * t8, t1 * t36 + t10 * t9 + t11 * t8 - t138 * t28 + t190 * t53 + t2 * t35 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t180 + t150, -t135 * t192 + t218, 0, 0, t124 * qJDD(2), 0.2e1 * t128 * t181, 0, t125 * qJDD(2), 0, 0, t217 * t130, -t217 * t128, qJ(3) * t167 + qJD(2) * t147 + t144, -t104 * t185 + (-t85 + t150) * pkin(2) + t144 * qJ(3) + t147 * t106, -t49 * t98 - t88 * t90, t49 * t97 - t50 * t98 + t86 * t90 - t88 * t91, -qJD(4) * t90 + qJDD(4) * t98, t50 * t97 + t86 * t91, -qJD(4) * t91 - qJDD(4) * t97, 0, qJD(4) * t203 + t59 * qJDD(4) - t116 * t50 + t119 * t150 - t185 * t86 + t66 * t97 + t92 * t91, -qJD(4) * t204 - t60 * qJDD(4) + t116 * t49 - t118 * t150 - t185 * t88 + t66 * t98 - t92 * t90, -t16 * t97 - t17 * t98 - t203 * t88 - t204 * t86 + t33 * t90 - t34 * t91 + t59 * t49 - t60 * t50 + t148, t16 * t60 + t17 * t59 - t66 * t116 - t92 * t185 - g(3) * (t116 * t138 + t135 * t205) + t204 * t34 + t203 * t33 + t165 * (t116 * t135 - t138 * t205), t153 * t52 + t161 * t20, t146 * t52 - t153 * t51 + t161 * t21 - t20 * t44, t122 * t52 - t127 * t20, -t146 * t51 - t21 * t44, -t122 * t51 - t127 * t21, 0, t115 * t150 + t18 * t122 + t127 * t213 - t146 * t70 - t163 * t44 + t53 * t21 + t28 * t51, -t114 * t150 - t19 * t122 - t127 * t214 + t153 * t70 - t161 * t163 - t53 * t20 + t28 * t52, -t1 * t51 + t10 * t20 - t11 * t21 + t146 * t19 - t153 * t18 + t161 * t213 - t2 * t52 + t214 * t44 + t148, t1 * t19 + t2 * t18 + t28 * t70 - g(3) * (-t123 * t135 + t138 * t99) + t163 * t53 + t214 * t11 + t213 * t10 + t165 * (t123 * t138 + t135 * t99); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t182, -t191 * t139, -qJD(2) * t106 * t191 + t145 - t197, 0, 0, 0, 0, 0, 0, 0.2e1 * t88 * qJD(4) + t162, (-t86 - t177) * qJD(4) + t178, -t83 - t216, t33 * t88 + t34 * t86 + t141, 0, 0, 0, 0, 0, 0, -t146 - t199, t153 + t198, -t208 - t209, -t10 * t161 - t11 * t44 + t141 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, -t83 + t216, (t86 - t177) * qJD(4) + t178, -t206, -t162, qJDD(4), -t92 * t88 + t140 + t172, t92 * t86 - g(1) * (-t118 * t129 - t119 * t195) - g(2) * (t118 * t131 - t119 * t196) + t119 * t210 + (t33 + t201) * qJD(4) + t179, 0, 0, t207, t224, t225, -t207, t226, t122, -t12 * t127 + (t122 * t136 - t127 * t187 + t44 * t88) * pkin(4) + t220, t13 * t127 + (-t122 * t133 - t127 * t186 + t161 * t88) * pkin(4) + t221, t10 * t44 - t11 * t161 - t12 * t161 - t13 * t44 + (t133 * t146 - t136 * t153 + (-t133 * t161 + t136 * t44) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t133 + t2 * t136 - t53 * t88 + (-t10 * t133 + t11 * t136) * qJD(5) + t140) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t224, t225, -t207, t226, t122, t11 * t127 + t220, t10 * t127 + t221, 0, 0;];
tau_reg = t3;
