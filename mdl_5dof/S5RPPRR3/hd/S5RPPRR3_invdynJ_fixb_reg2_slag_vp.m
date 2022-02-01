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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:40
% DurationCPUTime: 1.77s
% Computational Cost: add. (3153->267), mult. (6878->334), div. (0->0), fcn. (5142->16), ass. (0->151)
t141 = cos(pkin(8));
t119 = -pkin(1) * t141 - pkin(2);
t180 = qJDD(1) * t119;
t102 = qJDD(3) + t180;
t137 = qJ(1) + pkin(8);
t125 = sin(t137);
t127 = cos(t137);
t170 = -g(1) * t125 + g(2) * t127;
t215 = -t102 - t170;
t143 = sin(qJ(5));
t205 = cos(qJ(5));
t146 = cos(qJ(4));
t140 = cos(pkin(9));
t185 = qJD(1) * t140;
t172 = t146 * t185;
t138 = sin(pkin(9));
t144 = sin(qJ(4));
t188 = t138 * t144;
t174 = qJD(1) * t188;
t85 = -t172 + t174;
t96 = t138 * t146 + t140 * t144;
t87 = t96 * qJD(1);
t157 = t143 * t85 - t205 * t87;
t184 = qJD(4) * t144;
t173 = t138 * t184;
t178 = t140 * qJDD(1);
t179 = t138 * qJDD(1);
t175 = qJD(4) * t172 + t144 * t178 + t146 * t179;
t48 = qJD(1) * t173 - t175;
t162 = t144 * t179 - t146 * t178;
t90 = t96 * qJD(4);
t49 = qJD(1) * t90 + t162;
t153 = t157 * qJD(5) + t143 * t48 - t205 * t49;
t136 = qJD(4) + qJD(5);
t190 = t157 * t136;
t214 = t153 - t190;
t171 = qJD(5) * t205;
t182 = qJD(5) * t143;
t156 = -t143 * t49 - t85 * t171 - t87 * t182 - t205 * t48;
t40 = -t143 * t87 - t205 * t85;
t193 = t136 * t40;
t213 = t156 - t193;
t199 = t40 ^ 2;
t200 = t157 ^ 2;
t212 = -t199 + t200;
t198 = t40 * t157;
t167 = g(1) * t127 + g(2) * t125;
t139 = sin(pkin(8));
t113 = pkin(1) * t139 + qJ(3);
t104 = t113 * qJD(1);
t70 = t138 * qJD(2) + t140 * t104;
t64 = pkin(6) * t185 + t70;
t191 = t144 * t64;
t123 = t140 * qJD(2);
t63 = t123 + (-pkin(6) * qJD(1) - t104) * t138;
t29 = t146 * t63 - t191;
t23 = -pkin(7) * t87 + t29;
t22 = qJD(4) * pkin(4) + t23;
t30 = t144 * t63 + t146 * t64;
t24 = -pkin(7) * t85 + t30;
t121 = t140 * qJDD(2);
t93 = qJD(1) * qJD(3) + t113 * qJDD(1);
t59 = t121 + (-pkin(6) * qJDD(1) - t93) * t138;
t66 = t138 * qJDD(2) + t140 * t93;
t60 = pkin(6) * t178 + t66;
t168 = -t144 * t60 + t146 * t59;
t14 = -qJD(4) * t30 + t168;
t6 = qJDD(4) * pkin(4) + pkin(7) * t48 + t14;
t183 = qJD(4) * t146;
t177 = -t144 * t59 - t146 * t60 - t63 * t183;
t13 = -t64 * t184 - t177;
t7 = -pkin(7) * t49 + t13;
t1 = t143 * t6 + t22 * t171 - t24 * t182 + t205 * t7;
t135 = pkin(9) + qJ(4);
t128 = qJ(5) + t135;
t116 = sin(t128);
t117 = cos(t128);
t129 = t140 * pkin(3);
t209 = t119 - t129;
t82 = qJD(1) * t209 + qJD(3);
t47 = pkin(4) * t85 + t82;
t211 = g(3) * t116 + t167 * t117 - t47 * t40 - t1;
t176 = t205 * t24;
t9 = t143 * t22 + t176;
t2 = -t9 * qJD(5) - t143 * t7 + t205 * t6;
t210 = -g(3) * t117 + t167 * t116 + t47 * t157 + t2;
t189 = pkin(1) * qJDD(1);
t208 = t87 ^ 2;
t207 = pkin(4) * t49;
t206 = pkin(4) * t90;
t145 = sin(qJ(1));
t204 = pkin(1) * t145;
t197 = t87 * t85;
t142 = -pkin(6) - qJ(3);
t196 = pkin(6) + t113;
t118 = t129 + pkin(2);
t89 = -t140 * t183 + t173;
t187 = t146 * t140;
t95 = -t187 + t188;
t25 = t143 * t90 + t95 * t171 + t96 * t182 + t205 * t89;
t53 = -t143 * t95 + t205 * t96;
t195 = t153 * t53 - t25 * t40;
t194 = -t96 * t49 + t89 * t85;
t91 = t196 * t138;
t92 = t196 * t140;
t43 = -t144 * t91 + t146 * t92;
t192 = t143 * t24;
t133 = t138 ^ 2;
t134 = t140 ^ 2;
t186 = t133 + t134;
t42 = -t144 * t92 - t146 * t91;
t147 = cos(qJ(1));
t165 = g(1) * t145 - g(2) * t147;
t26 = t53 * qJD(5) - t143 * t89 + t205 * t90;
t52 = t143 * t96 + t205 * t95;
t164 = t156 * t52 - t157 * t26;
t163 = -t48 * t95 + t87 * t90;
t131 = qJDD(4) + qJDD(5);
t161 = t131 * t53 - t136 * t25;
t65 = -t138 * t93 + t121;
t160 = -t138 * t65 + t140 * t66;
t159 = t138 * (-t104 * t138 + t123) - t140 * t70;
t34 = -pkin(7) * t96 + t42;
t35 = -pkin(7) * t95 + t43;
t19 = -t143 * t35 + t205 * t34;
t20 = t143 * t34 + t205 * t35;
t155 = -t180 + t215;
t80 = qJDD(1) * t209 + qJDD(3);
t124 = sin(t135);
t126 = cos(t135);
t154 = -g(3) * t126 + t167 * t124;
t31 = qJD(3) * t187 - t91 * t183 + (-qJD(3) * t138 - qJD(4) * t92) * t144;
t152 = t170 + t80;
t32 = -t96 * qJD(3) - t43 * qJD(4);
t132 = pkin(7) - t142;
t130 = t147 * pkin(1);
t97 = pkin(4) * t126 + t118;
t83 = t85 ^ 2;
t62 = pkin(4) * t95 + t209;
t51 = -qJD(4) * t90 - qJDD(4) * t95;
t50 = -qJD(4) * t89 + qJDD(4) * t96;
t33 = t80 + t207;
t28 = pkin(7) * t89 + t32;
t27 = -pkin(7) * t90 + t31;
t15 = -t131 * t52 - t136 * t26;
t11 = t205 * t23 - t192;
t10 = -t143 * t23 - t176;
t8 = t205 * t22 - t192;
t4 = -t20 * qJD(5) - t143 * t27 + t205 * t28;
t3 = t19 * qJD(5) + t143 * t28 + t205 * t27;
t5 = [0, 0, 0, 0, 0, qJDD(1), t165, g(1) * t147 + g(2) * t145, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t141 * t189 - t170, -0.2e1 * t139 * t189 + t167, 0, (t165 + (t139 ^ 2 + t141 ^ 2) * t189) * pkin(1), t133 * qJDD(1), 0.2e1 * t138 * t178, 0, t134 * qJDD(1), 0, 0, t155 * t140, -t155 * t138, t93 * t186 + t160 - t167, t102 * t119 - g(1) * (-pkin(2) * t125 + qJ(3) * t127 - t204) - g(2) * (pkin(2) * t127 + qJ(3) * t125 + t130) + t160 * t113 - t159 * qJD(3), -t48 * t96 - t87 * t89, -t163 + t194, t50, t49 * t95 + t85 * t90, t51, 0, qJD(4) * t32 + qJDD(4) * t42 - t126 * t170 + t209 * t49 + t80 * t95 + t82 * t90, -qJD(4) * t31 - qJDD(4) * t43 + t124 * t170 - t209 * t48 + t80 * t96 - t82 * t89, -t13 * t95 - t14 * t96 + t29 * t89 - t30 * t90 - t31 * t85 - t32 * t87 + t42 * t48 - t43 * t49 - t167, t13 * t43 + t30 * t31 + t14 * t42 + t29 * t32 + t80 * t209 - g(1) * (-t118 * t125 - t127 * t142 - t204) - g(2) * (t118 * t127 - t125 * t142 + t130), t156 * t53 + t157 * t25, -t164 + t195, t161, -t153 * t52 - t26 * t40, t15, 0, -t117 * t170 + t131 * t19 + t136 * t4 - t153 * t62 - t206 * t40 + t26 * t47 + t33 * t52, t116 * t170 - t131 * t20 - t136 * t3 + t156 * t62 - t157 * t206 - t25 * t47 + t33 * t53, -t1 * t52 + t153 * t20 - t156 * t19 + t157 * t4 - t2 * t53 + t25 * t8 - t26 * t9 + t3 * t40 - t167, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t33 * t62 + t47 * t206 - g(1) * (-t125 * t97 + t127 * t132 - t204) - g(2) * (t125 * t132 + t127 * t97 + t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 * t66 + t140 * t65 - g(3), 0, 0, 0, 0, 0, 0, t51, -t50, t163 + t194, t13 * t96 - t14 * t95 - t29 * t90 - t30 * t89 - g(3), 0, 0, 0, 0, 0, 0, t15, -t161, t164 + t195, t1 * t53 - t2 * t52 - t25 * t9 - t26 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t179, -t186 * qJD(1) ^ 2, qJD(1) * t159 - t215, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t87 + t162, (-t85 - t174) * qJD(4) + t175, -t83 - t208, t29 * t87 + t30 * t85 + t152, 0, 0, 0, 0, 0, 0, -t153 - t190, t156 + t193, -t199 - t200, -t157 * t8 - t40 * t9 + t152 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t83 + t208, (t85 - t174) * qJD(4) + t175, -t197, -t162, qJDD(4), -t82 * t87 + t154 + t168, g(3) * t124 + t82 * t85 + t167 * t126 + (t29 + t191) * qJD(4) + t177, 0, 0, t198, t212, t213, -t198, t214, t131, -t10 * t136 + (t205 * t131 - t136 * t182 + t40 * t87) * pkin(4) + t210, t11 * t136 + (-t131 * t143 - t136 * t171 + t157 * t87) * pkin(4) + t211, -t10 * t157 - t11 * t40 - t9 * t157 + t8 * t40 + (-t205 * t156 + t143 * t153 + (-t143 * t157 + t205 * t40) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t205 * t2 + t1 * t143 - t47 * t87 + (-t143 * t8 + t205 * t9) * qJD(5) + t154) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t212, t213, -t198, t214, t131, t9 * t136 + t210, t8 * t136 + t211, 0, 0;];
tau_reg = t5;
