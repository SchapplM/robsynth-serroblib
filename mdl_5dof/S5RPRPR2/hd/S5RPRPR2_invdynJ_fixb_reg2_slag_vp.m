% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:12
% EndTime: 2020-01-03 11:34:17
% DurationCPUTime: 1.36s
% Computational Cost: add. (2291->210), mult. (3904->267), div. (0->0), fcn. (2533->16), ass. (0->141)
t137 = qJ(1) + pkin(8);
t127 = qJ(3) + t137;
t115 = sin(t127);
t116 = cos(t127);
t208 = g(2) * t116 + g(3) * t115;
t145 = sin(qJ(3));
t148 = cos(qJ(3));
t140 = sin(pkin(8));
t194 = pkin(1) * qJDD(1);
t178 = t140 * t194;
t142 = cos(pkin(8));
t118 = t142 * pkin(1) + pkin(2);
t94 = t118 * qJD(1);
t197 = t145 * t94;
t205 = pkin(1) * t140;
t179 = qJD(3) * t205;
t210 = qJD(1) * t179 - t118 * qJDD(1);
t171 = qJD(3) * t197 + t145 * t178 + t210 * t148;
t160 = qJDD(4) + t171;
t132 = qJDD(1) + qJDD(3);
t204 = t132 * pkin(3);
t31 = t160 - t204;
t211 = t31 + t208;
t180 = qJD(1) * t205;
t54 = -t145 * t180 + t148 * t94;
t165 = qJD(4) - t54;
t141 = cos(pkin(9));
t147 = cos(qJ(5));
t189 = t147 * t141;
t139 = sin(pkin(9));
t144 = sin(qJ(5));
t192 = t144 * t139;
t72 = -t189 + t192;
t73 = t147 * t139 + t144 * t141;
t136 = qJD(1) + qJD(3);
t60 = t73 * t136;
t133 = t139 ^ 2;
t134 = t141 ^ 2;
t184 = t133 + t134;
t209 = t136 * t184;
t207 = g(2) * t115 - g(3) * t116;
t206 = t60 ^ 2;
t203 = t141 * pkin(4);
t128 = t141 * pkin(7);
t177 = t136 * t192;
t58 = -t136 * t189 + t177;
t202 = t60 * t58;
t143 = -pkin(7) - qJ(4);
t85 = t143 * t139;
t86 = t141 * qJ(4) + t128;
t49 = -t144 * t86 + t147 * t85;
t201 = t49 * qJD(5) - t165 * t72;
t50 = t144 * t85 + t147 * t86;
t200 = -t50 * qJD(5) - t165 * t73;
t166 = t72 * t132;
t69 = t73 * qJD(5);
t39 = t136 * t69 + t166;
t175 = qJD(5) * t189;
t176 = qJD(5) * t192;
t68 = -t175 + t176;
t199 = -t73 * t39 + t68 * t58;
t117 = pkin(3) + t203;
t198 = t115 * t117 + t116 * t143;
t55 = t148 * t180 + t197;
t196 = t55 * t136;
t67 = t145 * t118 + t148 * t205;
t63 = t67 * qJD(3);
t195 = t63 * t136;
t183 = qJD(3) * t148;
t172 = -t210 * t145 + t148 * t178 + t94 * t183;
t26 = t132 * qJ(4) + t136 * qJD(4) + t172;
t21 = t139 * qJDD(2) + t141 * t26;
t48 = t136 * qJ(4) + t55;
t41 = t139 * qJD(2) + t141 * t48;
t193 = t141 * t132;
t188 = t116 * pkin(3) + t115 * qJ(4);
t124 = sin(t137);
t146 = sin(qJ(1));
t186 = t146 * pkin(1) + pkin(2) * t124;
t126 = cos(t137);
t149 = cos(qJ(1));
t185 = t149 * pkin(1) + pkin(2) * t126;
t182 = t211 * t139;
t181 = t73 * t132 + t136 * t175;
t120 = t141 * qJDD(2);
t15 = t120 + (-pkin(7) * t132 - t26) * t139;
t16 = pkin(7) * t193 + t21;
t174 = -t144 * t16 + t147 * t15;
t135 = pkin(9) + qJ(5);
t123 = sin(t135);
t22 = -t117 * t132 + t160;
t43 = -t117 * t136 + t165;
t173 = t208 * t123 + t22 * t73 - t43 * t68;
t170 = -t115 * t143 + t116 * t117;
t169 = t184 * t132;
t168 = t115 * pkin(3) - t116 * qJ(4);
t66 = t148 * t118 - t145 * t205;
t65 = -pkin(3) - t66;
t162 = -g(2) * t149 - g(3) * t146;
t38 = t136 * t176 - t181;
t161 = -t72 * t38 + t60 * t69;
t159 = -t196 - t204;
t122 = t141 * qJD(2);
t36 = t122 + (-pkin(7) * t136 - t48) * t139;
t37 = t136 * t128 + t41;
t10 = -t144 * t37 + t147 * t36;
t11 = t144 * t36 + t147 * t37;
t156 = t144 * t15 + t147 * t16;
t3 = t10 * qJD(5) + t156;
t4 = -t11 * qJD(5) + t174;
t158 = t10 * t68 - t11 * t69 - t3 * t72 - t4 * t73 - t207;
t157 = t132 * t65 + t195;
t64 = qJ(4) + t67;
t51 = (-pkin(7) - t64) * t139;
t52 = t141 * t64 + t128;
t23 = -t144 * t52 + t147 * t51;
t24 = t144 * t51 + t147 * t52;
t62 = t118 * t183 - t145 * t179;
t20 = -t139 * t26 + t120;
t155 = -t20 * t139 + t21 * t141 - t207;
t154 = -t172 + t207;
t153 = -t171 - t208;
t125 = cos(t135);
t151 = -t125 * t208 + t22 * t72 + t43 * t69;
t138 = qJDD(2) - g(1);
t109 = t134 * t132;
t108 = t133 * t132;
t80 = 0.2e1 * t139 * t193;
t57 = t58 ^ 2;
t56 = qJD(4) + t62;
t53 = t65 - t203;
t47 = -t136 * pkin(3) + t165;
t45 = -t69 * qJD(5) - t72 * qJDD(5);
t44 = -t68 * qJD(5) + t73 * qJDD(5);
t40 = -t139 * t48 + t122;
t13 = t39 * t72 + t58 * t69;
t12 = -t38 * t73 - t60 * t68;
t7 = -t24 * qJD(5) - t73 * t56;
t6 = t23 * qJD(5) - t72 * t56;
t5 = -t161 + t199;
t1 = [0, 0, 0, 0, 0, qJDD(1), t162, g(2) * t146 - g(3) * t149, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -g(2) * t126 - g(3) * t124 + 0.2e1 * t142 * t194, g(2) * t124 - g(3) * t126 - 0.2e1 * t178, 0, (t162 + (t140 ^ 2 + t142 ^ 2) * t194) * pkin(1), 0, 0, 0, 0, 0, t132, t66 * t132 + t153 - t195, -t67 * t132 - t62 * t136 + t154, 0, -g(2) * t185 - g(3) * t186 - t171 * t66 + t172 * t67 - t54 * t63 + t55 * t62, t108, t80, 0, t109, 0, 0, (-t211 - t157) * t141, t157 * t139 + t182, t64 * t169 + t56 * t209 + t155, t31 * t65 + t47 * t63 - g(2) * (t185 + t188) - g(3) * (t168 + t186) + (t21 * t64 + t41 * t56) * t141 + (-t20 * t64 - t40 * t56) * t139, t12, t5, t44, t13, t45, 0, t7 * qJD(5) + t23 * qJDD(5) + t53 * t39 + t63 * t58 + t151, -t6 * qJD(5) - t24 * qJDD(5) - t53 * t38 + t63 * t60 + t173, t23 * t38 - t24 * t39 - t6 * t58 - t7 * t60 + t158, t3 * t24 + t11 * t6 + t4 * t23 + t10 * t7 + t22 * t53 + t43 * t63 - g(2) * (t170 + t185) - g(3) * (t186 + t198); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t139 + t20 * t141 - g(1), 0, 0, 0, 0, 0, 0, t45, -t44, t161 + t199, -t10 * t69 - t11 * t68 + t3 * t73 - t4 * t72 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t153 + t196, t54 * t136 + t154, 0, 0, t108, t80, 0, t109, 0, 0, (-t211 - t159) * t141, t159 * t139 + t182, qJ(4) * t169 + t165 * t209 + t155, -t31 * pkin(3) - t47 * t55 - g(2) * t188 - g(3) * t168 + (t21 * qJ(4) + t165 * t41) * t141 + (-t20 * qJ(4) - t165 * t40) * t139, t12, t5, t44, t13, t45, 0, t200 * qJD(5) + t49 * qJDD(5) - t117 * t39 - t55 * t58 + t151, -t201 * qJD(5) - t50 * qJDD(5) + t117 * t38 - t55 * t60 + t173, -t200 * t60 - t201 * t58 + t49 * t38 - t50 * t39 + t158, -g(2) * t170 - g(3) * t198 + t200 * t10 + t201 * t11 - t22 * t117 + t3 * t50 + t4 * t49 - t43 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, t139 * t132, -t184 * t136 ^ 2, (t40 * t139 - t41 * t141) * t136 + t211, 0, 0, 0, 0, 0, 0, 0.2e1 * t60 * qJD(5) + t166, (-t58 - t177) * qJD(5) + t181, -t57 - t206, t10 * t60 + t11 * t58 + t208 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t57 + t206, (t58 - t177) * qJD(5) + t181, -t202, -t166, qJDD(5), -g(1) * t125 + t123 * t207 - t43 * t60 + t174, g(1) * t123 + t125 * t207 + t43 * t58 - t156, 0, 0;];
tau_reg = t1;
