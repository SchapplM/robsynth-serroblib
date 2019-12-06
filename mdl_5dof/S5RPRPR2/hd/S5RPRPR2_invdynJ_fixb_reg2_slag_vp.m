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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:49:44
% EndTime: 2019-12-05 17:49:48
% DurationCPUTime: 1.28s
% Computational Cost: add. (2291->214), mult. (3904->268), div. (0->0), fcn. (2533->16), ass. (0->141)
t134 = cos(pkin(8));
t112 = t134 * pkin(1) + pkin(2);
t132 = sin(pkin(8));
t196 = pkin(1) * t132;
t174 = qJD(3) * t196;
t200 = qJD(1) * t174 - t112 * qJDD(1);
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t175 = qJD(1) * t196;
t92 = t112 * qJD(1);
t54 = -t137 * t175 + t140 * t92;
t160 = qJD(4) - t54;
t129 = qJ(1) + pkin(8);
t121 = qJ(3) + t129;
t109 = sin(t121);
t110 = cos(t121);
t199 = g(2) * t110 + g(3) * t109;
t133 = cos(pkin(9));
t139 = cos(qJ(5));
t180 = t139 * t133;
t131 = sin(pkin(9));
t136 = sin(qJ(5));
t183 = t136 * t131;
t72 = -t180 + t183;
t73 = t139 * t131 + t136 * t133;
t128 = qJD(1) + qJD(3);
t60 = t73 * t128;
t185 = pkin(1) * qJDD(1);
t125 = t131 ^ 2;
t126 = t133 ^ 2;
t178 = t125 + t126;
t198 = t128 * t178;
t179 = g(2) * t109 - g(3) * t110;
t197 = t60 ^ 2;
t124 = qJDD(1) + qJDD(3);
t195 = t124 * pkin(3);
t194 = t133 * pkin(4);
t122 = t133 * pkin(7);
t172 = t128 * t183;
t58 = -t128 * t180 + t172;
t193 = t60 * t58;
t135 = -pkin(7) - qJ(4);
t83 = t135 * t131;
t84 = t133 * qJ(4) + t122;
t49 = -t136 * t84 + t139 * t83;
t192 = t49 * qJD(5) - t160 * t72;
t50 = t136 * t83 + t139 * t84;
t191 = -t50 * qJD(5) - t160 * t73;
t161 = t72 * t124;
t69 = t73 * qJD(5);
t39 = t128 * t69 + t161;
t170 = qJD(5) * t180;
t171 = qJD(5) * t183;
t68 = -t170 + t171;
t190 = -t73 * t39 + t68 * t58;
t189 = t199 * t133;
t188 = t137 * t92;
t55 = t140 * t175 + t188;
t187 = t55 * t128;
t67 = t137 * t112 + t140 * t196;
t63 = t67 * qJD(3);
t186 = t63 * t128;
t173 = t132 * t185;
t177 = qJD(3) * t140;
t166 = -t200 * t137 + t140 * t173 + t92 * t177;
t26 = t124 * qJ(4) + t128 * qJD(4) + t166;
t21 = t131 * qJDD(2) + t133 * t26;
t48 = t128 * qJ(4) + t55;
t41 = t131 * qJD(2) + t133 * t48;
t184 = t133 * t124;
t176 = t73 * t124 + t128 * t170;
t111 = pkin(3) + t194;
t169 = -t109 * pkin(3) + t110 * qJ(4);
t114 = t133 * qJDD(2);
t15 = t114 + (-pkin(7) * t124 - t26) * t131;
t16 = pkin(7) * t184 + t21;
t168 = -t136 * t16 + t139 * t15;
t127 = pkin(9) + qJ(5);
t119 = cos(t127);
t165 = qJD(3) * t188 + t137 * t173 + t200 * t140;
t153 = qJDD(4) + t165;
t22 = -t111 * t124 + t153;
t43 = -t111 * t128 + t160;
t167 = t199 * t119 + t22 * t72 + t43 * t69;
t164 = t109 * t135 - t110 * t111;
t163 = t178 * t124;
t66 = t140 * t112 - t137 * t196;
t65 = -pkin(3) - t66;
t118 = sin(t129);
t138 = sin(qJ(1));
t159 = -t138 * pkin(1) - pkin(2) * t118;
t120 = cos(t129);
t141 = cos(qJ(1));
t158 = -t141 * pkin(1) - pkin(2) * t120;
t155 = g(2) * t141 + g(3) * t138;
t38 = t128 * t171 - t176;
t154 = -t72 * t38 + t60 * t69;
t152 = t187 + t195;
t116 = t133 * qJD(2);
t36 = t116 + (-pkin(7) * t128 - t48) * t131;
t37 = t128 * t122 + t41;
t10 = -t136 * t37 + t139 * t36;
t11 = t136 * t36 + t139 * t37;
t148 = t136 * t15 + t139 * t16;
t3 = t10 * qJD(5) + t148;
t4 = -t11 * qJD(5) + t168;
t151 = t10 * t68 - t11 * t69 - t3 * t72 - t4 * t73 + t179;
t150 = -t110 * pkin(3) - t109 * qJ(4);
t149 = -t124 * t65 - t186;
t64 = qJ(4) + t67;
t51 = (-pkin(7) - t64) * t131;
t52 = t133 * t64 + t122;
t23 = -t136 * t52 + t139 * t51;
t24 = t136 * t51 + t139 * t52;
t147 = -t109 * t111 - t110 * t135;
t62 = t112 * t177 - t137 * t174;
t20 = -t131 * t26 + t114;
t146 = -t20 * t131 + t21 * t133 + t179;
t145 = -t166 - t179;
t144 = t165 - t199;
t31 = t153 - t195;
t117 = sin(t127);
t143 = -t117 * t199 + t22 * t73 - t43 * t68;
t130 = qJDD(2) - g(1);
t105 = t126 * t124;
t104 = t125 * t124;
t78 = 0.2e1 * t131 * t184;
t57 = t58 ^ 2;
t56 = qJD(4) + t62;
t53 = t65 - t194;
t47 = -t128 * pkin(3) + t160;
t45 = -t69 * qJD(5) - t72 * qJDD(5);
t44 = -t68 * qJD(5) + t73 * qJDD(5);
t40 = -t131 * t48 + t116;
t29 = t31 * t131;
t13 = t39 * t72 + t58 * t69;
t12 = -t38 * t73 - t60 * t68;
t7 = -t24 * qJD(5) - t73 * t56;
t6 = t23 * qJD(5) - t72 * t56;
t5 = -t154 + t190;
t1 = [0, 0, 0, 0, 0, qJDD(1), t155, -g(2) * t138 + g(3) * t141, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(2) * t120 + g(3) * t118 + 0.2e1 * t134 * t185, -g(2) * t118 + g(3) * t120 - 0.2e1 * t173, 0, (t155 + (t132 ^ 2 + t134 ^ 2) * t185) * pkin(1), 0, 0, 0, 0, 0, t124, t66 * t124 - t144 - t186, -t67 * t124 - t62 * t128 + t145, 0, -g(2) * t158 - g(3) * t159 - t165 * t66 + t166 * t67 - t54 * t63 + t55 * t62, t104, t78, 0, t105, 0, 0, (t149 - t31) * t133 + t189, t29 + (-t149 - t199) * t131, t64 * t163 + t56 * t198 + t146, t31 * t65 + t47 * t63 - g(2) * (t150 + t158) - g(3) * (t159 + t169) + (t21 * t64 + t41 * t56) * t133 + (-t20 * t64 - t40 * t56) * t131, t12, t5, t44, t13, t45, 0, t7 * qJD(5) + t23 * qJDD(5) + t53 * t39 + t63 * t58 + t167, -t6 * qJD(5) - t24 * qJDD(5) - t53 * t38 + t63 * t60 + t143, t23 * t38 - t24 * t39 - t6 * t58 - t7 * t60 + t151, t3 * t24 + t11 * t6 + t4 * t23 + t10 * t7 + t22 * t53 + t43 * t63 - g(2) * (t158 + t164) - g(3) * (t147 + t159); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t131 + t20 * t133 - g(1), 0, 0, 0, 0, 0, 0, t45, -t44, t154 + t190, -t10 * t69 - t11 * t68 + t3 * t73 - t4 * t72 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t144 + t187, t54 * t128 + t145, 0, 0, t104, t78, 0, t105, 0, 0, (t152 - t31) * t133 + t189, t29 + (-t152 - t199) * t131, qJ(4) * t163 + t160 * t198 + t146, -t31 * pkin(3) - t47 * t55 - g(2) * t150 - g(3) * t169 + (t21 * qJ(4) + t160 * t41) * t133 + (-t20 * qJ(4) - t160 * t40) * t131, t12, t5, t44, t13, t45, 0, t191 * qJD(5) + t49 * qJDD(5) - t111 * t39 - t55 * t58 + t167, -t192 * qJD(5) - t50 * qJDD(5) + t111 * t38 - t55 * t60 + t143, -t191 * t60 - t192 * t58 + t49 * t38 - t50 * t39 + t151, -g(2) * t164 - g(3) * t147 + t191 * t10 + t192 * t11 - t22 * t111 + t3 * t50 + t4 * t49 - t43 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t131 * t124, -t178 * t128 ^ 2, (t40 * t131 - t41 * t133) * t128 + t31 - t199, 0, 0, 0, 0, 0, 0, 0.2e1 * t60 * qJD(5) + t161, (-t58 - t172) * qJD(5) + t176, -t57 - t197, t10 * t60 + t11 * t58 - t199 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t57 + t197, (t58 - t172) * qJD(5) + t176, -t193, -t161, qJDD(5), -g(1) * t119 - t117 * t179 - t43 * t60 + t168, g(1) * t117 - t119 * t179 + t43 * t58 - t148, 0, 0;];
tau_reg = t1;
