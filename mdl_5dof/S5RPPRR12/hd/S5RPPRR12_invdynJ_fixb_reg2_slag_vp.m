% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:26
% EndTime: 2019-12-31 18:07:30
% DurationCPUTime: 2.15s
% Computational Cost: add. (3279->313), mult. (6448->393), div. (0->0), fcn. (4422->10), ass. (0->170)
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t66 = t113 * t105 + t110 * t106;
t57 = t66 * qJD(1);
t211 = qJD(5) + t57;
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t108 = -pkin(1) - qJ(3);
t74 = t108 * qJD(1) + qJD(2);
t151 = -pkin(6) * qJD(1) + t74;
t48 = t151 * t105;
t49 = t151 * t106;
t26 = t110 * t49 + t113 * t48;
t23 = qJD(4) * pkin(7) + t26;
t172 = t113 * t106;
t154 = qJD(1) * t172;
t168 = qJD(1) * t105;
t155 = t110 * t168;
t59 = t154 - t155;
t87 = qJD(1) * qJ(2) + qJD(3);
t71 = pkin(3) * t168 + t87;
t24 = t57 * pkin(4) - t59 * pkin(7) + t71;
t8 = t109 * t24 + t112 * t23;
t213 = t211 * t8;
t7 = -t109 * t23 + t112 * t24;
t212 = t7 * t211;
t111 = sin(qJ(1));
t114 = cos(qJ(1));
t208 = g(1) * t111 - g(2) * t114;
t101 = pkin(8) + qJ(4);
t88 = sin(t101);
t89 = cos(t101);
t120 = g(3) * t88 - t208 * t89;
t202 = -qJD(1) * qJD(3) + qJDD(1) * t108;
t68 = qJDD(2) + t202;
t148 = -pkin(6) * qJDD(1) + t68;
t46 = t148 * t105;
t47 = t148 * t106;
t149 = t110 * t46 - t113 * t47;
t10 = -t26 * qJD(4) - t149;
t6 = -qJDD(4) * pkin(4) - t10;
t118 = -t6 + t120;
t210 = -pkin(7) * qJD(5) * t211 + t118;
t146 = t112 * t211;
t121 = -qJD(4) * t155 + t66 * qJDD(1);
t35 = qJD(4) * t154 + t121;
t31 = qJDD(5) + t35;
t183 = t109 * t31;
t209 = -t146 * t211 - t183;
t140 = g(1) * t114 + g(2) * t111;
t207 = t140 * t89;
t100 = t106 ^ 2;
t99 = t105 ^ 2;
t178 = t100 + t99;
t206 = t178 * t74;
t102 = qJDD(1) * qJ(2);
t103 = qJD(1) * qJD(2);
t204 = t102 + t103;
t72 = qJDD(3) + t204;
t122 = -t140 + t72;
t44 = t109 * qJD(4) + t112 * t59;
t176 = qJD(5) * t44;
t160 = t106 * qJDD(1);
t161 = t105 * qJDD(1);
t130 = -t110 * t161 + t113 * t160;
t61 = t66 * qJD(4);
t34 = qJD(1) * t61 - t130;
t19 = -t112 * qJDD(4) - t109 * t34 + t176;
t67 = -t110 * t105 + t172;
t83 = pkin(3) * t161;
t64 = t83 + t72;
t14 = t35 * pkin(4) + t34 * pkin(7) + t64;
t11 = t112 * t14;
t184 = qJD(5) * t8;
t166 = qJD(4) * t113;
t159 = -t110 * t47 - t113 * t46 - t49 * t166;
t167 = qJD(4) * t110;
t9 = -t48 * t167 - t159;
t5 = qJDD(4) * pkin(7) + t9;
t2 = -t109 * t5 + t11 - t184;
t203 = t2 + t213;
t198 = g(3) * t89;
t119 = -t208 * t88 - t198;
t201 = t59 ^ 2;
t200 = 0.2e1 * t103;
t62 = -t105 * t167 + t106 * t166;
t197 = t7 * t62;
t196 = t8 * t62;
t93 = t105 * pkin(3);
t163 = t112 * qJD(4);
t42 = t109 * t59 - t163;
t195 = t42 * t57;
t194 = t44 * t42;
t193 = t44 * t59;
t192 = t59 * t42;
t191 = t59 * t57;
t190 = -pkin(6) + t108;
t164 = qJD(5) * t112;
t189 = -t109 * t19 - t42 * t164;
t188 = -t61 * qJD(4) + t67 * qJDD(4);
t187 = t114 * pkin(1) + t111 * qJ(2);
t185 = qJD(5) * t7;
t182 = t109 * t44;
t181 = t110 * t48;
t27 = t112 * t31;
t165 = qJD(5) * t109;
t18 = -qJD(5) * t163 - t109 * qJDD(4) + t112 * t34 + t59 * t165;
t180 = t18 * t109;
t179 = t19 * t112;
t81 = qJ(2) + t93;
t177 = pkin(1) * qJDD(1);
t174 = t111 * t109;
t173 = t111 * t112;
t171 = t114 * t109;
t170 = t114 * t112;
t157 = t67 * t165;
t156 = t67 * t164;
t95 = t114 * qJ(2);
t153 = -t111 * pkin(1) + t95;
t152 = t178 * t68;
t147 = t109 * t211;
t145 = qJD(5) * t66 + qJD(1);
t144 = qJDD(2) - t177;
t143 = t88 * pkin(4) - t89 * pkin(7);
t25 = t113 * t49 - t181;
t22 = -qJD(4) * pkin(4) - t25;
t142 = -t22 * t61 + t6 * t67;
t1 = t109 * t14 + t112 * t5 + t185;
t141 = t1 - t212;
t138 = t109 * t8 + t112 * t7;
t137 = -t67 * t18 - t61 * t44;
t136 = -t18 * t66 + t44 * t62;
t135 = t67 * t19 - t61 * t42;
t134 = -t19 * t66 - t42 * t62;
t133 = -t211 * t61 + t31 * t67;
t132 = -t67 * t34 - t59 * t61;
t131 = t66 * t35 + t62 * t57;
t32 = t66 * pkin(4) - t67 * pkin(7) + t81;
t69 = t190 * t105;
t70 = t190 * t106;
t37 = t110 * t70 + t113 * t69;
t15 = -t109 * t37 + t112 * t32;
t16 = t109 * t32 + t112 * t37;
t36 = t110 * t69 - t113 * t70;
t107 = -pkin(6) - qJ(3);
t129 = t111 * t107 + t114 * t93 + t153;
t128 = -t114 * t107 + t111 * t93 + t187;
t127 = -t62 * qJD(4) - t66 * qJDD(4);
t126 = t27 + (-t109 * t57 - t165) * t211;
t125 = -qJD(5) * t24 + t198 - t5;
t124 = -pkin(7) * t31 + t211 * t22;
t117 = t10 * t67 - t25 * t61 + t26 * t62 + t9 * t66 - t208;
t116 = t204 + t122;
t115 = qJD(1) ^ 2;
t56 = t57 ^ 2;
t55 = t88 * t170 - t174;
t54 = t88 * t171 + t173;
t53 = t88 * t173 + t171;
t52 = -t88 * t174 + t170;
t33 = t59 * pkin(4) + t57 * pkin(7);
t29 = t62 * pkin(4) + t61 * pkin(7) + qJD(2);
t21 = t67 * qJD(3) + t37 * qJD(4);
t20 = -t66 * qJD(3) - t36 * qJD(4);
t13 = t109 * t33 + t112 * t25;
t12 = -t109 * t25 + t112 * t33;
t4 = -t16 * qJD(5) - t109 * t20 + t112 * t29;
t3 = t15 * qJD(5) + t109 * t29 + t112 * t20;
t17 = [0, 0, 0, 0, 0, qJDD(1), t208, t140, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t177 - t208, 0.2e1 * t102 + t200 - t140, -t144 * pkin(1) - g(1) * t153 - g(2) * t187 + (t102 + t200) * qJ(2), t100 * qJDD(1), -0.2e1 * t105 * t160, 0, t99 * qJDD(1), 0, 0, t116 * t105, t116 * t106, t208 + t178 * (-t202 - t68), t72 * qJ(2) + t87 * qJD(2) - g(1) * (t108 * t111 + t95) - g(2) * (t114 * qJ(3) + t187) + t108 * t152 - qJD(3) * t206, t132, t34 * t66 - t67 * t35 + t61 * t57 - t59 * t62, t188, t131, t127, 0, qJD(2) * t57 - t21 * qJD(4) - t36 * qJDD(4) - t140 * t88 + t81 * t35 + t71 * t62 + t64 * t66, qJD(2) * t59 - t20 * qJD(4) - t37 * qJDD(4) - t81 * t34 - t71 * t61 + t64 * t67 - t207, -t20 * t57 + t21 * t59 - t36 * t34 - t37 * t35 - t117, -g(1) * t129 - g(2) * t128 + t71 * qJD(2) - t10 * t36 + t26 * t20 - t25 * t21 + t9 * t37 + t64 * t81, t112 * t137 - t44 * t157, (t112 * t42 + t182) * t61 + (t180 - t179 + (t109 * t42 - t112 * t44) * qJD(5)) * t67, t112 * t133 - t157 * t211 + t136, t109 * t135 + t42 * t156, -t109 * t133 - t156 * t211 + t134, t211 * t62 + t31 * t66, -g(1) * t55 - g(2) * t53 + t109 * t142 + t15 * t31 + t156 * t22 + t36 * t19 + t2 * t66 + t21 * t42 + t211 * t4 + t197, g(1) * t54 - g(2) * t52 - t1 * t66 + t112 * t142 - t157 * t22 - t16 * t31 - t36 * t18 + t21 * t44 - t211 * t3 - t196, t15 * t18 - t16 * t19 - t3 * t42 - t4 * t44 + t207 + t138 * t61 + (-t1 * t109 - t2 * t112 + (t109 * t7 - t112 * t8) * qJD(5)) * t67, t1 * t16 + t8 * t3 + t2 * t15 + t7 * t4 + t6 * t36 + t22 * t21 - g(1) * (t114 * t143 + t129) - g(2) * (t111 * t143 + t128); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t115, -t115 * qJ(2) + t144 - t208, 0, 0, 0, 0, 0, 0, -t115 * t105, -t115 * t106, -t178 * qJDD(1), -t87 * qJD(1) + t152 - t208, 0, 0, 0, 0, 0, 0, -qJD(1) * t57 + t188, -qJD(1) * t59 + t127, -t131 - t132, -t71 * qJD(1) + t117, 0, 0, 0, 0, 0, 0, -t66 * t183 + (-t109 * t62 - t112 * t145) * t211 - t135, -t66 * t27 + (t145 * t109 - t112 * t62) * t211 - t137, (t145 * t44 + t134) * t112 + (t145 * t42 + t136) * t109, (-qJD(1) * t7 + t196 + (t1 - t185) * t66) * t112 + (-qJD(1) * t8 - t197 + (-t2 - t184) * t66) * t109 - t142 - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t160, -t178 * t115, qJD(1) * t206 + t122, 0, 0, 0, 0, 0, 0, (t59 + t154) * qJD(4) + t121, -0.2e1 * t57 * qJD(4) + t130, -t56 - t201, t25 * t59 + t26 * t57 + t122 + t83, 0, 0, 0, 0, 0, 0, t126 - t192, -t193 + t209, (t18 - t195) * t112 + t44 * t147 + t189, t141 * t109 + t203 * t112 - t22 * t59 - t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t56 + t201, t130, -t191, (t59 - t154) * qJD(4) - t121, qJDD(4), -t71 * t59 + t120 - t149, t71 * t57 + (t25 + t181) * qJD(4) + t159 - t119, 0, 0, t44 * t146 - t180, (-t18 - t195) * t112 - t211 * t182 + t189, -t193 - t209, t42 * t147 - t179, t126 + t192, -t211 * t59, -pkin(4) * t19 + t124 * t109 + t210 * t112 - t12 * t211 - t26 * t42 - t7 * t59, pkin(4) * t18 - t210 * t109 + t124 * t112 + t13 * t211 - t26 * t44 + t8 * t59, t12 * t44 + t13 * t42 + ((-t19 + t176) * pkin(7) + t141) * t112 + ((qJD(5) * t42 - t18) * pkin(7) - t203) * t109 + t119, -t7 * t12 - t8 * t13 - t22 * t26 + t118 * pkin(4) + (-qJD(5) * t138 + t1 * t112 - t2 * t109 + t119) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, -t42 ^ 2 + t44 ^ 2, t211 * t42 - t18, -t194, t211 * t44 - t19, t31, -g(1) * t52 - g(2) * t54 + t109 * t125 - t164 * t23 - t22 * t44 + t11 + t213, g(1) * t53 - g(2) * t55 + t22 * t42 + t212 + (qJD(5) * t23 - t14) * t109 + t125 * t112, 0, 0;];
tau_reg = t17;
