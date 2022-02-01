% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:52
% EndTime: 2022-01-20 11:30:56
% DurationCPUTime: 1.46s
% Computational Cost: add. (3105->246), mult. (5261->303), div. (0->0), fcn. (3029->16), ass. (0->162)
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t112 = sin(qJ(2));
t161 = qJDD(1) * t112;
t116 = cos(qJ(2));
t165 = qJD(2) * t116;
t130 = (qJD(1) * t165 + t161) * pkin(1);
t103 = qJDD(1) + qJDD(2);
t181 = pkin(1) * qJD(1);
t158 = t112 * t181;
t193 = pkin(1) * t116;
t94 = qJDD(1) * t193;
t49 = pkin(2) * t103 - qJD(2) * t158 + t94;
t104 = qJD(1) + qJD(2);
t64 = pkin(2) * t104 + t116 * t181;
t146 = qJD(3) * t158;
t67 = t111 * t146;
t21 = t111 * t49 + (qJD(3) * t64 + t130) * t115 - t67;
t107 = qJ(1) + qJ(2);
t101 = qJ(3) + t107;
t89 = pkin(9) + t101;
t77 = cos(t89);
t197 = g(2) * t77;
t108 = sin(pkin(9));
t109 = cos(pkin(9));
t163 = qJD(3) * t115;
t154 = t112 * t163;
t164 = qJD(3) * t111;
t43 = t115 * t49;
t22 = -t64 * t164 + t43 + (-t111 * t161 + (-t111 * t165 - t154) * qJD(1)) * pkin(1);
t96 = qJDD(3) + t103;
t15 = t96 * pkin(3) + t22;
t153 = t108 * t21 - t109 * t15;
t5 = -pkin(4) * t96 + t153;
t207 = t5 + t197;
t118 = qJD(5) ^ 2;
t174 = t109 * t111;
t180 = pkin(2) * qJD(3);
t171 = t112 * t115;
t138 = -t111 * t116 - t171;
t57 = t138 * t181;
t172 = t111 * t112;
t137 = t115 * t116 - t172;
t58 = t137 * t181;
t185 = -t108 * t58 + t109 * t57 + (t108 * t115 + t174) * t180;
t97 = qJD(3) + t104;
t155 = t185 * t97;
t175 = t108 * t111;
t92 = pkin(2) * t115 + pkin(3);
t55 = -pkin(2) * t175 + t109 * t92;
t50 = -pkin(4) - t55;
t56 = pkin(2) * t174 + t108 * t92;
t51 = pkin(8) + t56;
t205 = -t118 * t51 - t50 * t96 - t155;
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t41 = -t111 * t158 + t115 * t64;
t38 = pkin(3) * t97 + t41;
t42 = t111 * t64 + t115 * t158;
t39 = t109 * t42;
t24 = t108 * t38 + t39;
t20 = pkin(8) * t97 + t24;
t13 = qJD(4) * t114 - t110 * t20;
t178 = t114 * t20;
t14 = qJD(4) * t110 + t178;
t170 = t13 * qJD(5);
t8 = t108 * t15 + t109 * t21;
t6 = pkin(8) * t96 + t8;
t2 = t110 * qJDD(4) + t114 * t6 + t170;
t169 = t14 * qJD(5);
t98 = t114 * qJDD(4);
t3 = -t110 * t6 - t169 + t98;
t124 = -t3 * t110 + t2 * t114 + (-t110 * t14 - t114 * t13) * qJD(5);
t76 = sin(t89);
t183 = g(1) * t77 + g(2) * t76;
t90 = sin(t101);
t91 = cos(t101);
t204 = g(1) * t91 + g(2) * t90;
t184 = -t108 * t57 - t109 * t58 + (t109 * t115 - t175) * t180;
t203 = t184 * t97;
t202 = g(1) * t90 - g(2) * t91;
t100 = cos(t107);
t99 = sin(t107);
t201 = g(1) * t99 - g(2) * t100;
t105 = t110 ^ 2;
t106 = t114 ^ 2;
t166 = t105 + t106;
t200 = pkin(2) * t96;
t199 = pkin(2) * t99;
t198 = pkin(3) * t90;
t74 = g(1) * t76;
t93 = pkin(2) + t193;
t34 = t93 * t163 + (t137 * qJD(2) - t112 * t164) * pkin(1);
t35 = -t93 * t164 + (t138 * qJD(2) - t154) * pkin(1);
t9 = t108 * t34 - t109 * t35;
t195 = t9 * t97;
t113 = sin(qJ(1));
t194 = pkin(1) * t113;
t192 = pkin(3) * t108;
t191 = pkin(3) * t109;
t10 = t108 * t35 + t109 * t34;
t189 = t10 * t97;
t25 = t108 * t41 + t39;
t188 = t25 * t97;
t179 = t108 * t42;
t26 = t109 * t41 - t179;
t187 = t26 * t97;
t23 = t109 * t38 - t179;
t19 = -pkin(4) * t97 - t23;
t186 = t19 * qJD(5) * t110 + t114 * t74;
t59 = -pkin(1) * t172 + t115 * t93;
t54 = pkin(3) + t59;
t60 = pkin(1) * t171 + t111 * t93;
t32 = t108 * t54 + t109 * t60;
t182 = g(1) * t100 + g(2) * t99;
t117 = cos(qJ(1));
t102 = t117 * pkin(1);
t88 = pkin(2) * t100;
t177 = t88 + t102;
t173 = t110 * t114;
t168 = qJDD(4) - g(3);
t167 = t105 - t106;
t162 = qJD(5) * t114;
t160 = t207 * t110 + t19 * t162;
t81 = pkin(3) * t91;
t159 = t77 * pkin(4) + t76 * pkin(8) + t81;
t95 = t97 ^ 2;
t157 = t95 * t173;
t151 = t166 * t96;
t150 = t88 + t159;
t149 = t183 - t8;
t148 = qJD(1) * (-qJD(2) + t104);
t147 = qJD(2) * (-qJD(1) - t104);
t145 = t110 * t97 * t162;
t144 = t94 + t201;
t143 = -t194 - t199;
t142 = g(1) * t113 - g(2) * t117;
t31 = -t108 * t60 + t109 * t54;
t139 = t110 * t13 - t114 * t14;
t136 = -pkin(4) * t76 + t77 * pkin(8) - t198;
t135 = -t153 + t74 - t197;
t27 = -pkin(4) - t31;
t28 = pkin(8) + t32;
t133 = t118 * t28 + t27 * t96 + t195;
t83 = pkin(8) + t192;
t84 = -pkin(4) - t191;
t132 = t118 * t83 + t84 * t96 - t188;
t131 = t136 - t199;
t129 = -qJDD(5) * t28 + (t27 * t97 - t10) * qJD(5);
t128 = -qJDD(5) * t83 + (t84 * t97 + t26) * qJD(5);
t127 = -qJDD(5) * t51 + (t50 * t97 - t184) * qJD(5);
t125 = -qJD(4) * qJD(5) - t19 * t97 + t183 - t6;
t123 = -t183 + t124;
t122 = (-pkin(2) * t97 - t64) * qJD(3) - t130;
t121 = t204 - t21;
t120 = t22 + t202;
t66 = qJDD(5) * t114 - t110 * t118;
t65 = qJDD(5) * t110 + t114 * t118;
t45 = t106 * t96 - 0.2e1 * t145;
t44 = t105 * t96 + 0.2e1 * t145;
t36 = -0.2e1 * t167 * t97 * qJD(5) + 0.2e1 * t96 * t173;
t1 = [0, 0, 0, 0, 0, qJDD(1), t142, g(1) * t117 + g(2) * t113, 0, 0, 0, 0, 0, 0, 0, t103, (t103 * t116 + t112 * t147) * pkin(1) + t144, ((-qJDD(1) - t103) * t112 + t116 * t147) * pkin(1) + t182, 0, (t142 + (t112 ^ 2 + t116 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t96, t35 * t97 + t59 * t96 + t120, -t34 * t97 - t60 * t96 + t121, 0, -g(1) * t143 - g(2) * t177 + t21 * t60 + t22 * t59 + t42 * t34 + t41 * t35, 0, 0, 0, 0, 0, t96, t31 * t96 + t135 - t195, -t32 * t96 + t149 - t189, 0, t8 * t32 + t24 * t10 - t153 * t31 - t23 * t9 - g(1) * (t143 - t198) - g(2) * (t81 + t177), t44, t36, t65, t45, t66, 0, t129 * t110 + (-t133 - t207) * t114 + t186, t129 * t114 + (t133 - t74) * t110 + t160, t28 * t151 + t166 * t189 + t123, t5 * t27 + t19 * t9 - g(1) * (t131 - t194) - g(2) * (t102 + t150) - t139 * t10 + t124 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t112 * pkin(1) * t148 + t144, (t116 * t148 - t161) * pkin(1) + t182, 0, 0, 0, 0, 0, 0, 0, t96, -t57 * t97 + t43 + (-t146 + t200) * t115 + t122 * t111 + t202, t58 * t97 + t67 + (-t49 - t200) * t111 + t122 * t115 + t204, 0, -t41 * t57 - t42 * t58 + (t111 * t21 + t115 * t22 + (-t111 * t41 + t115 * t42) * qJD(3) + t201) * pkin(2), 0, 0, 0, 0, 0, t96, t55 * t96 + t135 - t155, -t56 * t96 + t149 - t203, 0, t8 * t56 - t153 * t55 - g(1) * (-t198 - t199) - g(2) * (t81 + t88) + t184 * t24 - t185 * t23, t44, t36, t65, t45, t66, 0, t127 * t110 + (-t207 + t205) * t114 + t186, t127 * t114 + (-t74 - t205) * t110 + t160, t51 * t151 + t166 * t203 + t123, t5 * t50 - g(1) * t131 - g(2) * t150 + t185 * t19 + ((t2 - t170) * t51 + t184 * t14) * t114 + ((-t3 - t169) * t51 - t184 * t13) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t42 * t97 + t120, t41 * t97 + t121, 0, 0, 0, 0, 0, 0, 0, t96, t96 * t191 + t135 + t188, -t96 * t192 + t149 + t187, 0, t23 * t25 - t24 * t26 + (t108 * t8 - t109 * t153 + t202) * pkin(3), t44, t36, t65, t45, t66, 0, t128 * t110 + (-t132 - t207) * t114 + t186, t128 * t114 + (t132 - t74) * t110 + t160, t123 + t166 * (t83 * t96 - t187), -g(1) * t136 - g(2) * t159 + t124 * t83 + t139 * t26 - t19 * t25 + t5 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, 0, 0, 0, 0, 0, 0, t66, -t65, 0, -qJD(5) * t139 + t2 * t110 + t3 * t114 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t167 * t95, t110 * t96, t157, t114 * t96, qJDD(5), -g(3) * t114 + t98 + (t14 - t178) * qJD(5) + t125 * t110, t170 + (qJD(5) * t20 - t168) * t110 + t125 * t114, 0, 0;];
tau_reg = t1;
