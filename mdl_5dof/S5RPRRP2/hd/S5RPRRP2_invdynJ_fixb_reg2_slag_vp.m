% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:13
% EndTime: 2022-01-23 09:28:16
% DurationCPUTime: 1.29s
% Computational Cost: add. (1782->251), mult. (3137->297), div. (0->0), fcn. (1761->12), ass. (0->153)
t101 = qJDD(1) + qJDD(3);
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t107 = sin(pkin(8));
t172 = pkin(1) * qJDD(1);
t152 = t107 * t172;
t108 = cos(pkin(8));
t87 = t108 * pkin(1) + pkin(2);
t65 = t87 * qJD(1);
t193 = qJD(3) * t65 + t152;
t190 = pkin(1) * t107;
t154 = qJD(1) * t190;
t194 = qJD(3) * t154 - t87 * qJDD(1);
t144 = -t194 * t111 + t193 * t114;
t13 = t101 * pkin(7) + t144;
t195 = qJD(2) * qJD(4) + t13;
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t102 = qJD(1) + qJD(3);
t37 = t111 * t65 + t114 * t154;
t29 = t102 * pkin(7) + t37;
t96 = t113 * qJD(2);
t21 = -t110 * t29 + t96;
t160 = t110 * qJD(2);
t22 = t113 * t29 + t160;
t162 = qJD(4) * t110;
t5 = t110 * qJDD(2) + t195 * t113 - t29 * t162;
t93 = t113 * qJDD(2);
t6 = -t22 * qJD(4) - t110 * t13 + t93;
t119 = -t6 * t110 + t5 * t113 + (-t110 * t22 - t113 * t21) * qJD(4);
t103 = qJ(1) + pkin(8);
t94 = qJ(3) + t103;
t85 = sin(t94);
t77 = g(2) * t85;
t86 = cos(t94);
t79 = g(1) * t86;
t183 = -t77 - t79;
t149 = t102 * t162;
t186 = t113 * pkin(4);
t89 = pkin(3) + t186;
t180 = t101 * t89;
t192 = -pkin(4) * t149 + t180;
t44 = -t111 * t190 + t114 * t87;
t142 = qJ(5) * t102 + t29;
t128 = t142 * t113;
t78 = g(1) * t85;
t191 = g(2) * t86;
t104 = t110 ^ 2;
t189 = pkin(4) * t104;
t188 = g(3) * t113;
t187 = t101 * pkin(3);
t109 = -qJ(5) - pkin(7);
t16 = -t142 * t110 + t96;
t181 = qJD(4) * pkin(4);
t15 = t16 + t181;
t185 = t15 - t16;
t184 = t86 * pkin(3) + t85 * pkin(7);
t45 = t111 * t87 + t114 * t190;
t115 = cos(qJ(1));
t91 = cos(t103);
t182 = t115 * pkin(1) + pkin(2) * t91;
t179 = t110 * t86;
t178 = t113 * t85;
t176 = t37 * t102;
t38 = t44 * qJD(3);
t175 = t38 * t102;
t39 = t45 * qJD(3);
t174 = t39 * t102;
t43 = pkin(7) + t45;
t173 = -qJ(5) - t43;
t170 = qJDD(4) * pkin(4);
t100 = t102 ^ 2;
t169 = t100 * t113;
t168 = t102 * t110;
t167 = t102 * t113;
t82 = t110 * t101;
t83 = t113 * t101;
t36 = -t111 * t154 + t114 * t65;
t20 = -t89 * t102 + qJD(5) - t36;
t166 = -qJD(5) - t20;
t105 = t113 ^ 2;
t165 = t104 - t105;
t164 = t104 + t105;
t163 = qJD(4) * t102;
t161 = qJD(4) * t113;
t95 = t113 * qJD(5);
t3 = t102 * t95 + (-t149 + t83) * qJ(5) + t5;
t158 = t3 * t113 + t183;
t58 = g(2) * t179;
t143 = t193 * t111 + t194 * t114;
t8 = qJDD(5) + t143 - t192;
t157 = t8 * t110 + t20 * t161 + t58;
t14 = t143 - t187;
t28 = -t102 * pkin(3) - t36;
t156 = t14 * t110 + t28 * t161 + t58;
t59 = g(1) * t178;
t155 = t36 * t162 + t37 * t167 + t59;
t153 = -t8 - t191;
t151 = -t85 * pkin(3) + t86 * pkin(7);
t150 = -t14 - t191;
t148 = t102 * t161;
t147 = -t85 * t109 + t86 * t89;
t146 = qJD(4) * t109;
t145 = t164 * t36;
t141 = qJD(4) * t173;
t140 = t164 * t101;
t42 = -pkin(3) - t44;
t138 = t110 * t148;
t112 = sin(qJ(1));
t90 = sin(t103);
t137 = -t112 * pkin(1) - pkin(2) * t90;
t116 = qJD(4) ^ 2;
t136 = -pkin(7) * t116 + t187;
t135 = g(1) * t112 - g(2) * t115;
t134 = -t176 - t78;
t133 = -t86 * t109 - t85 * t89;
t34 = pkin(4) * t162 + t39;
t35 = t42 - t186;
t132 = t101 * t35 + t102 * t34;
t17 = t160 + t128;
t131 = t15 * t110 - t17 * t113;
t130 = t21 * t110 - t22 * t113;
t129 = g(1) * t179 + t110 * t77 - t188 + t93;
t127 = -t144 - t183;
t126 = g(2) * t178 + g(3) * t110 + t113 * t79 - t5;
t125 = -pkin(3) * t163 - pkin(7) * qJDD(4);
t124 = -qJ(5) * t101 - t195;
t123 = t143 - t78 + t191;
t122 = t101 * t42 + t116 * t43 + t174;
t121 = -qJDD(4) * t43 + (t102 * t42 - t38) * qJD(4);
t118 = t183 + t119;
t106 = qJDD(2) - g(3);
t97 = t113 * qJ(5);
t68 = t110 * t169;
t67 = t113 * pkin(7) + t97;
t66 = t109 * t110;
t55 = qJDD(4) * t113 - t116 * t110;
t54 = qJDD(4) * t110 + t116 * t113;
t48 = t165 * t100;
t47 = -t110 * qJD(5) + t113 * t146;
t46 = t110 * t146 + t95;
t41 = t105 * t101 - 0.2e1 * t138;
t40 = t104 * t101 + 0.2e1 * t138;
t33 = t113 * t43 + t97;
t32 = t173 * t110;
t31 = t36 * t161;
t26 = 0.2e1 * t110 * t83 - 0.2e1 * t165 * t163;
t23 = t28 * t162;
t18 = t20 * t162;
t10 = (-qJD(5) - t38) * t110 + t113 * t141;
t9 = t110 * t141 + t113 * t38 + t95;
t2 = t170 + t93 - qJD(4) * t128 + (-qJD(5) * t102 + t124) * t110;
t1 = [0, 0, 0, 0, 0, qJDD(1), t135, g(1) * t115 + g(2) * t112, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t90 - g(2) * t91 + 0.2e1 * t108 * t172, g(1) * t91 + g(2) * t90 - 0.2e1 * t152, 0, (t135 + (t107 ^ 2 + t108 ^ 2) * t172) * pkin(1), 0, 0, 0, 0, 0, t101, t44 * t101 - t123 - t174, -t45 * t101 + t127 - t175, 0, -g(1) * t137 - g(2) * t182 - t143 * t44 + t144 * t45 - t36 * t39 + t37 * t38, t40, t26, t54, t41, t55, 0, t23 + t59 + t121 * t110 + (-t122 + t150) * t113, t121 * t113 + (t122 - t78) * t110 + t156, t43 * t140 + t164 * t175 + t118, t14 * t42 + t28 * t39 - g(1) * (t137 + t151) - g(2) * (t182 + t184) - t130 * t38 + t119 * t43, t40, t26, t54, t41, t55, 0, t32 * qJDD(4) + t18 + t59 + (t35 * t168 + t10) * qJD(4) + (-t132 + t153) * t113, -t33 * qJDD(4) + (t35 * t167 - t9) * qJD(4) + (t132 - t78) * t110 + t157, (t101 * t33 + t102 * t9 + (-t102 * t32 - t15) * qJD(4)) * t113 + (-t10 * t102 - t101 * t32 - t2 + (-t102 * t33 - t17) * qJD(4)) * t110 + t158, t3 * t33 + t17 * t9 + t2 * t32 + t15 * t10 + t8 * t35 + t20 * t34 - g(1) * (t133 + t137) - g(2) * (t147 + t182); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, 0, t55, -t54, 0, -t130 * qJD(4) + t5 * t110 + t6 * t113 - g(3), 0, 0, 0, 0, 0, 0, t55, -t54, 0, -t131 * qJD(4) + t3 * t110 + t2 * t113 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t123 + t176, t36 * t102 + t127, 0, 0, t40, t26, t54, t41, t55, 0, t23 + t125 * t110 + (t136 + t150) * t113 + t155, t31 + t125 * t113 + (t134 - t136) * t110 + t156, pkin(7) * t140 - t102 * t145 + t118, -t14 * pkin(3) + t119 * pkin(7) - g(1) * t151 - g(2) * t184 + t130 * t36 - t28 * t37, t40, t26, t54, t41, t55, 0, t66 * qJDD(4) + t18 + (-t89 * t168 + t47) * qJD(4) + (t153 + t192) * t113 + t155, -t67 * qJDD(4) + t31 + (t134 - t180) * t110 + (-t46 + (-t113 * t89 + t189) * t102) * qJD(4) + t157, (-qJD(4) * t15 + t101 * t67) * t113 + (-t17 * qJD(4) - t101 * t66 - t2) * t110 + (-t110 * t47 + t113 * t46 - t145 + (-t110 * t67 - t113 * t66) * qJD(4)) * t102 + t158, t3 * t67 + t2 * t66 + t15 * t47 - t8 * t89 - t20 * t37 - g(1) * t133 - g(2) * t147 + (-t113 * t36 + t46) * t17 + (t15 * t36 + t20 * t181) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t48, t82, t68, t83, qJDD(4), (-t102 * t28 - t13) * t110 + t129, t21 * qJD(4) - t28 * t167 + t126, 0, 0, -t68, t48, t82, t68, t83, qJDD(4), 0.2e1 * t170 + (t17 - t128) * qJD(4) + (pkin(4) * t169 + t166 * t102 + t124) * t110 + t129, -t100 * t189 - qJ(5) * t83 + t16 * qJD(4) + (qJ(5) * t162 + t166 * t113) * t102 + t126, -pkin(4) * t82 + (-t181 + t185) * t167, t185 * t17 + (-t188 + t2 + (-t102 * t20 - t183) * t110) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 + 0.2e1 * t149, t82 + 0.2e1 * t148, -t164 * t100, t131 * t102 - t153 - t78;];
tau_reg = t1;
