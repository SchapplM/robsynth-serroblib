% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:05
% EndTime: 2019-12-05 18:00:08
% DurationCPUTime: 1.81s
% Computational Cost: add. (2496->286), mult. (4790->320), div. (0->0), fcn. (3015->8), ass. (0->160)
t128 = qJ(3) + qJ(4);
t111 = sin(t128);
t112 = cos(t128);
t131 = sin(qJ(1));
t134 = cos(qJ(1));
t186 = g(1) * t131 - g(2) * t134;
t208 = g(3) * t111 - t112 * t186;
t138 = qJD(1) ^ 2;
t213 = -t138 * qJ(2) - t186;
t133 = cos(qJ(3));
t176 = qJD(1) * qJD(3);
t164 = t133 * t176;
t130 = sin(qJ(3));
t174 = t130 * qJDD(1);
t212 = t164 + t174;
t173 = t133 * qJDD(1);
t211 = t130 * t176 - t173;
t123 = qJDD(1) * qJ(2);
t132 = cos(qJ(4));
t129 = sin(qJ(4));
t182 = qJD(1) * t130;
t136 = -pkin(1) - pkin(6);
t91 = t136 * qJD(1) + qJD(2);
t58 = -pkin(7) * t182 + t130 * t91;
t52 = t129 * t58;
t181 = qJD(1) * t133;
t59 = -pkin(7) * t181 + t133 * t91;
t55 = qJD(3) * pkin(3) + t59;
t29 = t132 * t55 - t52;
t166 = t129 * t182;
t68 = t132 * t181 - t166;
t60 = t68 * qJ(5);
t19 = t29 - t60;
t119 = qJDD(3) + qJDD(4);
t121 = qJD(3) + qJD(4);
t177 = qJD(4) * t132;
t169 = pkin(3) * t177;
t204 = pkin(3) * t129;
t210 = -t119 * t204 - t121 * t169;
t76 = t129 * t133 + t130 * t132;
t66 = t76 * qJD(1);
t194 = t121 * t66;
t126 = t130 ^ 2;
t127 = t133 ^ 2;
t184 = t126 + t127;
t90 = t136 * qJDD(1) + qJDD(2);
t161 = t184 * t90;
t109 = t119 * pkin(4);
t153 = t129 * t174 - t132 * t173;
t42 = t121 * t76;
t26 = t42 * qJD(1) + t153;
t197 = qJ(5) * t26;
t209 = t109 + t197;
t156 = g(1) * t134 + g(2) * t131;
t124 = qJD(1) * qJD(2);
t168 = 0.2e1 * t124;
t207 = 0.2e1 * t123 + t168 - t156;
t206 = t68 ^ 2;
t135 = -pkin(7) - pkin(6);
t203 = pkin(3) * t132;
t202 = g(3) * t130;
t115 = t130 * pkin(3);
t201 = t68 * t66;
t200 = pkin(7) - t136;
t16 = pkin(4) * t121 + t19;
t199 = t16 - t19;
t152 = -qJD(4) * t166 - t211 * t129;
t160 = t121 * t133;
t27 = (qJD(1) * t160 + t174) * t132 + t152;
t198 = -t66 * t169 - t27 * t204;
t77 = -t129 * t130 + t132 * t133;
t24 = t77 * t119 - t42 * t121;
t36 = t132 * t59 - t52;
t82 = t200 * t130;
t83 = t200 * t133;
t45 = -t129 * t83 - t132 * t82;
t196 = qJ(5) * t27;
t195 = qJ(5) * t66;
t193 = t121 * t68;
t53 = t132 * t58;
t192 = pkin(1) * qJDD(1);
t84 = pkin(3) * t182 + qJD(1) * qJ(2);
t191 = qJD(1) * t84;
t99 = qJ(2) + t115;
t163 = -pkin(4) * t66 - qJD(5);
t46 = -t163 + t84;
t189 = qJD(5) + t46;
t188 = (t168 + t123) * qJ(2);
t187 = t134 * pkin(1) + t131 * qJ(2);
t185 = t126 - t127;
t137 = qJD(3) ^ 2;
t183 = -t137 - t138;
t180 = qJD(3) * t130;
t179 = qJD(3) * t133;
t178 = qJD(4) * t129;
t93 = pkin(3) * t179 + qJD(2);
t175 = qJDD(3) * t130;
t171 = pkin(3) * t181;
t170 = pkin(3) * t178;
t167 = t133 * t138 * t130;
t78 = t133 * t90;
t38 = qJDD(3) * pkin(3) + t211 * pkin(7) - t91 * t180 + t78;
t40 = -t212 * pkin(7) + t130 * t90 + t91 * t179;
t162 = -t129 * t40 + t132 * t38;
t35 = -t129 * t59 - t53;
t44 = t129 * t82 - t132 * t83;
t5 = t129 * t38 + t132 * t40 + t55 * t177 - t58 * t178;
t159 = t184 * qJDD(1);
t56 = t212 * pkin(3) + t123 + t124;
t158 = qJDD(2) - t192;
t157 = t130 * t164;
t9 = -t26 * t77 - t42 * t68;
t43 = -t129 * t180 - t130 * t178 + t132 * t160;
t10 = t27 * t76 + t43 * t66;
t25 = -t119 * t76 - t121 * t43;
t30 = t129 * t55 + t53;
t151 = t156 * t111;
t150 = t156 * t112;
t74 = t200 * t180;
t75 = qJD(3) * t83;
t17 = t129 * t74 - t132 * t75 - t83 * t177 + t82 * t178;
t149 = 0.2e1 * qJ(2) * t176 + qJDD(3) * t136;
t11 = pkin(4) * t27 + qJDD(5) + t56;
t147 = g(3) * t112 + t186 * t111 - t5;
t6 = -qJD(4) * t30 + t162;
t18 = -t45 * qJD(4) + t129 * t75 + t132 * t74;
t20 = t30 - t195;
t3 = -qJD(5) * t68 + t209 + t6;
t4 = -qJD(5) * t66 - t196 + t5;
t146 = -t16 * t42 + t20 * t43 + t3 * t77 + t4 * t76 - t186;
t145 = -t29 * t42 + t30 * t43 + t5 * t76 + t6 * t77 - t186;
t144 = t66 * t84 + t147;
t143 = -t136 * t137 + t207;
t142 = t189 * t66 + t147 + t196;
t141 = t6 + t208;
t140 = -t68 * t84 + t141;
t139 = -t153 - t194;
t120 = -qJ(5) + t135;
t114 = t134 * qJ(2);
t110 = qJDD(3) * t133;
t103 = pkin(4) + t203;
t80 = pkin(4) * t111 + t115;
t65 = t66 ^ 2;
t57 = pkin(4) * t76 + t99;
t49 = pkin(4) * t68 + t171;
t39 = pkin(4) * t43 + t93;
t34 = -qJ(5) * t76 + t45;
t33 = -qJ(5) * t77 + t44;
t28 = -t65 + t206;
t22 = -t60 + t36;
t21 = t35 + t195;
t15 = -qJD(1) * t68 + t25;
t14 = -qJD(1) * t66 + t24;
t13 = t193 + (-t121 * t181 - t174) * t132 - t152;
t12 = t139 + t194;
t8 = qJ(5) * t42 - qJD(5) * t77 + t18;
t7 = -qJ(5) * t43 - qJD(5) * t76 + t17;
t2 = -t10 - t9;
t1 = t26 * t76 - t27 * t77 + t42 * t66 - t43 * t68;
t23 = [0, 0, 0, 0, 0, qJDD(1), t186, t156, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t186 - 0.2e1 * t192, t207, -t158 * pkin(1) - g(1) * (-pkin(1) * t131 + t114) - g(2) * t187 + t188, qJDD(1) * t127 - 0.2e1 * t157, -0.2e1 * t130 * t173 + 0.2e1 * t185 * t176, -t130 * t137 + t110, qJDD(1) * t126 + 0.2e1 * t157, -t133 * t137 - t175, 0, t130 * t143 + t133 * t149, -t130 * t149 + t133 * t143, -t136 * t159 - t161 + t186, -g(1) * (t136 * t131 + t114) - g(2) * (pkin(6) * t134 + t187) + t136 * t161 + t188, t9, t1, t24, t10, t25, 0, t119 * t44 + t121 * t18 + t27 * t99 + t43 * t84 + t56 * t76 + t66 * t93 - t151, -t119 * t45 - t121 * t17 - t26 * t99 - t42 * t84 + t56 * t77 + t68 * t93 - t150, -t17 * t66 - t18 * t68 + t26 * t44 - t27 * t45 - t145, t5 * t45 + t30 * t17 + t6 * t44 + t29 * t18 + t56 * t99 + t84 * t93 - g(1) * (t134 * t115 + t114 + (-pkin(1) + t135) * t131) - g(2) * (t115 * t131 - t134 * t135 + t187), t9, t1, t24, t10, t25, 0, t11 * t76 + t119 * t33 + t121 * t8 + t27 * t57 + t39 * t66 + t43 * t46 - t151, t11 * t77 - t119 * t34 - t121 * t7 - t26 * t57 + t39 * t68 - t42 * t46 - t150, t26 * t33 - t27 * t34 - t66 * t7 - t68 * t8 - t146, t4 * t34 + t20 * t7 + t3 * t33 + t16 * t8 + t11 * t57 + t46 * t39 - g(1) * (t134 * t80 + t114 + (-pkin(1) + t120) * t131) - g(2) * (-t120 * t134 + t131 * t80 + t187); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t138, t213 + t158, 0, 0, 0, 0, 0, 0, t183 * t130 + t110, t183 * t133 - t175, -t159, t161 + t213, 0, 0, 0, 0, 0, 0, t14, t15, t2, t145 - t191, 0, 0, 0, 0, 0, 0, t14, t15, t2, -qJD(1) * t46 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t185 * t138, t173, -t167, -t174, qJDD(3), t133 * t213 + t202 + t78, g(3) * t133 + (-t213 - t90) * t130, 0, 0, t201, t28, t12, -t201, t13, t119, -t121 * t35 + (t119 * t132 - t121 * t178 - t66 * t181) * pkin(3) + t140, t121 * t36 - t68 * t171 + t144 + t210, t26 * t203 + (-t29 + t36) * t66 + (t30 + t35 + t170) * t68 + t198, -t29 * t35 - t30 * t36 + (t202 + t129 * t5 + t132 * t6 + (-t129 * t29 + t132 * t30) * qJD(4) + (-t186 - t191) * t133) * pkin(3), t201, t28, t12, -t201, t13, t119, t103 * t119 - t121 * t21 - t49 * t66 - t189 * t68 + (-t53 + (-pkin(3) * t121 - t55) * t129) * qJD(4) + t162 + t208 + t209, t121 * t22 - t49 * t68 + t142 + t210, t103 * t26 + (-t16 + t22) * t66 + (t20 + t21 + t170) * t68 + t198, g(3) * t80 + t103 * t3 - t16 * t21 - t20 * t22 - t46 * t49 - t186 * (pkin(3) * t133 + pkin(4) * t112) + (t129 * t4 + (-t129 * t16 + t132 * t20) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t28, t12, -t201, t13, t119, t121 * t30 + t140, t121 * t29 + t144, 0, 0, t201, t28, t12, -t201, t13, t119, t197 + t121 * t20 + 0.2e1 * t109 + (t163 - t46) * t68 + t141, -pkin(4) * t206 + t121 * t19 + t142, t26 * pkin(4) - t199 * t66, t199 * t20 + (-t46 * t68 + t208 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 + t193, t139 - t194, -t65 - t206, t16 * t68 + t20 * t66 + t11 - t156;];
tau_reg = t23;
