% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:05
% DurationCPUTime: 1.58s
% Computational Cost: add. (1918->281), mult. (3139->329), div. (0->0), fcn. (1791->12), ass. (0->174)
t111 = qJDD(1) + qJDD(2);
t116 = sin(pkin(8));
t117 = cos(pkin(8));
t123 = cos(qJ(2));
t197 = pkin(1) * qJD(2);
t163 = qJD(1) * t197;
t120 = sin(qJ(2));
t173 = qJDD(1) * t120;
t213 = pkin(1) * t173 + t123 * t163;
t202 = t123 * pkin(1);
t99 = qJDD(1) * t202;
t44 = t111 * pkin(2) - t120 * t163 + t99;
t21 = t116 * t44 + t213 * t117;
t15 = t111 * pkin(7) + t21;
t214 = qJD(3) * qJD(4) + t15;
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t104 = t122 * qJD(3);
t112 = qJD(1) + qJD(2);
t198 = pkin(1) * qJD(1);
t166 = t120 * t198;
t165 = t123 * t198;
t60 = t112 * pkin(2) + t165;
t35 = t116 * t60 + t117 * t166;
t30 = t112 * pkin(7) + t35;
t23 = -t119 * t30 + t104;
t175 = t119 * qJD(3);
t24 = t122 * t30 + t175;
t177 = qJD(4) * t119;
t5 = t119 * qJDD(3) + t214 * t122 - t30 * t177;
t101 = t122 * qJDD(3);
t6 = -t24 * qJD(4) - t119 * t15 + t101;
t128 = -t6 * t119 + t5 * t122 + (-t119 * t24 - t122 * t23) * qJD(4);
t115 = qJ(1) + qJ(2);
t102 = pkin(8) + t115;
t88 = sin(t102);
t81 = g(2) * t88;
t89 = cos(t102);
t83 = g(1) * t89;
t200 = -t81 - t83;
t105 = sin(t115);
t106 = cos(t115);
t212 = g(1) * t105 - g(2) * t106;
t157 = qJ(5) * t112 + t30;
t140 = t157 * t122;
t82 = g(1) * t88;
t211 = g(2) * t89;
t210 = pkin(2) * t105;
t113 = t119 ^ 2;
t209 = pkin(4) * t113;
t207 = g(3) * t122;
t206 = t116 * pkin(2);
t205 = t117 * pkin(2);
t121 = sin(qJ(1));
t204 = t121 * pkin(1);
t203 = t122 * pkin(4);
t16 = -t157 * t119 + t104;
t196 = qJD(4) * pkin(4);
t13 = t16 + t196;
t201 = t13 - t16;
t182 = t117 * t120;
t98 = pkin(2) + t202;
t54 = pkin(1) * t182 + t116 * t98;
t199 = g(1) * t106 + g(2) * t105;
t97 = pkin(3) + t203;
t61 = -t97 - t205;
t195 = t111 * t61;
t194 = t119 * t89;
t193 = t122 * t88;
t133 = pkin(1) * (t116 * t123 + t182);
t49 = qJD(1) * t133;
t192 = t49 * t112;
t50 = qJD(2) * t133;
t191 = t50 * t112;
t183 = t116 * t120;
t52 = (t117 * t123 - t183) * t197;
t190 = t52 * t112;
t48 = pkin(7) + t54;
t189 = -qJ(5) - t48;
t90 = pkin(7) + t206;
t188 = -qJ(5) - t90;
t187 = qJDD(4) * pkin(4);
t110 = t112 ^ 2;
t186 = t110 * t122;
t185 = t112 * t119;
t184 = t112 * t122;
t86 = t119 * t111;
t87 = t122 * t111;
t77 = t116 * t166;
t34 = t117 * t60 - t77;
t22 = -t97 * t112 + qJD(5) - t34;
t181 = -qJD(5) - t22;
t114 = t122 ^ 2;
t180 = t113 - t114;
t179 = t113 + t114;
t178 = qJD(4) * t112;
t176 = qJD(4) * t122;
t103 = t122 * qJD(5);
t160 = t112 * t177;
t3 = t112 * t103 + (-t160 + t87) * qJ(5) + t5;
t172 = t3 * t122 + t200;
t66 = g(2) * t194;
t168 = t213 * t116 - t117 * t44;
t74 = pkin(4) * t160;
t8 = -t97 * t111 + qJDD(5) + t168 + t74;
t171 = t8 * t119 + t22 * t176 + t66;
t14 = -t111 * pkin(3) + t168;
t29 = -t112 * pkin(3) - t34;
t170 = t14 * t119 + t29 * t176 + t66;
t51 = t117 * t165 - t77;
t67 = g(1) * t193;
t169 = t51 * t177 + t49 * t184 + t67;
t96 = pkin(2) * t106;
t167 = t89 * pkin(3) + t88 * pkin(7) + t96;
t164 = -t8 - t211;
t161 = -t14 - t211;
t159 = t112 * t176;
t53 = -pkin(1) * t183 + t117 * t98;
t158 = t179 * t51;
t156 = qJD(4) * t189;
t155 = qJD(4) * t188;
t154 = t111 * t179;
t153 = qJD(1) * (-qJD(2) + t112);
t152 = qJD(2) * (-qJD(1) - t112);
t47 = -pkin(3) - t53;
t150 = t99 + t212;
t149 = t119 * t159;
t118 = -qJ(5) - pkin(7);
t148 = -t88 * t118 + t89 * t97 + t96;
t147 = -t21 - t200;
t124 = cos(qJ(1));
t146 = g(1) * t121 - g(2) * t124;
t145 = -t192 - t82;
t36 = pkin(4) * t177 + t50;
t41 = t47 - t203;
t144 = t111 * t41 + t112 * t36;
t125 = qJD(4) ^ 2;
t91 = -pkin(3) - t205;
t143 = -t111 * t91 - t125 * t90;
t17 = t175 + t140;
t142 = t13 * t119 - t17 * t122;
t141 = t23 * t119 - t24 * t122;
t139 = g(1) * t194 + t119 * t81 + t101 - t207;
t138 = -t88 * pkin(3) + t89 * pkin(7) - t210;
t137 = t168 - t82 + t211;
t136 = g(2) * t193 + g(3) * t119 + t122 * t83 - t5;
t135 = -qJDD(4) * t90 + t91 * t178;
t134 = -qJ(5) * t111 - t214;
t132 = -t89 * t118 - t88 * t97 - t210;
t131 = t111 * t47 + t125 * t48 + t191;
t130 = -qJDD(4) * t48 + (t112 * t47 - t52) * qJD(4);
t127 = t200 + t128;
t108 = t124 * pkin(1);
t107 = t122 * qJ(5);
t72 = t119 * t186;
t63 = qJDD(4) * t122 - t125 * t119;
t62 = qJDD(4) * t119 + t125 * t122;
t58 = t180 * t110;
t57 = t122 * t90 + t107;
t56 = t188 * t119;
t46 = t114 * t111 - 0.2e1 * t149;
t45 = t113 * t111 + 0.2e1 * t149;
t43 = -t119 * qJD(5) + t122 * t155;
t42 = t119 * t155 + t103;
t38 = t51 * t176;
t32 = t122 * t48 + t107;
t31 = t189 * t119;
t28 = 0.2e1 * t119 * t87 - 0.2e1 * t180 * t178;
t25 = t29 * t177;
t18 = t22 * t177;
t10 = (-qJD(5) - t52) * t119 + t122 * t156;
t9 = t119 * t156 + t122 * t52 + t103;
t2 = t187 + t101 - qJD(4) * t140 + (-qJD(5) * t112 + t134) * t119;
t1 = [0, 0, 0, 0, 0, qJDD(1), t146, g(1) * t124 + g(2) * t121, 0, 0, 0, 0, 0, 0, 0, t111, (t111 * t123 + t120 * t152) * pkin(1) + t150, ((-qJDD(1) - t111) * t120 + t123 * t152) * pkin(1) + t199, 0, (t146 + (t120 ^ 2 + t123 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t111, t53 * t111 - t137 - t191, -t54 * t111 + t147 - t190, 0, t21 * t54 + t35 * t52 - t168 * t53 - t34 * t50 - g(1) * (-t204 - t210) - g(2) * (t96 + t108), t45, t28, t62, t46, t63, 0, t25 + t67 + t130 * t119 + (-t131 + t161) * t122, t130 * t122 + (t131 - t82) * t119 + t170, t48 * t154 + t179 * t190 + t127, t14 * t47 + t29 * t50 - g(1) * (t138 - t204) - g(2) * (t108 + t167) - t141 * t52 + t128 * t48, t45, t28, t62, t46, t63, 0, t31 * qJDD(4) + t18 + t67 + (t41 * t185 + t10) * qJD(4) + (-t144 + t164) * t122, -t32 * qJDD(4) + (t41 * t184 - t9) * qJD(4) + (t144 - t82) * t119 + t171, (t111 * t32 + t112 * t9 + (-t112 * t31 - t13) * qJD(4)) * t122 + (-t10 * t112 - t111 * t31 - t2 + (-t112 * t32 - t17) * qJD(4)) * t119 + t172, t3 * t32 + t17 * t9 + t2 * t31 + t13 * t10 + t8 * t41 + t22 * t36 - g(1) * (t132 - t204) - g(2) * (t108 + t148); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t120 * pkin(1) * t153 + t150, (t123 * t153 - t173) * pkin(1) + t199, 0, 0, 0, 0, 0, 0, 0, t111, t111 * t205 - t137 + t192, -t111 * t206 + t51 * t112 + t147, 0, t34 * t49 - t35 * t51 + (t116 * t21 - t117 * t168 + t212) * pkin(2), t45, t28, t62, t46, t63, 0, t25 + t135 * t119 + (t143 + t161) * t122 + t169, t38 + t135 * t122 + (-t143 + t145) * t119 + t170, -t112 * t158 + t90 * t154 + t127, -g(1) * t138 - g(2) * t167 + t128 * t90 + t14 * t91 + t141 * t51 - t29 * t49, t45, t28, t62, t46, t63, 0, t56 * qJDD(4) + t18 + (t61 * t185 + t43) * qJD(4) + (t164 - t74 - t195) * t122 + t169, -t57 * qJDD(4) + t38 + (t145 + t195) * t119 + (-t42 + (t122 * t61 + t209) * t112) * qJD(4) + t171, (-qJD(4) * t13 + t111 * t57) * t122 + (-t17 * qJD(4) - t111 * t56 - t2) * t119 + (-t119 * t43 + t122 * t42 - t158 + (-t119 * t57 - t122 * t56) * qJD(4)) * t112 + t172, t3 * t57 + t2 * t56 + t13 * t43 + t8 * t61 - t22 * t49 - g(1) * t132 - g(2) * t148 + (-t122 * t51 + t42) * t17 + (t13 * t51 + t22 * t196) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0, 0, 0, 0, t63, -t62, 0, -t141 * qJD(4) + t5 * t119 + t6 * t122 - g(3), 0, 0, 0, 0, 0, 0, t63, -t62, 0, -qJD(4) * t142 + t3 * t119 + t2 * t122 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t58, t86, t72, t87, qJDD(4), (-t112 * t29 - t15) * t119 + t139, t23 * qJD(4) - t29 * t184 + t136, 0, 0, -t72, t58, t86, t72, t87, qJDD(4), 0.2e1 * t187 + (t17 - t140) * qJD(4) + (pkin(4) * t186 + t181 * t112 + t134) * t119 + t139, -t110 * t209 - qJ(5) * t87 + t16 * qJD(4) + (qJ(5) * t177 + t181 * t122) * t112 + t136, -pkin(4) * t86 + (-t196 + t201) * t184, t201 * t17 + (-t207 + t2 + (-t112 * t22 - t200) * t119) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 + 0.2e1 * t160, t86 + 0.2e1 * t159, -t179 * t110, t112 * t142 - t164 - t82;];
tau_reg = t1;
