% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:14
% DurationCPUTime: 1.49s
% Computational Cost: add. (2611->213), mult. (6331->296), div. (0->0), fcn. (4442->8), ass. (0->140)
t128 = sin(qJ(5));
t124 = sin(pkin(9));
t126 = cos(pkin(9));
t130 = cos(qJ(3));
t160 = qJD(1) * t130;
t150 = t126 * t160;
t129 = sin(qJ(3));
t161 = qJD(1) * t129;
t100 = t124 * t161 - t150;
t108 = t124 * t130 + t126 * t129;
t137 = qJD(1) * t108;
t181 = cos(qJ(5));
t140 = t128 * t100 - t137 * t181;
t86 = qJD(3) * t137;
t153 = qJD(1) * qJD(3);
t148 = t129 * t153;
t112 = t124 * t148;
t147 = t130 * t153;
t87 = t126 * t147 - t112;
t134 = qJD(5) * t140 - t128 * t87 - t181 * t86;
t121 = qJD(3) + qJD(5);
t169 = t140 * t121;
t191 = t134 - t169;
t149 = qJD(5) * t181;
t158 = qJD(5) * t128;
t139 = -t100 * t149 - t128 * t86 - t137 * t158 + t181 * t87;
t41 = -t181 * t100 - t128 * t137;
t172 = t121 * t41;
t190 = t139 - t172;
t175 = t41 ^ 2;
t176 = t140 ^ 2;
t189 = -t175 + t176;
t174 = t41 * t140;
t154 = t130 * qJD(4);
t159 = qJD(3) * t129;
t116 = sin(pkin(8)) * pkin(1) + pkin(6);
t110 = t116 * qJD(1);
t120 = t130 * qJD(2);
t119 = qJD(3) * t120;
t76 = -t110 * t159 + t119;
t48 = (-qJ(4) * t159 + t154) * qJD(1) + t76;
t155 = t129 * qJD(4);
t144 = qJ(4) * qJD(1) + t110;
t156 = t129 * qJD(2);
t68 = t130 * t144 + t156;
t49 = -qJD(1) * t155 - qJD(3) * t68;
t22 = -t124 * t48 + t126 * t49;
t14 = -pkin(7) * t87 + t22;
t23 = t124 * t49 + t126 * t48;
t15 = -pkin(7) * t86 + t23;
t178 = pkin(7) * t137;
t58 = t124 * t68;
t67 = -t129 * t144 + t120;
t62 = qJD(3) * pkin(3) + t67;
t30 = t126 * t62 - t58;
t19 = qJD(3) * pkin(4) - t178 + t30;
t179 = pkin(7) * t100;
t171 = t126 * t68;
t31 = t124 * t62 + t171;
t21 = t31 - t179;
t6 = t128 * t19 + t181 * t21;
t2 = -qJD(5) * t6 - t128 * t15 + t181 * t14;
t118 = -cos(pkin(8)) * pkin(1) - pkin(2);
t109 = -pkin(3) * t130 + t118;
t162 = qJD(1) * t109;
t98 = qJD(4) + t162;
t50 = t100 * pkin(4) + t98;
t188 = t140 * t50 + t2;
t135 = -t128 * t14 - t19 * t149 - t181 * t15 + t21 * t158;
t187 = -t50 * t41 + t135;
t186 = 0.2e1 * t137;
t185 = t137 ^ 2;
t131 = qJD(3) ^ 2;
t107 = t124 * t129 - t126 * t130;
t104 = t107 * qJD(3);
t138 = t108 * qJD(3);
t24 = t181 * t104 + t107 * t149 + t108 * t158 + t128 * t138;
t52 = -t128 * t107 + t181 * t108;
t184 = t134 * t52 - t24 * t41;
t32 = -t124 * t67 - t171;
t26 = t32 + t179;
t33 = t126 * t67 - t58;
t27 = t33 - t178;
t117 = pkin(3) * t126 + pkin(4);
t180 = pkin(3) * t124;
t96 = t128 * t117 + t181 * t180;
t183 = t96 * qJD(5) - t128 * t27 + t181 * t26;
t95 = t181 * t117 - t128 * t180;
t182 = t95 * qJD(5) - t128 * t26 - t181 * t27;
t177 = t129 * pkin(3);
t173 = t104 * t100 - t108 * t86;
t164 = qJ(4) + t116;
t145 = qJD(3) * t164;
t71 = -t129 * t145 + t154;
t72 = -t130 * t145 - t155;
t35 = t124 * t72 + t126 * t71;
t105 = t164 * t129;
t106 = t164 * t130;
t47 = -t124 * t105 + t126 * t106;
t170 = t24 * t121;
t168 = t137 * t100;
t167 = t110 * t129;
t166 = t131 * t129;
t165 = t131 * t130;
t163 = t129 ^ 2 - t130 ^ 2;
t111 = qJD(1) * t118;
t157 = t104 * qJD(3);
t152 = pkin(3) * t161;
t132 = qJD(1) ^ 2;
t151 = t129 * t132 * t130;
t114 = pkin(3) * t148;
t57 = pkin(4) * t86 + t114;
t34 = -t124 * t71 + t126 * t72;
t46 = -t126 * t105 - t106 * t124;
t143 = t129 * t147;
t25 = t52 * qJD(5) - t128 * t104 + t181 * t138;
t51 = t181 * t107 + t108 * t128;
t142 = t139 * t51 - t140 * t25;
t89 = t110 * t130 + t156;
t141 = 0.2e1 * qJD(3) * t111;
t36 = -pkin(7) * t108 + t46;
t37 = -pkin(7) * t107 + t47;
t11 = -t128 * t37 + t181 * t36;
t12 = t128 * t36 + t181 * t37;
t136 = -t87 * t107 - t137 * t138;
t77 = t89 * qJD(3);
t88 = t120 - t167;
t133 = t77 * t129 + t76 * t130 + (-t129 * t89 - t130 * t88) * qJD(3);
t97 = t100 ^ 2;
t93 = t108 * t131;
t70 = (pkin(4) * t108 + t177) * qJD(3);
t69 = pkin(4) * t137 + t152;
t66 = pkin(4) * t107 + t109;
t29 = -pkin(7) * t138 + t35;
t28 = pkin(7) * t104 + t34;
t20 = t25 * t121;
t5 = -t128 * t21 + t181 * t19;
t4 = -t12 * qJD(5) - t128 * t29 + t181 * t28;
t3 = t11 * qJD(5) + t128 * t28 + t181 * t29;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t143, -0.2e1 * t163 * t153, t165, -0.2e1 * t143, -t166, 0, -t116 * t165 + t129 * t141, t116 * t166 + t130 * t141, t133, t133 * t116, -t104 * t137 + t108 * t87, t136 + t173, -t157, t100 * t138 + t86 * t107, -t93, 0, t109 * t86 + (t98 * t108 + t34 + (qJD(1) * t107 + t100) * t177) * qJD(3), -t98 * t104 + t109 * t87 + (t177 * t186 - t35) * qJD(3), -t35 * t100 + t30 * t104 - t23 * t107 - t22 * t108 - t137 * t34 - t138 * t31 - t46 * t87 - t47 * t86, t22 * t46 + t23 * t47 + t30 * t34 + t31 * t35 + (t98 + t162) * pkin(3) * t159, t139 * t52 + t140 * t24, -t142 + t184, -t170, -t134 * t51 - t25 * t41, -t20, 0, t121 * t4 - t134 * t66 + t25 * t50 - t41 * t70 + t51 * t57, -t121 * t3 + t139 * t66 - t140 * t70 - t24 * t50 + t52 * t57, -t11 * t139 + t12 * t134 + t135 * t51 + t140 * t4 - t2 * t52 + t24 * t5 - t25 * t6 + t3 * t41, t11 * t2 - t12 * t135 + t3 * t6 + t4 * t5 + t50 * t70 + t57 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t165, 0, t76 * t129 - t77 * t130 + (-t129 * t88 + t130 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, -t93, t157, -t136 + t173, -t31 * t104 - t22 * t107 + t23 * t108 - t138 * t30, 0, 0, 0, 0, 0, 0, -t20, t170, t142 + t184, -t135 * t52 - t2 * t51 - t24 * t6 - t25 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t163 * t132, 0, t151, 0, 0, -t111 * t161, -t111 * t160 - t119 + (t88 + t167) * qJD(3), 0, 0, t168, -t97 + t185, -t112 + (t100 + t150) * qJD(3), -t168, 0, 0, -qJD(3) * t32 - t100 * t152 - t137 * t98 + t22, qJD(3) * t33 + t100 * t98 - t137 * t152 - t23, (t31 + t32) * t137 + (-t30 + t33) * t100 + (-t124 * t86 - t126 * t87) * pkin(3), -t30 * t32 - t31 * t33 + (t124 * t23 + t126 * t22 - t98 * t161) * pkin(3), t174, t189, t190, -t174, t191, 0, -t183 * t121 + t41 * t69 + t188, -t182 * t121 + t140 * t69 + t187, -t139 * t95 + t96 * t134 + (t182 + t5) * t41 + (-t183 - t6) * t140, -t135 * t96 + t182 * t6 - t183 * t5 + t2 * t95 - t50 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186 * qJD(3), -t112 + (-t100 + t150) * qJD(3), -t97 - t185, t100 * t31 + t137 * t30 + t114, 0, 0, 0, 0, 0, 0, -t134 - t169, t139 + t172, -t175 - t176, -t140 * t5 - t41 * t6 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t189, t190, -t174, t191, 0, t6 * t121 + t188, t5 * t121 + t187, 0, 0;];
tauc_reg = t1;
