% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:22
% DurationCPUTime: 1.83s
% Computational Cost: add. (2107->270), mult. (5054->373), div. (0->0), fcn. (3024->6), ass. (0->155)
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t143 = qJD(4) * t102;
t144 = qJD(4) * t100;
t103 = cos(qJ(3));
t101 = sin(qJ(3));
t142 = t101 * qJD(2);
t90 = sin(pkin(8)) * pkin(1) + pkin(6);
t84 = t90 * qJD(1);
t59 = t103 * t84 + t142;
t51 = qJD(3) * pkin(7) + t59;
t188 = qJD(2) * t103 - t101 * t84;
t52 = t188 * qJD(3);
t91 = -cos(pkin(8)) * pkin(1) - pkin(2);
t73 = -pkin(3) * t103 - pkin(7) * t101 + t91;
t55 = t73 * qJD(1);
t120 = pkin(3) * t101 - pkin(7) * t103;
t83 = t120 * qJD(3);
t72 = qJD(1) * t83;
t127 = -t100 * t72 - t102 * t52 - t55 * t143 + t51 * t144;
t26 = -t100 * t51 + t102 * t55;
t149 = qJD(1) * t103;
t89 = -qJD(4) + t149;
t192 = t26 * t89 - t127;
t27 = t100 * t55 + t102 * t51;
t7 = -qJD(4) * t27 - t100 * t52 + t102 * t72;
t191 = t27 * t89 - t7;
t141 = t102 * qJD(3);
t150 = qJD(1) * t101;
t77 = t100 * t150 - t141;
t180 = t77 * t89;
t133 = t101 * t144;
t134 = t103 * t141;
t138 = qJD(3) * qJD(4);
t46 = -t102 * t138 + (t133 - t134) * qJD(1);
t190 = -t46 + t180;
t147 = qJD(3) * t100;
t79 = t102 * t150 + t147;
t178 = t79 * t89;
t132 = t101 * t143;
t145 = qJD(3) * t103;
t109 = t100 * t145 + t132;
t47 = qJD(1) * t109 + t100 * t138;
t189 = t47 - t178;
t96 = t101 ^ 2;
t115 = qJD(1) * t96 - t103 * t89;
t187 = -t115 * t141 - t89 * t133;
t186 = t79 ^ 2;
t185 = pkin(4) * t77;
t184 = pkin(4) * t100;
t20 = -qJ(5) * t77 + t27;
t183 = t20 * t89;
t179 = t79 * t77;
t177 = -qJ(5) - pkin(7);
t19 = -qJ(5) * t79 + t26;
t14 = -pkin(4) * t89 + t19;
t176 = t14 - t19;
t153 = t102 * t103;
t112 = pkin(4) * t101 - qJ(5) * t153;
t128 = qJD(4) * t177;
t82 = t120 * qJD(1);
t35 = -t100 * t188 + t102 * t82;
t175 = qJD(1) * t112 + t100 * qJD(5) - t102 * t128 + t35;
t140 = t102 * qJD(5);
t36 = t100 * t82 + t102 * t188;
t174 = -t140 + t36 + (-qJ(5) * t149 - t128) * t100;
t154 = t101 * t102;
t173 = -t77 * t134 - t47 * t154;
t172 = t100 * t83 + t73 * t143;
t146 = qJD(3) * t101;
t167 = t100 * t90;
t171 = t102 * t83 + t146 * t167;
t81 = t90 * t153;
t39 = t100 * t73 + t81;
t170 = -t103 ^ 2 + t96;
t50 = -qJD(3) * pkin(3) - t188;
t169 = t100 * t50;
t168 = t100 * t89;
t166 = t101 * t77;
t165 = t102 * t50;
t164 = t102 * t89;
t163 = t103 * t47;
t53 = t59 * qJD(3);
t28 = t47 * pkin(4) + t53;
t162 = t28 * t100;
t161 = t28 * t102;
t160 = t53 * t100;
t159 = t53 * t101;
t158 = t53 * t102;
t157 = qJ(5) * t101;
t85 = qJD(1) * t91;
t156 = qJD(4) * t77;
t155 = t100 * t103;
t104 = qJD(3) ^ 2;
t152 = t104 * t101;
t151 = t104 * t103;
t139 = qJD(1) * qJD(3);
t137 = t79 * t145;
t105 = qJD(1) ^ 2;
t136 = t101 * t105 * t103;
t135 = t89 * t150;
t92 = t101 * t139;
t130 = -qJD(5) - t185;
t129 = t46 * t103 + t79 * t146;
t126 = -t46 + t156;
t124 = t89 * t132;
t123 = pkin(4) * t92;
t122 = t79 * t132;
t121 = t103 * t92;
t119 = -t100 * t20 - t102 * t14;
t118 = t100 * t14 - t102 * t20;
t117 = -t100 * t27 - t102 * t26;
t116 = t100 * t26 - t102 * t27;
t114 = 0.2e1 * qJD(3) * t85;
t111 = t47 * qJ(5) + t127;
t110 = t115 * t100;
t108 = qJD(4) * t117 - t7 * t100 - t102 * t127;
t107 = t159 + t52 * t103 + (-t101 * t59 - t103 * t188) * qJD(3);
t106 = t46 * qJ(5) + t7;
t94 = -pkin(4) * t102 - pkin(3);
t87 = t177 * t102;
t86 = t177 * t100;
t75 = t77 ^ 2;
t67 = (t90 + t184) * t101;
t63 = t102 * t73;
t54 = (-t89 - t149) * t146;
t42 = t142 + (qJD(1) * t184 + t84) * t103;
t41 = pkin(4) * t109 + t145 * t90;
t38 = -t155 * t90 + t63;
t37 = -t130 + t50;
t34 = -t100 * t157 + t39;
t33 = -t75 + t186;
t32 = -qJ(5) * t154 + t63 + (-pkin(4) - t167) * t103;
t31 = -t178 - t47;
t30 = -t46 - t180;
t24 = -t89 * t143 + (t89 * t153 + (-t79 + t147) * t101) * qJD(1);
t23 = t89 * t144 + (-t89 * t155 + (t77 + t141) * t101) * qJD(1);
t22 = -t39 * qJD(4) + t171;
t21 = (-t101 * t141 - t103 * t144) * t90 + t172;
t18 = -t47 * t102 - t168 * t77;
t17 = -t46 * t100 - t164 * t79;
t16 = t77 * t132 + (t101 * t47 + t145 * t77) * t100;
t15 = t79 * t134 + (-t46 * t102 - t144 * t79) * t101;
t13 = t124 + t163 + (-t110 - t166) * qJD(3);
t12 = t124 - t163 + (-t110 + t166) * qJD(3);
t11 = t129 - t187;
t10 = t129 + t187;
t9 = (-qJ(5) * qJD(4) - qJD(3) * t90) * t154 + (-qJD(5) * t101 + (-qJ(5) * qJD(3) - qJD(4) * t90) * t103) * t100 + t172;
t8 = -t101 * t140 + t112 * qJD(3) + (-t81 + (-t73 + t157) * t100) * qJD(4) + t171;
t5 = -t189 * t100 + t190 * t102;
t4 = -qJD(5) * t77 - t111;
t3 = -t122 + (-t137 + (t46 + t156) * t101) * t100 + t173;
t2 = t122 + (t101 * t126 + t137) * t100 + t173;
t1 = -t79 * qJD(5) + t106 + t123;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, -0.2e1 * t170 * t139, t151, -0.2e1 * t121, -t152, 0, t101 * t114 - t151 * t90, t103 * t114 + t152 * t90, t107, t107 * t90, t15, t3, t11, t16, t13, t54, -t22 * t89 + (-t7 + (t77 * t90 + t169) * qJD(3)) * t103 + (t50 * t143 + t160 + t47 * t90 + (qJD(1) * t38 + t26) * qJD(3)) * t101, t21 * t89 + (-t127 + (t79 * t90 + t165) * qJD(3)) * t103 + (-t50 * t144 + t158 - t46 * t90 + (-qJD(1) * t39 - t27) * qJD(3)) * t101, -t21 * t77 - t22 * t79 + t38 * t46 - t39 * t47 + t117 * t145 + (qJD(4) * t116 + t100 * t127 - t102 * t7) * t101, t27 * t21 + t26 * t22 + t7 * t38 - t127 * t39 + (t145 * t50 + t159) * t90, t15, t3, t11, t16, t13, t54, t41 * t77 + t67 * t47 - t8 * t89 + (t147 * t37 - t1) * t103 + (t37 * t143 + t162 + (qJD(1) * t32 + t14) * qJD(3)) * t101, t41 * t79 - t67 * t46 + t9 * t89 + (t141 * t37 + t4) * t103 + (-t37 * t144 + t161 + (-qJD(1) * t34 - t20) * qJD(3)) * t101, t32 * t46 - t34 * t47 - t9 * t77 - t8 * t79 + t119 * t145 + (qJD(4) * t118 - t1 * t102 - t100 * t4) * t101, t1 * t32 + t14 * t8 + t20 * t9 + t28 * t67 + t34 * t4 + t37 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151, 0, t52 * t101 - t53 * t103 + (-t101 * t188 + t103 * t59) * qJD(3), 0, 0, 0, 0, 0, 0, t12, t10, t2, (-qJD(3) * t116 - t53) * t103 + (qJD(3) * t50 + t108) * t101, 0, 0, 0, 0, 0, 0, t12, t10, t2, (-qJD(3) * t118 - t28) * t103 + (qJD(3) * t37 + qJD(4) * t119 - t1 * t100 + t4 * t102) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t170 * t105, 0, t136, 0, 0, -t85 * t150, -t85 * t149, 0, 0, t17, t5, t24, t18, t23, t135, -pkin(3) * t47 - t158 + t35 * t89 - t59 * t77 + (pkin(7) * t164 + t169) * qJD(4) + (-t101 * t26 + (-pkin(7) * t146 - t103 * t50) * t100) * qJD(1), pkin(3) * t46 + t160 - t36 * t89 - t59 * t79 + (-pkin(7) * t168 + t165) * qJD(4) + (-t50 * t153 + (-pkin(7) * t141 + t27) * t101) * qJD(1), t35 * t79 + t36 * t77 + ((qJD(4) * t79 - t47) * pkin(7) + t192) * t102 + (pkin(7) * t126 + t191) * t100, -t53 * pkin(3) + pkin(7) * t108 - t26 * t35 - t27 * t36 - t50 * t59, t17, t5, t24, t18, t23, t135, -t161 - t42 * t77 + t47 * t94 + t175 * t89 + (t37 + t185) * t144 + (-t37 * t155 + (qJD(3) * t86 - t14) * t101) * qJD(1), t162 - t42 * t79 - t46 * t94 - t174 * t89 + (t102 * t37 + t79 * t184) * qJD(4) + (-t37 * t153 + (qJD(3) * t87 + t20) * t101) * qJD(1), t46 * t86 + t47 * t87 + t175 * t79 + t174 * t77 + (t14 * t89 + t4) * t102 + (-t1 + t183) * t100, t1 * t86 + t28 * t94 - t4 * t87 + (pkin(4) * t144 - t42) * t37 - t174 * t20 - t175 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t33, t30, -t179, t31, t92, -t50 * t79 - t191, t50 * t77 - t192, 0, 0, t179, t33, t30, -t179, t31, t92, 0.2e1 * t123 - t183 + (t130 - t37) * t79 + t106, -t186 * pkin(4) - t19 * t89 + (qJD(5) + t37) * t77 + t111, t46 * pkin(4) - t176 * t77, t176 * t20 + (-t37 * t79 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t190, -t75 - t186, t14 * t79 + t20 * t77 + t28;];
tauc_reg = t6;
