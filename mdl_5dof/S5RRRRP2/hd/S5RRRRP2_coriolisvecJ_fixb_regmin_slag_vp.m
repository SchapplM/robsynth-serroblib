% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:12
% EndTime: 2019-12-05 18:48:17
% DurationCPUTime: 1.10s
% Computational Cost: add. (1873->187), mult. (3202->245), div. (0->0), fcn. (2038->6), ass. (0->146)
t112 = cos(qJ(4));
t109 = sin(qJ(4));
t113 = cos(qJ(3));
t106 = qJD(1) + qJD(2);
t111 = sin(qJ(2));
t171 = pkin(1) * qJD(1);
t151 = t111 * t171;
t193 = pkin(7) + pkin(8);
t138 = t193 * t106 + t151;
t53 = t138 * t113;
t46 = t109 * t53;
t110 = sin(qJ(3));
t52 = t138 * t110;
t49 = qJD(3) * pkin(3) - t52;
t134 = t112 * t49 - t46;
t77 = t109 * t113 + t112 * t110;
t61 = t77 * t106;
t168 = t61 * qJ(5);
t194 = t168 - t134;
t156 = qJD(4) * t112;
t159 = qJD(3) * t113;
t192 = -t112 * t159 - t113 * t156;
t114 = cos(qJ(2));
t150 = t114 * t171;
t157 = qJD(4) * t109;
t163 = t112 * t113;
t165 = t109 * t110;
t76 = -t163 + t165;
t146 = qJD(3) * t193;
t78 = t110 * t146;
t79 = t113 * t146;
t89 = t193 * t110;
t103 = t113 * pkin(8);
t90 = t113 * pkin(7) + t103;
t191 = t109 * t79 + t112 * t78 - t76 * t150 + t89 * t156 + t90 * t157;
t126 = t109 * t89 - t112 * t90;
t190 = t126 * qJD(4) + t109 * t78 - t112 * t79 + t77 * t150;
t130 = qJD(3) * t138;
t170 = pkin(1) * qJD(2);
t148 = qJD(1) * t170;
t131 = t114 * t148;
t32 = -t110 * t130 + t113 * t131;
t189 = (qJD(4) * t49 + t32) * t112;
t105 = qJD(3) + qJD(4);
t188 = t61 ^ 2;
t186 = t76 * pkin(4);
t96 = t111 * pkin(1) + pkin(7);
t185 = -pkin(8) - t96;
t184 = t114 * pkin(1);
t147 = t106 * t165;
t59 = -t106 * t163 + t147;
t99 = -t113 * pkin(3) - pkin(2);
t63 = t99 * t106 - t150;
t28 = t59 * pkin(4) + qJD(5) + t63;
t183 = t28 * t61;
t182 = t61 * t59;
t181 = t63 * t61;
t38 = t105 * t77;
t176 = -t38 * qJ(5) - t76 * qJD(5);
t180 = t176 - t191;
t129 = t105 * t165;
t37 = t129 + t192;
t125 = t37 * qJ(5) - t77 * qJD(5);
t179 = t125 + t190;
t13 = t105 * pkin(4) - t194;
t178 = t13 + t194;
t160 = qJD(3) * t110;
t145 = t106 * t160;
t94 = t111 * t148;
t65 = pkin(3) * t145 + t94;
t177 = t63 * t38 + t65 * t76;
t175 = -t63 * t37 + t65 * t77;
t174 = -t112 * t52 - t46;
t83 = -t106 * pkin(2) - t150;
t173 = t110 * t94 + t83 * t159;
t172 = t192 * t106;
t48 = t112 * t53;
t169 = t59 * qJ(5);
t167 = t77 * qJ(5);
t166 = t106 * t110;
t164 = t111 * t113;
t115 = qJD(3) ^ 2;
t162 = t115 * t110;
t102 = t115 * t113;
t161 = t110 ^ 2 - t113 ^ 2;
t158 = qJD(3) * t114;
t155 = -qJD(1) - t106;
t154 = -qJD(2) + t106;
t153 = pkin(3) * t166;
t152 = t114 * t170;
t101 = t111 * t170;
t100 = pkin(3) * t160;
t128 = -t109 * t49 - t48;
t15 = -t128 - t169;
t33 = -t110 * t131 - t113 * t130;
t135 = t109 * t33 - t53 * t157;
t27 = t38 * t106;
t4 = -t27 * qJ(5) - t59 * qJD(5) + t135 + t189;
t136 = -t109 * t32 + t112 * t33;
t118 = t128 * qJD(4) + t136;
t26 = t106 * t129 + t172;
t5 = t26 * qJ(5) - t61 * qJD(5) + t118;
t149 = t13 * t37 - t15 * t38 - t4 * t76 - t5 * t77;
t144 = t106 * t159;
t143 = t110 * t158;
t140 = -pkin(3) * t105 - t49;
t139 = t38 * pkin(4) + t100;
t137 = qJD(3) * t185;
t133 = t109 * t52 - t48;
t20 = t27 * pkin(4) + t65;
t73 = t185 * t110;
t74 = t113 * t96 + t103;
t127 = -t109 * t73 - t112 * t74;
t86 = t99 - t184;
t124 = t63 * t59 - t135;
t50 = t110 * t137 + t113 * t152;
t51 = -t110 * t152 + t113 * t137;
t123 = t109 * t51 + t112 * t50 + t73 * t156 - t74 * t157;
t121 = -t151 + t100;
t120 = -t106 * t83 - t131;
t119 = -t111 * t166 + t113 * t158;
t117 = t127 * qJD(4) - t109 * t50 + t112 * t51;
t104 = t106 ^ 2;
t98 = -pkin(2) - t184;
t97 = t112 * pkin(3) + pkin(4);
t81 = 0.2e1 * t110 * t144;
t80 = t101 + t100;
t72 = t76 * qJ(5);
t66 = t83 * t160;
t58 = -0.2e1 * t161 * t106 * qJD(3);
t57 = t59 ^ 2;
t35 = t38 * t105;
t34 = t37 * t105;
t30 = -t126 - t72;
t29 = -t109 * t90 - t112 * t89 - t167;
t23 = -t127 - t72;
t22 = -t109 * t74 + t112 * t73 - t167;
t21 = -t57 + t188;
t18 = -t172 + (-t147 + t59) * t105;
t17 = -t168 + t174;
t16 = t133 + t169;
t10 = -t26 * t77 - t61 * t37;
t7 = t117 + t125;
t6 = t123 + t176;
t3 = t26 * t76 - t77 * t27 + t37 * t59 - t61 * t38;
t1 = [0, 0, 0, 0, -t106 * t101 - t94, t155 * t152, t81, t58, t102, -t162, 0, t98 * t145 - t96 * t102 + t66 + (t155 * t164 - t143) * t170, -t119 * t170 + t98 * t144 + t96 * t162 + t173, t10, t3, -t34, -t35, 0, t117 * t105 + t86 * t27 + t80 * t59 + t177, -t123 * t105 - t86 * t26 + t80 * t61 + t175, t22 * t26 - t23 * t27 - t6 * t59 - t7 * t61 + t149, t4 * t23 + t15 * t6 + t5 * t22 + t13 * t7 + t20 * (t86 + t186) + t28 * (t101 + t139); 0, 0, 0, 0, t106 * t151 - t94, t154 * t150, t81, t58, t102, -t162, 0, -pkin(2) * t145 - pkin(7) * t102 + t66 + (t154 * t164 + t143) * t171, -pkin(2) * t144 + pkin(7) * t162 + t119 * t171 + t173, t10, t3, -t34, -t35, 0, t105 * t190 + t121 * t59 + t99 * t27 + t177, t105 * t191 + t121 * t61 - t99 * t26 + t175, -t179 * t61 - t180 * t59 + t29 * t26 - t30 * t27 + t149, t4 * t30 + t5 * t29 + t20 * (t99 + t186) + (t139 - t151) * t28 + t180 * t15 + t179 * t13; 0, 0, 0, 0, 0, 0, -t110 * t104 * t113, t161 * t104, 0, 0, 0, t120 * t110, t120 * t113, t182, t21, t18, 0, 0, -t59 * t153 - t181 - t133 * t105 + (t140 * t109 - t48) * qJD(4) + t136, -t61 * t153 + t174 * t105 + (t140 * qJD(4) - t32) * t112 + t124, t97 * t26 + (t15 + t16) * t61 + (-t13 + t17) * t59 + (-t109 * t27 + (t109 * t61 - t112 * t59) * qJD(4)) * pkin(3), -pkin(4) * t183 - t13 * t16 - t15 * t17 + t5 * t97 + (-t28 * t166 + t4 * t109 + (-t109 * t13 + t112 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t21, t18, 0, 0, -t128 * t105 + t118 - t181, t134 * t105 + t124 - t189, pkin(4) * t26 - t178 * t59, t178 * t15 + (t5 - t183) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 - t188, t13 * t61 + t15 * t59 + t20;];
tauc_reg = t1;
