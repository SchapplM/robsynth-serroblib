% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:47
% EndTime: 2019-12-31 17:42:49
% DurationCPUTime: 1.27s
% Computational Cost: add. (2059->212), mult. (3764->277), div. (0->0), fcn. (2812->14), ass. (0->142)
t102 = sin(pkin(8));
t104 = cos(pkin(8));
t131 = g(1) * t104 + g(2) * t102;
t100 = qJ(2) + qJ(3);
t90 = pkin(9) + t100;
t81 = sin(t90);
t187 = t131 * t81;
t91 = sin(t100);
t92 = cos(t100);
t181 = -g(3) * t92 + t131 * t91;
t82 = cos(t90);
t176 = g(3) * t82;
t101 = sin(pkin(9));
t103 = cos(pkin(9));
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t107 = sin(qJ(2));
t110 = cos(qJ(2));
t149 = qJD(1) * qJD(2);
t127 = -t107 * qJDD(1) - t110 * t149;
t76 = qJD(2) * pkin(2) + qJD(1) * t110;
t118 = qJD(3) * t76 - t127;
t151 = qJD(1) * t107;
t138 = qJD(3) * t151;
t89 = t110 * qJDD(1);
t54 = qJDD(2) * pkin(2) - t107 * t149 + t89;
t46 = t109 * t54;
t22 = -t118 * t106 - t109 * t138 + t46;
t96 = qJDD(2) + qJDD(3);
t15 = t96 * pkin(3) + t22;
t180 = t106 * t54 + t118 * t109;
t75 = t106 * t138;
t21 = t180 - t75;
t137 = t101 * t21 - t103 * t15;
t5 = -pkin(4) * t96 + t137;
t185 = t5 + t176;
t157 = t101 * t106;
t162 = pkin(2) * qJD(3);
t58 = t106 * t110 + t107 * t109;
t52 = t58 * qJD(1);
t57 = -t106 * t107 + t109 * t110;
t53 = t57 * qJD(1);
t165 = t101 * t52 - t103 * t53 + (t103 * t109 - t157) * t162;
t97 = qJD(2) + qJD(3);
t183 = t165 * t97;
t111 = qJD(5) ^ 2;
t156 = t103 * t106;
t166 = -t101 * t53 - t103 * t52 + (t101 * t109 + t156) * t162;
t141 = t166 * t97;
t87 = pkin(2) * t109 + pkin(3);
t49 = -pkin(2) * t157 + t103 * t87;
t44 = -pkin(4) - t49;
t50 = pkin(2) * t156 + t101 * t87;
t45 = pkin(7) + t50;
t182 = t111 * t45 + t44 * t96 + t141;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t39 = -t106 * t151 + t109 * t76;
t36 = pkin(3) * t97 + t39;
t40 = t106 * t76 + t109 * t151;
t37 = t103 * t40;
t24 = t101 * t36 + t37;
t20 = pkin(7) * t97 + t24;
t11 = qJD(4) * t108 - t105 * t20;
t12 = qJD(4) * t105 + t108 * t20;
t154 = t11 * qJD(5);
t8 = t101 * t15 + t103 * t21;
t6 = pkin(7) * t96 + t8;
t2 = t105 * qJDD(4) + t108 * t6 + t154;
t153 = t12 * qJD(5);
t88 = t108 * qJDD(4);
t3 = -t105 * t6 - t153 + t88;
t116 = -t3 * t105 + t2 * t108 + (-t105 * t12 - t108 * t11) * qJD(5);
t139 = -g(1) * t102 + g(2) * t104;
t179 = pkin(2) * t96;
t178 = pkin(3) * t91;
t177 = pkin(4) * t81;
t79 = g(3) * t81;
t174 = pkin(3) * t101;
t173 = pkin(3) * t103;
t32 = t97 * t57;
t33 = t97 * t58;
t10 = -t101 * t33 + t103 * t32;
t169 = t10 * t97;
t25 = t101 * t39 + t37;
t168 = t25 * t97;
t161 = t101 * t40;
t26 = t103 * t39 - t161;
t167 = t26 * t97;
t98 = t105 ^ 2;
t99 = t108 ^ 2;
t164 = t98 - t99;
t163 = t98 + t99;
t160 = t102 * t82;
t159 = t104 * t82;
t155 = t105 * t108;
t152 = qJDD(1) - g(3);
t150 = qJD(5) * t108;
t23 = t103 * t36 - t161;
t19 = -pkin(4) * t97 - t23;
t148 = t185 * t105 + t19 * t150;
t147 = t19 * qJD(5) * t105 + t108 * t187;
t146 = g(1) * t159 + g(2) * t160 + t79;
t86 = pkin(3) * t92;
t145 = t82 * pkin(4) + t81 * pkin(7) + t86;
t95 = t97 ^ 2;
t144 = t95 * t155;
t60 = -pkin(2) * t107 - t178;
t142 = t60 - t177;
t140 = t163 * t96;
t136 = g(3) * t91 + t131 * t92 + t75;
t135 = t105 * t97 * t150;
t134 = t146 - t8;
t133 = -t177 - t178;
t30 = t101 * t58 - t103 * t57;
t9 = t101 * t32 + t103 * t33;
t132 = -t30 * t96 - t9 * t97;
t129 = t105 * t11 - t108 * t12;
t31 = t101 * t57 + t103 * t58;
t126 = t111 * t31 - t132;
t83 = pkin(7) + t174;
t84 = -pkin(4) - t173;
t125 = t111 * t83 + t84 * t96 - t168;
t124 = -qJD(5) * t20 + t139;
t123 = -t137 - t176 + t187;
t122 = -qJDD(5) * t31 + (t30 * t97 - t10) * qJD(5);
t121 = -qJDD(5) * t83 + (t84 * t97 + t26) * qJD(5);
t120 = -qJDD(5) * t45 + (t44 * t97 - t165) * qJD(5);
t119 = -g(3) * t110 + t131 * t107;
t115 = (-pkin(2) * t97 - t76) * qJD(3) + t127;
t114 = -t146 + t116;
t113 = -qJD(5) * qJD(4) + t131 * t82 - t19 * t97 - t6 + t79;
t112 = qJD(2) ^ 2;
t94 = t110 * pkin(2);
t68 = qJDD(5) * t108 - t105 * t111;
t67 = qJDD(5) * t105 + t108 * t111;
t62 = pkin(7) * t159;
t61 = pkin(7) * t160;
t43 = t96 * t99 - 0.2e1 * t135;
t42 = t96 * t98 + 0.2e1 * t135;
t34 = -0.2e1 * t164 * t97 * qJD(5) + 0.2e1 * t96 * t155;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, 0, 0, 0, 0, 0, qJDD(2) * t110 - t107 * t112, -qJDD(2) * t107 - t110 * t112, 0, -g(3) + (t107 ^ 2 + t110 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t33 * t97 + t57 * t96, -t32 * t97 - t58 * t96, 0, t21 * t58 + t22 * t57 + t32 * t40 - t33 * t39 - g(3), 0, 0, 0, 0, 0, 0, t132, -t31 * t96 - t169, 0, t10 * t24 + t137 * t30 - t23 * t9 + t31 * t8 - g(3), 0, 0, 0, 0, 0, 0, t122 * t105 - t126 * t108, t126 * t105 + t122 * t108, t31 * t140 + t163 * t169, -t129 * t10 + t116 * t31 + t19 * t9 + t5 * t30 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t89 + t119, -t152 * t107 + t131 * t110, 0, 0, 0, 0, 0, 0, 0, t96, t52 * t97 + t46 + (-t138 + t179) * t109 + t115 * t106 + t181, t53 * t97 + (-t54 - t179) * t106 + t115 * t109 + t136, 0, t39 * t52 - t40 * t53 + (t106 * t21 + t109 * t22 + (-t106 * t39 + t109 * t40) * qJD(3) + t119) * pkin(2), 0, 0, 0, 0, 0, t96, t49 * t96 + t123 - t141, -t50 * t96 + t134 - t183, 0, t8 * t50 - t137 * t49 - g(3) * (t86 + t94) - t131 * t60 + t165 * t24 - t166 * t23, t42, t34, t67, t43, t68, 0, t120 * t105 + (-t185 - t182) * t108 + t147, t120 * t108 + (-t187 + t182) * t105 + t148, t45 * t140 + t163 * t183 + t114, t5 * t44 - g(1) * (t142 * t104 + t62) - g(2) * (t142 * t102 + t61) - g(3) * (t94 + t145) + t166 * t19 + ((t2 - t154) * t45 + t165 * t12) * t108 + ((-t3 - t153) * t45 - t165 * t11) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t40 * t97 + t181 + t22, t39 * t97 + t136 - t180, 0, 0, 0, 0, 0, 0, 0, t96, t96 * t173 + t123 + t168, -t96 * t174 + t134 + t167, 0, t23 * t25 - t24 * t26 + (t101 * t8 - t103 * t137 + t181) * pkin(3), t42, t34, t67, t43, t68, 0, t121 * t105 + (-t125 - t185) * t108 + t147, t121 * t108 + (-t187 + t125) * t105 + t148, t83 * t140 - t163 * t167 + t114, t5 * t84 - t19 * t25 - g(1) * (t133 * t104 + t62) - g(2) * (t133 * t102 + t61) - g(3) * t145 + t129 * t26 + t116 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + t139, 0, 0, 0, 0, 0, 0, t68, -t67, 0, -t129 * qJD(5) + t2 * t105 + t3 * t108 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t164 * t95, t105 * t96, t144, t108 * t96, qJDD(5), t113 * t105 + t124 * t108 + t153 + t88, t154 + (-qJDD(4) - t124) * t105 + t113 * t108, 0, 0;];
tau_reg = t1;
