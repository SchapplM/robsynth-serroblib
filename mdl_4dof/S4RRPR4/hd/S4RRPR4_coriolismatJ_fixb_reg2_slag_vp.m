% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:36
% EndTime: 2019-12-31 17:02:38
% DurationCPUTime: 0.84s
% Computational Cost: add. (1500->117), mult. (3044->155), div. (0->0), fcn. (2900->6), ass. (0->107)
t102 = sin(pkin(7));
t100 = t102 ^ 2;
t103 = cos(pkin(7));
t101 = t103 ^ 2;
t105 = sin(qJ(2));
t154 = t105 * pkin(1);
t96 = qJ(3) + t154;
t174 = (qJ(3) + t96) * (t101 / 0.2e1 + t100 / 0.2e1);
t156 = cos(qJ(4));
t124 = t156 * t103;
t104 = sin(qJ(4));
t141 = t104 * t102;
t109 = t124 - t141;
t133 = qJD(1) + qJD(2);
t125 = t156 * t102;
t140 = t104 * t103;
t79 = t125 + t140;
t166 = t133 * t79;
t172 = t109 * t166;
t106 = cos(qJ(2));
t153 = t106 * pkin(1);
t97 = -t103 * pkin(3) - pkin(2);
t88 = t97 - t153;
t126 = t97 / 0.2e1 + t88 / 0.2e1;
t170 = t126 * t79;
t163 = t79 ^ 2;
t75 = t109 ^ 2;
t43 = t75 - t163;
t169 = t133 * t43;
t50 = t75 + t163;
t168 = t133 * t50;
t167 = t133 * t109;
t92 = t101 + t100;
t165 = t133 * t92;
t157 = -pkin(6) - t96;
t155 = pkin(2) * t105;
t99 = t103 * pkin(6);
t74 = t103 * t96 + t99;
t46 = t104 * t74 - t157 * t125;
t152 = t46 * t79;
t47 = t157 * t141 + t156 * t74;
t151 = t47 * t109;
t122 = (-pkin(6) - qJ(3)) * t102;
t89 = t103 * qJ(3) + t99;
t53 = t104 * t89 - t156 * t122;
t150 = t53 * t79;
t54 = t104 * t122 + t156 * t89;
t149 = t54 * t109;
t61 = t79 * t153;
t62 = t109 * t153;
t17 = t109 * t62 + t61 * t79;
t49 = t50 * qJD(3);
t148 = t17 * qJD(2) + t49;
t120 = t92 * t106;
t73 = pkin(1) * t120;
t90 = t92 * qJD(3);
t147 = t73 * qJD(2) + t90;
t146 = pkin(1) * qJD(1);
t145 = pkin(1) * qJD(2);
t144 = qJD(1) * t88;
t143 = qJD(2) * t97;
t10 = t88 * t154 + t46 * t61 + t47 * t62;
t142 = t10 * qJD(1);
t12 = t151 + t152;
t139 = t12 * qJD(1);
t138 = t17 * qJD(1);
t63 = t92 * t96;
t48 = (-t155 + (t63 - t154) * t106) * pkin(1);
t137 = t48 * qJD(1);
t136 = t63 * qJD(1);
t135 = t73 * qJD(1);
t69 = t109 * qJD(4);
t134 = t79 * qJD(4);
t132 = t105 * t145;
t131 = t105 * t146;
t130 = t109 * t134;
t129 = t109 * t144;
t128 = t79 * t144;
t127 = t154 / 0.2e1;
t121 = pkin(1) * t133;
t117 = t109 * t131;
t116 = t79 * t131;
t115 = t102 * t131;
t114 = t105 * t121;
t15 = t149 + t150;
t6 = t127 + (-t53 / 0.2e1 - t46 / 0.2e1) * t79 - (t54 / 0.2e1 + t47 / 0.2e1) * t109;
t113 = t6 * qJD(1) - t15 * qJD(2);
t44 = t127 - t174;
t87 = t92 * qJ(3);
t112 = t44 * qJD(1) - t87 * qJD(2);
t108 = (-t140 / 0.2e1 - t125 / 0.2e1) * t153;
t18 = t108 - t170;
t111 = t18 * qJD(1) - t79 * t143;
t107 = (-t124 / 0.2e1 + t141 / 0.2e1) * t153;
t19 = -t109 * t126 + t107;
t110 = t19 * qJD(1) - t109 * t143;
t91 = t102 * t132;
t72 = t79 * qJD(3);
t68 = t109 * qJD(3);
t65 = t79 * t132;
t64 = t109 * t132;
t45 = t127 + t174;
t36 = t43 * qJD(4);
t21 = t108 + t170;
t20 = t107 + (t88 + t97) * t109 / 0.2e1;
t7 = t149 / 0.2e1 + t151 / 0.2e1 + t150 / 0.2e1 + t152 / 0.2e1 + t127;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t106 * t145, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * t132, t91, t147, t48 * qJD(2) + t63 * qJD(3), t130, t36, 0, -t130, 0, 0, t88 * t134 - t64, t88 * t69 + t65, t148, t10 * qJD(2) + t12 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t106 * t121, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * t114, t91 + t115, t135 + t147, t137 + t45 * qJD(3) + (qJ(3) * t120 - t155) * t145, t130, t36, 0, -t130, 0, 0, t21 * qJD(4) - t117 - t64, t20 * qJD(4) + t116 + t65, t138 + t148, t142 + (t97 * t154 + t61 * t53 + t62 * t54) * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t45 * qJD(2) + t136, 0, 0, 0, 0, 0, 0, 0, 0, t168, t7 * qJD(2) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t169, t69, -t172, -t134, 0, t21 * qJD(2) - t47 * qJD(4) + t128, t20 * qJD(2) + t46 * qJD(4) + t129, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t106 * t146, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t131, -t115, t90 - t135, -t44 * qJD(3) - t137, t130, t36, 0, -t130, 0, 0, -t18 * qJD(4) + t117, -t19 * qJD(4) - t116, t49 - t138, -t6 * qJD(3) - t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t87 * qJD(3), t130, t36, 0, -t130, 0, 0, t97 * t134, t97 * t69, t49, t15 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t112, 0, 0, 0, 0, 0, 0, 0, 0, t168, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t169, t69, -t172, -t134, 0, -t54 * qJD(4) - t111, t53 * qJD(4) - t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t44 * qJD(2) - t136, 0, 0, 0, 0, 0, 0, t134, t69, -t168, t6 * qJD(2) - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, t112, 0, 0, 0, 0, 0, 0, t134, t69, -t168, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t167, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t169, 0, t172, 0, 0, t18 * qJD(2) - t128 - t72, t19 * qJD(2) - t129 - t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t169, 0, t172, 0, 0, t111 - t72, t110 - t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t167, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
