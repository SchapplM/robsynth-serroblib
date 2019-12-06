% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:42
% EndTime: 2019-12-05 15:54:49
% DurationCPUTime: 1.84s
% Computational Cost: add. (2319->220), mult. (6038->302), div. (0->0), fcn. (4646->8), ass. (0->123)
t110 = cos(qJ(2));
t131 = t110 * qJD(1);
t118 = qJD(3) - t131;
t103 = sin(pkin(9));
t106 = sin(qJ(4));
t104 = cos(pkin(9));
t109 = cos(qJ(4));
t140 = t109 * t104;
t141 = t103 * t106;
t86 = -t140 + t141;
t113 = t86 * t110;
t135 = qJD(4) * t109;
t153 = pkin(6) + qJ(3);
t89 = t153 * t103;
t90 = t153 * t104;
t152 = -t89 * t135 + qJD(3) * t140 + (-qJD(3) * t103 - qJD(4) * t90) * t106 + qJD(1) * t113;
t49 = -t106 * t89 + t109 * t90;
t87 = t103 * t109 + t104 * t106;
t151 = -t49 * qJD(4) - t118 * t87;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t143 = qJD(2) * t87;
t128 = qJD(2) * t140;
t130 = qJD(2) * t141;
t75 = -t128 + t130;
t116 = t105 * t75 - t108 * t143;
t129 = qJD(4) * t141;
t92 = qJD(4) * t128;
t67 = qJD(2) * t129 - t92;
t80 = t87 * qJD(4);
t68 = qJD(2) * t80;
t112 = t116 * qJD(5) + t105 * t67 - t108 * t68;
t102 = qJD(4) + qJD(5);
t144 = t116 * t102;
t169 = t112 - t144;
t133 = qJD(5) * t108;
t134 = qJD(5) * t105;
t115 = -t105 * t68 - t108 * t67 - t75 * t133 - t134 * t143;
t43 = -t105 * t143 - t108 * t75;
t148 = t102 * t43;
t168 = t115 - t148;
t156 = t43 ^ 2;
t157 = t116 ^ 2;
t167 = -t156 + t157;
t98 = -pkin(3) * t104 - pkin(2);
t142 = qJD(2) * t98;
t81 = t118 + t142;
t47 = t75 * pkin(4) + t81;
t166 = t43 * t47;
t165 = -t80 * pkin(7) + t152;
t79 = -t104 * t135 + t129;
t164 = -t79 * pkin(7) - t151;
t163 = t116 * t47;
t155 = t43 * t116;
t137 = t103 ^ 2 + t104 ^ 2;
t107 = sin(qJ(2));
t132 = t107 * qJD(1);
t93 = qJD(2) * qJ(3) + t132;
t162 = t137 * t93;
t70 = t86 * t107;
t126 = pkin(6) * qJD(2) + t93;
t65 = t126 * t103;
t88 = (qJD(3) + t131) * qJD(2);
t150 = -t65 * t135 + t88 * t140;
t66 = t126 * t104;
t18 = (-qJD(4) * t66 - t103 * t88) * t106 + t150;
t12 = -t68 * pkin(7) + t18;
t114 = t87 * t88;
t33 = -t106 * t65 + t109 * t66;
t19 = -t33 * qJD(4) - t114;
t13 = t67 * pkin(7) + t19;
t27 = -pkin(7) * t75 + t33;
t124 = t105 * t13 - t27 * t134;
t146 = t106 * t66;
t32 = -t109 * t65 - t146;
t26 = -pkin(7) * t143 + t32;
t23 = qJD(4) * pkin(4) + t26;
t1 = (qJD(5) * t23 + t12) * t108 + t124;
t161 = t143 ^ 2;
t160 = pkin(4) * t143;
t48 = -t106 * t90 - t109 * t89;
t36 = -pkin(7) * t87 + t48;
t37 = -pkin(7) * t86 + t49;
t15 = t105 * t36 + t108 * t37;
t159 = t15 * qJD(5) + t165 * t105 + t164 * t108;
t14 = -t105 * t37 + t108 * t36;
t158 = -t14 * qJD(5) + t164 * t105 - t165 * t108;
t154 = t143 * t75;
t149 = qJD(2) * pkin(2);
t147 = t105 * t27;
t145 = t108 * t27;
t111 = qJD(2) ^ 2;
t139 = t111 * t107;
t138 = t111 * t110;
t136 = qJD(2) * t107;
t99 = qJD(2) * t132;
t50 = pkin(4) * t68 + t99;
t127 = -pkin(4) * t102 - t23;
t125 = -t105 * t12 + t108 * t13;
t122 = t137 * t88;
t120 = t137 * t110;
t119 = t137 * qJD(3);
t117 = pkin(4) * t80 - t132;
t8 = t105 * t23 + t145;
t69 = t87 * t107;
t34 = t105 * t70 - t108 * t69;
t35 = -t105 * t69 - t108 * t70;
t46 = -t105 * t86 + t108 * t87;
t2 = -qJD(5) * t8 + t125;
t91 = t118 - t149;
t72 = t75 ^ 2;
t58 = pkin(4) * t86 + t98;
t45 = t105 * t87 + t108 * t86;
t39 = qJD(4) * t70 - t110 * t143;
t38 = -qJD(2) * t113 - qJD(4) * t69;
t21 = t46 * qJD(5) - t105 * t79 + t108 * t80;
t20 = t105 * t80 + t108 * t79 + t86 * t133 + t87 * t134;
t10 = t108 * t26 - t147;
t9 = -t105 * t26 - t145;
t7 = t108 * t23 - t147;
t6 = -t35 * qJD(5) - t105 * t38 + t108 * t39;
t5 = t34 * qJD(5) + t105 * t39 + t108 * t38;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t138, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t139, t103 * t139, t137 * t138, t107 * t122 + (t107 * t91 + (-t132 + t162) * t110) * qJD(2), 0, 0, 0, 0, 0, 0, qJD(4) * t39 - t110 * t68 + t75 * t136, -qJD(4) * t38 + t110 * t67 + t136 * t143, -t143 * t39 - t38 * t75 - t67 * t69 + t68 * t70, -t18 * t70 - t19 * t69 + t32 * t39 + t33 * t38 + (t81 - t131) * t136, 0, 0, 0, 0, 0, 0, t102 * t6 + t110 * t112 - t136 * t43, -t102 * t5 - t110 * t115 - t116 * t136, t112 * t35 - t115 * t34 + t116 * t6 + t43 * t5, t1 * t35 - t110 * t50 + t47 * t136 + t2 * t34 + t5 * t8 + t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 + (-qJD(1) * t120 + t119) * qJD(2), t93 * t119 + qJ(3) * t122 + (-t93 * t120 + (-t91 - t149) * t107) * qJD(1), -t143 * t79 - t67 * t87, -t143 * t80 + t67 * t86 - t68 * t87 + t75 * t79, -t79 * qJD(4), t68 * t86 + t75 * t80, -t80 * qJD(4), 0, t98 * t68 + t81 * t80 + t151 * qJD(4) + (qJD(2) * t86 - t75) * t132, -t152 * qJD(4) - t98 * t67 - t81 * t79, -t143 * t151 - t152 * t75 - t18 * t86 - t19 * t87 + t32 * t79 - t33 * t80 + t48 * t67 - t49 * t68, t18 * t49 + t19 * t48 + t152 * t33 + t151 * t32 + (-t81 + t142) * t132, t115 * t46 + t116 * t20, t112 * t46 - t115 * t45 + t116 * t21 - t20 * t43, -t20 * t102, -t112 * t45 - t21 * t43, -t21 * t102, 0, -t159 * t102 - t112 * t58 - t117 * t43 + t21 * t47 + t45 * t50, t158 * t102 + t115 * t58 - t116 * t117 - t20 * t47 + t46 * t50, -t1 * t45 + t112 * t15 - t115 * t14 - t116 * t159 - t158 * t43 - t2 * t46 + t20 * t7 - t21 * t8, t1 * t15 + t117 * t47 + t14 * t2 - t158 * t8 - t159 * t7 + t50 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t111, -qJD(2) * t162 + t99, 0, 0, 0, 0, 0, 0, 0.2e1 * t143 * qJD(4), t92 + (-t75 - t130) * qJD(4), -t72 - t161, t143 * t32 + t33 * t75 + t99, 0, 0, 0, 0, 0, 0, -t112 - t144, t115 + t148, -t156 - t157, -t116 * t7 - t43 * t8 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t72 + t161, t92 + (t75 - t130) * qJD(4), -t154, 0, 0, -t143 * t81 - t114, t88 * t141 + t81 * t75 + (t32 + t146) * qJD(4) - t150, 0, 0, t155, t167, t168, -t155, t169, 0, t43 * t160 - t102 * t9 + t163 + (t105 * t127 - t145) * qJD(5) + t125, t116 * t160 + t10 * t102 - t166 + (qJD(5) * t127 - t12) * t108 - t124, -t10 * t43 - t116 * t8 + t43 * t7 - t116 * t9 + (t105 * t112 - t108 * t115 + (-t105 * t116 + t108 * t43) * qJD(5)) * pkin(4), -t10 * t8 - t7 * t9 + (t1 * t105 + t108 * t2 - t47 * t143 + (-t105 * t7 + t108 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t167, t168, -t155, t169, 0, t8 * t102 + t163 + t2, t7 * t102 - t1 - t166, 0, 0;];
tauc_reg = t3;
