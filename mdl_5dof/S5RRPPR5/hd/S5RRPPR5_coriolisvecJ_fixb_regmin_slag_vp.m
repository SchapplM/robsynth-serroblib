% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:30:01
% DurationCPUTime: 1.04s
% Computational Cost: add. (1197->193), mult. (3180->260), div. (0->0), fcn. (2221->6), ass. (0->123)
t100 = qJD(2) - qJD(5);
t107 = sin(qJ(5));
t109 = cos(qJ(5));
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t110 = cos(qJ(2));
t147 = t106 * t110;
t135 = qJD(1) * t147;
t108 = sin(qJ(2));
t139 = t108 * qJD(1);
t70 = t105 * t139 - t135;
t83 = t105 * t110 + t106 * t108;
t73 = t83 * qJD(1);
t165 = t107 * t70 + t109 * t73;
t150 = t165 * t100;
t72 = t83 * qJD(2);
t63 = qJD(1) * t72;
t138 = qJD(1) * qJD(2);
t134 = t108 * t138;
t92 = t105 * t134;
t64 = qJD(2) * t135 - t92;
t3 = qJD(5) * t165 + t107 * t64 - t109 * t63;
t174 = t3 + t150;
t136 = t110 * pkin(2) + pkin(1);
t124 = t136 * qJD(1);
t89 = qJD(3) - t124;
t173 = -t73 * qJ(4) + t89;
t121 = t107 * t73 - t109 * t70;
t130 = t121 * qJD(5) - t107 * t63 - t109 * t64;
t149 = t121 * t100;
t172 = t130 + t149;
t171 = -t121 ^ 2 + t165 ^ 2;
t162 = -pkin(3) - pkin(4);
t10 = t162 * t70 - t173;
t101 = qJD(2) * qJD(4);
t157 = -qJ(3) - pkin(6);
t131 = qJD(2) * t157;
t66 = t110 * qJD(3) + t108 * t131;
t53 = t66 * qJD(1);
t67 = -t108 * qJD(3) + t110 * t131;
t54 = t67 * qJD(1);
t21 = t105 * t54 + t106 * t53;
t18 = t101 + t21;
t7 = t63 * pkin(7) + t18;
t20 = t105 * t53 - t106 * t54;
t8 = -t64 * pkin(7) + t20;
t170 = t10 * t165 + t107 * t7 - t109 * t8;
t168 = -0.2e1 * t138;
t69 = t73 ^ 2;
t167 = -t70 ^ 2 - t69;
t166 = t165 * t121;
t91 = t157 * t110;
t88 = qJD(1) * t91;
t156 = t105 * t88;
t90 = t157 * t108;
t87 = qJD(1) * t90;
t44 = t106 * t87 + t156;
t142 = qJD(4) - t44;
t164 = qJD(5) + t100;
t163 = t10 * t121 - t107 * t8 - t109 * t7;
t161 = t63 * pkin(3);
t160 = t70 * pkin(7);
t159 = t73 * pkin(7);
t45 = -t105 * t91 - t106 * t90;
t158 = t20 * t45;
t28 = t105 * t67 + t106 * t66;
t155 = t106 * t88;
t81 = qJD(2) * pkin(2) + t87;
t40 = t105 * t81 - t155;
t46 = t105 * t90 - t106 * t91;
t112 = qJD(1) ^ 2;
t146 = t110 * t112;
t111 = qJD(2) ^ 2;
t145 = t111 * t108;
t144 = t111 * t110;
t143 = -t159 + t142;
t141 = t108 ^ 2 - t110 ^ 2;
t140 = qJD(2) * t108;
t37 = qJD(2) * qJ(4) + t40;
t137 = pkin(2) * t140;
t99 = -t106 * pkin(2) - pkin(3);
t96 = pkin(2) * t134;
t133 = t64 * qJ(4) - t96;
t27 = t105 * t66 - t106 * t67;
t43 = t105 * t87 - t155;
t39 = t106 * t81 + t156;
t129 = pkin(1) * t168;
t128 = t100 ^ 2;
t127 = qJD(4) - t39;
t11 = t162 * qJD(2) + t127 - t159;
t17 = t37 + t160;
t123 = t107 * t17 - t109 * t11;
t122 = -t107 * t11 - t109 * t17;
t82 = t105 * t108 - t147;
t41 = t107 * t83 - t109 * t82;
t42 = t107 * t82 + t109 * t83;
t119 = t83 * qJ(4) + t136;
t25 = t70 * pkin(3) + t173;
t118 = t25 * t73 + t20;
t117 = -pkin(2) * t139 - t70 * qJ(4);
t116 = t73 * qJD(4) + t133;
t75 = qJD(2) * t147 - t105 * t140;
t115 = t75 * qJ(4) + t83 * qJD(4) - t137;
t114 = t20 * t83 + t27 * t73 - t28 * t70 + t45 * t64 - t46 * t63;
t97 = t105 * pkin(2) + qJ(4);
t95 = -pkin(4) + t99;
t38 = t82 * pkin(3) - t119;
t31 = -qJD(2) * pkin(3) + t127;
t30 = t82 * pkin(7) + t46;
t29 = -t83 * pkin(7) + t45;
t26 = t73 * pkin(3) - t117;
t22 = t43 + t160;
t19 = t162 * t82 + t119;
t16 = t72 * pkin(3) - t115;
t14 = t72 * pkin(7) + t28;
t13 = -t75 * pkin(7) + t27;
t12 = t162 * t73 + t117;
t9 = -t116 + t161;
t6 = t162 * t72 + t115;
t5 = t42 * qJD(5) + t107 * t75 - t109 * t72;
t4 = -t41 * qJD(5) + t107 * t72 + t109 * t75;
t1 = t162 * t63 + t116;
t2 = [0, 0, 0, 0.2e1 * t110 * t134, t141 * t168, t144, -t145, 0, -pkin(6) * t144 + t108 * t129, pkin(6) * t145 + t110 * t129, -t21 * t82 - t39 * t75 - t40 * t72 + t114, t158 + t21 * t46 - t39 * t27 + t40 * t28 + (t89 - t124) * t137, -t27 * qJD(2) + t16 * t70 + t25 * t72 + t38 * t63 + t9 * t82, -t18 * t82 + t31 * t75 - t37 * t72 + t114, t28 * qJD(2) - t16 * t73 - t25 * t75 - t38 * t64 - t9 * t83, t25 * t16 + t18 * t46 + t31 * t27 + t37 * t28 + t9 * t38 + t158, -t130 * t42 + t165 * t4, -t121 * t4 + t130 * t41 - t165 * t5 - t42 * t3, -t4 * t100, t5 * t100, 0, t6 * t121 + t19 * t3 + t1 * t41 + t10 * t5 - (-t107 * t14 + t109 * t13 + (-t107 * t29 - t109 * t30) * qJD(5)) * t100, t6 * t165 - t19 * t130 + t1 * t42 + t10 * t4 + (t107 * t13 + t109 * t14 + (-t107 * t30 + t109 * t29) * qJD(5)) * t100; 0, 0, 0, -t108 * t146, t141 * t112, 0, 0, 0, t112 * pkin(1) * t108, pkin(1) * t146, (t40 - t43) * t73 + (-t39 + t44) * t70 + (-t105 * t63 - t106 * t64) * pkin(2), t39 * t43 - t40 * t44 + (t105 * t21 - t106 * t20 - t89 * t139) * pkin(2), t43 * qJD(2) - t26 * t70 - t118, -t97 * t63 + t99 * t64 + (t37 - t43) * t73 + (t31 - t142) * t70, -t44 * qJD(2) - t25 * t70 + t26 * t73 + 0.2e1 * t101 + t21, t142 * t37 + t18 * t97 + t20 * t99 - t25 * t26 - t31 * t43, -t166, -t171, t172, t174, 0, -t12 * t121 + (t143 * t107 + t109 * t22) * t100 + (-(-t107 * t95 - t109 * t97) * t100 - t122) * qJD(5) + t170, -t12 * t165 + (-t107 * t22 + t143 * t109) * t100 + ((-t107 * t97 + t109 * t95) * t100 - t123) * qJD(5) - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t39 * t73 + t40 * t70 + t96, 0.2e1 * t73 * qJD(2), t167, t92 + (t70 - t135) * qJD(2), t161 + t37 * t70 + (-qJD(4) - t31) * t73 - t133, 0, 0, 0, 0, 0, -t3 + t150, t130 - t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t70, -t92 + (t70 + t135) * qJD(2), -t69 - t111, -t37 * qJD(2) + t118, 0, 0, 0, 0, 0, -t107 * t128 - t121 * t73, -t109 * t128 - t165 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t171, -t172, -t174, 0, t164 * t122 - t170, t164 * t123 + t163;];
tauc_reg = t2;
