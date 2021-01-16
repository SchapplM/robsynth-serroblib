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
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 19:36:30
% EndTime: 2021-01-15 19:36:36
% DurationCPUTime: 1.15s
% Computational Cost: add. (1265->206), mult. (3374->275), div. (0->0), fcn. (2341->6), ass. (0->131)
t101 = qJD(2) - qJD(5);
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t83 = t106 * t111 + t107 * t109;
t154 = qJD(1) * t83;
t153 = t107 * t111;
t140 = qJD(1) * t153;
t146 = qJD(1) * t109;
t70 = t106 * t146 - t140;
t174 = t108 * t70 + t110 * t154;
t157 = t174 * t101;
t72 = t83 * qJD(2);
t64 = qJD(1) * t72;
t144 = qJD(1) * qJD(2);
t138 = t111 * t144;
t139 = t109 * t144;
t92 = t106 * t139;
t65 = t107 * t138 - t92;
t3 = qJD(5) * t174 + t108 * t65 - t110 * t64;
t184 = t3 + t157;
t141 = t111 * pkin(2) + pkin(1);
t181 = t141 * qJD(1);
t89 = qJD(3) - t181;
t183 = -t154 * qJ(4) + t89;
t124 = t108 * t154 - t110 * t70;
t134 = t124 * qJD(5) - t108 * t64 - t110 * t65;
t156 = t124 * t101;
t182 = t134 + t156;
t180 = -t124 ^ 2 + t174 ^ 2;
t171 = -pkin(3) - pkin(4);
t10 = t171 * t70 - t183;
t102 = qJD(2) * qJD(4);
t164 = -qJ(3) - pkin(6);
t135 = qJD(2) * t164;
t67 = t111 * qJD(3) + t109 * t135;
t54 = t67 * qJD(1);
t118 = -t109 * qJD(3) + t111 * t135;
t55 = t118 * qJD(1);
t21 = t106 * t55 + t107 * t54;
t18 = t102 + t21;
t7 = t64 * pkin(7) + t18;
t20 = t106 * t54 - t107 * t55;
t8 = -t65 * pkin(7) + t20;
t179 = t10 * t174 + t108 * t7 - t110 * t8;
t177 = -0.2e1 * t144;
t69 = t154 ^ 2;
t176 = -t70 ^ 2 - t69;
t175 = t174 * t124;
t91 = t164 * t111;
t88 = qJD(1) * t91;
t163 = t106 * t88;
t90 = t164 * t109;
t87 = qJD(1) * t90;
t45 = t107 * t87 + t163;
t148 = qJD(4) - t45;
t173 = qJD(5) + t101;
t172 = t10 * t124 - t108 * t8 - t110 * t7;
t170 = t64 * pkin(3);
t169 = t70 * pkin(7);
t168 = t154 * pkin(7);
t167 = pkin(2) * t109;
t46 = -t106 * t91 - t107 * t90;
t166 = t20 * t46;
t25 = t70 * pkin(3) + t183;
t165 = t25 * t154;
t28 = t106 * t118 + t107 * t67;
t162 = t107 * t88;
t81 = qJD(2) * pkin(2) + t87;
t41 = t106 * t81 - t162;
t47 = t106 * t90 - t107 * t91;
t113 = qJD(1) ^ 2;
t152 = t111 * t113;
t112 = qJD(2) ^ 2;
t151 = t112 * t109;
t150 = t112 * t111;
t149 = -t168 + t148;
t147 = t109 ^ 2 - t111 ^ 2;
t145 = qJD(2) * t109;
t37 = qJD(2) * qJ(4) + t41;
t143 = pkin(2) * t145;
t142 = pkin(2) * t146;
t99 = -t107 * pkin(2) - pkin(3);
t96 = pkin(2) * t139;
t137 = t65 * qJ(4) - t96;
t27 = t106 * t67 - t107 * t118;
t44 = t106 * t87 - t162;
t40 = t107 * t81 + t163;
t133 = 0.2e1 * t154;
t132 = pkin(1) * t177;
t131 = t101 ^ 2;
t130 = t45 * qJD(2) - t21;
t129 = qJD(4) - t40;
t11 = t171 * qJD(2) + t129 - t168;
t17 = t37 + t169;
t126 = t108 * t17 - t110 * t11;
t125 = -t108 * t11 - t110 * t17;
t82 = t106 * t109 - t153;
t42 = t108 * t83 - t110 * t82;
t43 = t108 * t82 + t110 * t83;
t122 = t83 * qJ(4) + t141;
t121 = -t70 * qJ(4) - t142;
t120 = qJD(4) * t154 + t137;
t119 = t44 * qJD(2) - t20;
t75 = qJD(2) * t153 - t106 * t145;
t117 = t75 * qJ(4) + t83 * qJD(4) - t143;
t116 = t154 * t27 + t20 * t83 - t28 * t70 + t46 * t65 - t47 * t64;
t115 = t133 * qJD(2);
t97 = t106 * pkin(2) + qJ(4);
t95 = -pkin(4) + t99;
t39 = t82 * pkin(3) - t122;
t38 = -t92 + (-t70 + t140) * qJD(2);
t31 = -qJD(2) * pkin(3) + t129;
t30 = t82 * pkin(7) + t47;
t29 = -t83 * pkin(7) + t46;
t26 = pkin(3) * t154 - t121;
t22 = t44 + t169;
t19 = t171 * t82 + t122;
t16 = t72 * pkin(3) - t117;
t14 = t72 * pkin(7) + t28;
t13 = -t75 * pkin(7) + t27;
t12 = t154 * t171 + t121;
t9 = -t120 + t170;
t6 = t171 * t72 + t117;
t5 = t43 * qJD(5) + t108 * t75 - t110 * t72;
t4 = -t42 * qJD(5) + t108 * t72 + t110 * t75;
t1 = t171 * t64 + t120;
t2 = [0, 0, 0, 0.2e1 * t109 * t138, t147 * t177, t150, -t151, 0, -pkin(6) * t150 + t109 * t132, pkin(6) * t151 + t111 * t132, -t141 * t64 + t89 * t72 + (-t27 + (qJD(1) * t82 + t70) * t167) * qJD(2), -t141 * t65 + t89 * t75 + (t133 * t167 - t28) * qJD(2), -t21 * t82 - t40 * t75 - t41 * t72 + t116, t166 + t21 * t47 - t40 * t27 + t41 * t28 + (t89 - t181) * t143, -t27 * qJD(2) + t16 * t70 + t25 * t72 + t39 * t64 + t9 * t82, -t18 * t82 + t31 * t75 - t37 * t72 + t116, t28 * qJD(2) - t154 * t16 - t25 * t75 - t39 * t65 - t9 * t83, t25 * t16 + t18 * t47 + t31 * t27 + t37 * t28 + t9 * t39 + t166, -t134 * t43 + t174 * t4, -t124 * t4 + t134 * t42 - t174 * t5 - t43 * t3, -t4 * t101, t5 * t101, 0, t6 * t124 + t19 * t3 + t1 * t42 + t10 * t5 - (-t108 * t14 + t110 * t13 + (-t108 * t29 - t110 * t30) * qJD(5)) * t101, t6 * t174 - t19 * t134 + t1 * t43 + t10 * t4 + (t108 * t13 + t110 * t14 + (-t108 * t30 + t110 * t29) * qJD(5)) * t101; 0, 0, 0, -t109 * t152, t147 * t113, 0, 0, 0, t113 * pkin(1) * t109, pkin(1) * t152, -t70 * t142 - t154 * t89 + t119, -t142 * t154 + t89 * t70 + t130, (t41 - t44) * t154 + (-t40 + t45) * t70 + (-t106 * t64 - t107 * t65) * pkin(2), t40 * t44 - t41 * t45 + (t106 * t21 - t107 * t20 - t89 * t146) * pkin(2), -t26 * t70 + t119 - t165, -t97 * t64 + t99 * t65 + (t37 - t44) * t154 + (t31 - t148) * t70, t154 * t26 - t25 * t70 + 0.2e1 * t102 - t130, t148 * t37 + t18 * t97 + t20 * t99 - t25 * t26 - t31 * t44, -t175, -t180, t182, t184, 0, -t12 * t124 + (t108 * t149 + t110 * t22) * t101 + (-(-t108 * t95 - t110 * t97) * t101 - t125) * qJD(5) + t179, -t12 * t174 + (-t108 * t22 + t110 * t149) * t101 + ((-t108 * t97 + t110 * t95) * t101 - t126) * qJD(5) - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t38, t176, t154 * t40 + t41 * t70 + t96, t115, t176, -t38, t170 + t37 * t70 + (-qJD(4) - t31) * t154 - t137, 0, 0, 0, 0, 0, -t3 + t157, t134 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154 * t70, -t92 + (t70 + t140) * qJD(2), -t69 - t112, -t37 * qJD(2) + t165 + t20, 0, 0, 0, 0, 0, -t108 * t131 - t124 * t154, -t110 * t131 - t154 * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t180, -t182, -t184, 0, t173 * t125 - t179, t173 * t126 + t172;];
tauc_reg = t2;
