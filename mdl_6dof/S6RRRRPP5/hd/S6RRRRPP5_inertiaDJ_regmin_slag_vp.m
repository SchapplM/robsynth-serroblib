% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:07
% EndTime: 2019-03-09 21:09:15
% DurationCPUTime: 2.55s
% Computational Cost: add. (3648->310), mult. (8818->499), div. (0->0), fcn. (7426->6), ass. (0->141)
t117 = sin(qJ(2));
t116 = sin(qJ(3));
t119 = cos(qJ(2));
t168 = t119 * qJD(2);
t155 = t116 * t168;
t118 = cos(qJ(3));
t171 = qJD(3) * t118;
t193 = t117 * t171 + t155;
t192 = -0.4e1 * t117;
t115 = sin(qJ(4));
t187 = cos(qJ(4));
t151 = qJD(2) * t187;
t142 = t119 * t151;
t172 = qJD(3) * t116;
t160 = t117 * t172;
t177 = t116 * t117;
t163 = t115 * t177;
t149 = t187 * qJD(4);
t190 = t187 * qJD(3) + t149;
t27 = t116 * t142 - t115 * t160 - qJD(4) * t163 + (t115 * t168 + t190 * t117) * t118;
t178 = t115 * t118;
t67 = t187 * t116 + t178;
t53 = t67 * t117;
t191 = t27 * qJ(6) + t53 * qJD(6);
t112 = t117 ^ 2;
t146 = (-t119 ^ 2 + t112) * qJD(2);
t113 = t118 ^ 2;
t174 = t116 ^ 2 - t113;
t147 = t174 * qJD(3);
t189 = qJD(3) + qJD(4);
t107 = t117 * qJD(2);
t170 = qJD(3) * t119;
t159 = t116 * t170;
t129 = t107 * t118 + t159;
t138 = pkin(2) * t117 - pkin(8) * t119;
t70 = t138 * qJD(2);
t139 = -t119 * pkin(2) - t117 * pkin(8);
t74 = -pkin(1) + t139;
t32 = pkin(7) * t129 - t116 * t70 - t74 * t171;
t121 = 2 * qJD(5);
t120 = -pkin(4) - pkin(5);
t188 = -pkin(9) - pkin(8);
t186 = pkin(3) * t117;
t185 = pkin(7) * t116;
t176 = t117 * t118;
t65 = t118 * t74;
t42 = -pkin(9) * t176 + t65 + (-pkin(3) - t185) * t119;
t175 = t118 * t119;
t92 = pkin(7) * t175;
t180 = t116 * t74 + t92;
t47 = -pkin(9) * t177 + t180;
t184 = t115 * t42 + t187 * t47;
t156 = t116 * t107;
t182 = pkin(7) * t156 + t118 * t70;
t179 = t115 * t116;
t79 = t188 * t118;
t49 = t188 * t179 - t187 * t79;
t103 = pkin(3) * t149;
t87 = t103 + qJD(5);
t96 = t115 * pkin(3) + qJ(5);
t181 = t87 * qJ(5) + t96 * qJD(5);
t71 = pkin(3) * t177 + t117 * pkin(7);
t169 = qJD(4) * t115;
t167 = t119 * qJD(5);
t166 = -0.2e1 * pkin(1) * qJD(2);
t165 = -0.2e1 * pkin(2) * qJD(3);
t101 = pkin(7) * t168;
t50 = t193 * pkin(3) + t101;
t164 = pkin(3) * t172;
t102 = pkin(3) * t169;
t100 = -t118 * pkin(3) - pkin(2);
t162 = t188 * qJD(3);
t161 = t187 * t118;
t157 = t118 * t170;
t154 = t116 * t171;
t153 = t117 * t168;
t152 = t118 * t168;
t20 = (-pkin(9) * t175 + t186) * qJD(2) + (-t92 + (pkin(9) * t117 - t74) * t116) * qJD(3) + t182;
t23 = -pkin(9) * t193 - t32;
t148 = t115 * t23 + t47 * t149 + t42 * t169 - t187 * t20;
t145 = 0.2e1 * t153;
t144 = t188 * t187;
t143 = t116 * t152;
t99 = -t187 * pkin(3) - pkin(4);
t54 = t117 * t161 - t163;
t141 = t54 * qJ(5) - t71;
t17 = -t119 * qJ(5) + t184;
t140 = -t115 * t47 + t187 * t42;
t105 = pkin(4) * t107;
t5 = -t105 + t148;
t137 = -qJ(5) * t27 - qJD(5) * t53;
t46 = t189 * t67;
t66 = -t161 + t179;
t136 = -qJ(5) * t46 - qJD(5) * t66;
t135 = t116 * t144;
t134 = qJD(3) * t144;
t18 = t119 * pkin(4) - t140;
t133 = t67 * qJ(5) - t100;
t28 = -qJD(4) * t135 - t116 * t134 - t162 * t178 - t169 * t79;
t132 = t107 * t49 + t28 * t119;
t29 = -t79 * t149 - t118 * t134 + (qJD(4) * t188 + t162) * t179;
t48 = -t115 * t79 - t135;
t131 = -t107 * t48 + t29 * t119;
t6 = -t115 * t20 - t42 * t149 + t169 * t47 - t187 * t23;
t97 = qJ(5) * t107;
t130 = -t6 + t97;
t26 = t115 * t155 + t117 * t46 - t118 * t142;
t1 = -pkin(5) * t107 + t26 * qJ(6) - t54 * qJD(6) + t5;
t127 = -t26 * qJ(5) + t54 * qJD(5) - t50;
t126 = -0.2e1 * t167 - t6 + 0.2e1 * t97;
t45 = -t190 * t118 + t189 * t179;
t125 = -t45 * qJ(5) + t67 * qJD(5) - t164;
t124 = -t102 * t54 + t96 * t27 + t87 * t53;
t123 = -t102 * t67 + t96 * t46 + t87 * t66;
t4 = t130 - t167;
t122 = t96 * t107 + (-qJD(5) - t87) * t119 + t130;
t110 = qJ(5) * t121;
t95 = -0.2e1 * t102;
t93 = -pkin(5) + t99;
t88 = t121 + t103;
t86 = -0.2e1 * t153;
t83 = t119 * t102;
t73 = 0.2e1 * t87;
t62 = t96 * t87;
t39 = t66 * pkin(4) - t133;
t35 = t66 * qJ(6) + t49;
t34 = -t67 * qJ(6) + t48;
t33 = -t180 * qJD(3) + t182;
t31 = t120 * t66 + t133;
t30 = t53 * pkin(4) - t141;
t19 = t120 * t53 + t141;
t14 = t46 * pkin(4) - t125;
t13 = t53 * qJ(6) + t17;
t12 = t119 * pkin(5) - t54 * qJ(6) + t18;
t11 = t45 * qJ(6) - t67 * qJD(6) + t29;
t10 = t46 * qJ(6) + t66 * qJD(6) - t28;
t9 = t120 * t46 + t125;
t8 = t27 * pkin(4) - t127;
t3 = t120 * t27 + t127;
t2 = t4 + t191;
t7 = [0, 0, 0, t145, -0.2e1 * t146, 0, 0, 0, t117 * t166, t119 * t166, -0.2e1 * t112 * t154 + 0.2e1 * t113 * t153, 0.2e1 * t112 * t147 + t143 * t192, 0.2e1 * t117 * t159 + 0.2e1 * t118 * t146, -0.2e1 * t116 * t146 + 0.2e1 * t117 * t157, t86, 0.2e1 * t65 * t107 - 0.2e1 * t33 * t119 + 0.2e1 * (t112 * t171 + t116 * t153) * pkin(7), -0.2e1 * t32 * t119 - 0.2e1 * t180 * t107 + 0.2e1 * (-t112 * t172 + t118 * t145) * pkin(7), -0.2e1 * t54 * t26, 0.2e1 * t26 * t53 - 0.2e1 * t54 * t27, 0.2e1 * t107 * t54 + 0.2e1 * t26 * t119, -0.2e1 * t107 * t53 + 0.2e1 * t27 * t119, t86, 0.2e1 * t107 * t140 + 0.2e1 * t119 * t148 + 0.2e1 * t71 * t27 + 0.2e1 * t50 * t53, -0.2e1 * t184 * t107 - 0.2e1 * t6 * t119 - 0.2e1 * t71 * t26 + 0.2e1 * t50 * t54, -0.2e1 * t107 * t18 + 0.2e1 * t5 * t119 + 0.2e1 * t30 * t27 + 0.2e1 * t8 * t53, -0.2e1 * t17 * t27 - 0.2e1 * t18 * t26 - 0.2e1 * t4 * t53 + 0.2e1 * t5 * t54, 0.2e1 * t107 * t17 - 0.2e1 * t4 * t119 + 0.2e1 * t30 * t26 - 0.2e1 * t8 * t54, 0.2e1 * t17 * t4 + 0.2e1 * t18 * t5 + 0.2e1 * t30 * t8, 0.2e1 * t1 * t119 - 0.2e1 * t107 * t12 - 0.2e1 * t19 * t27 - 0.2e1 * t3 * t53, 0.2e1 * t107 * t13 - 0.2e1 * t2 * t119 - 0.2e1 * t19 * t26 + 0.2e1 * t3 * t54, -0.2e1 * t1 * t54 + 0.2e1 * t12 * t26 + 0.2e1 * t13 * t27 + 0.2e1 * t2 * t53, 0.2e1 * t12 * t1 + 0.2e1 * t13 * t2 + 0.2e1 * t19 * t3; 0, 0, 0, 0, 0, t168, -t107, 0, -t101, pkin(7) * t107, -t117 * t147 + t143, t154 * t192 - t174 * t168, t156 - t157, t129, 0 (pkin(8) * t175 + (-pkin(2) * t118 + t185) * t117) * qJD(3) + (t116 * t139 - t92) * qJD(2) (pkin(7) * t176 + t116 * t138) * qJD(3) + (t118 * t139 + t119 * t185) * qJD(2), -t26 * t67 - t54 * t45, t26 * t66 - t67 * t27 + t45 * t53 - t54 * t46, t107 * t67 + t45 * t119, -t107 * t66 + t46 * t119, 0, t100 * t27 + t164 * t53 + t71 * t46 + t50 * t66 + t131, -t100 * t26 + t164 * t54 - t71 * t45 + t50 * t67 - t132, t14 * t53 + t39 * t27 + t30 * t46 + t8 * t66 + t131, -t17 * t46 - t18 * t45 - t48 * t26 - t49 * t27 + t28 * t53 + t29 * t54 - t4 * t66 + t5 * t67, -t14 * t54 + t39 * t26 + t30 * t45 - t8 * t67 + t132, t30 * t14 - t17 * t28 + t18 * t29 + t8 * t39 + t4 * t49 + t5 * t48, -t107 * t34 + t11 * t119 - t19 * t46 - t31 * t27 - t3 * t66 - t9 * t53, -t10 * t119 + t107 * t35 - t19 * t45 - t31 * t26 + t3 * t67 + t9 * t54, -t1 * t67 + t10 * t53 - t11 * t54 + t12 * t45 + t13 * t46 + t2 * t66 + t34 * t26 + t35 * t27, t1 * t34 + t13 * t10 + t12 * t11 + t19 * t9 + t2 * t35 + t3 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t154, -0.2e1 * t147, 0, 0, 0, t116 * t165, t118 * t165, -0.2e1 * t67 * t45, 0.2e1 * t45 * t66 - 0.2e1 * t67 * t46, 0, 0, 0, 0.2e1 * t100 * t46 + 0.2e1 * t164 * t66, -0.2e1 * t100 * t45 + 0.2e1 * t164 * t67, 0.2e1 * t14 * t66 + 0.2e1 * t39 * t46, 0.2e1 * t28 * t66 + 0.2e1 * t29 * t67 - 0.2e1 * t48 * t45 - 0.2e1 * t49 * t46, -0.2e1 * t14 * t67 + 0.2e1 * t39 * t45, 0.2e1 * t39 * t14 - 0.2e1 * t49 * t28 + 0.2e1 * t48 * t29, -0.2e1 * t31 * t46 - 0.2e1 * t9 * t66, -0.2e1 * t31 * t45 + 0.2e1 * t9 * t67, 0.2e1 * t10 * t66 - 0.2e1 * t11 * t67 + 0.2e1 * t34 * t45 + 0.2e1 * t35 * t46, 0.2e1 * t35 * t10 + 0.2e1 * t34 * t11 + 0.2e1 * t31 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152 - t160, -t193, t107, t33, t32, 0, 0, -t26, -t27, t107, t151 * t186 - t148 + t83 (-t107 * t115 + t119 * t149) * pkin(3) + t6, -t107 * t99 - t5 + t83, -t99 * t26 - t124, t122, t102 * t18 + t17 * t87 + t4 * t96 + t5 * t99, -t107 * t93 - t1 + t83, t122 + t191, t93 * t26 + t124, t1 * t93 + t102 * t12 + t13 * t87 + t2 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t172, 0, -pkin(8) * t171, pkin(8) * t172, 0, 0, -t45, -t46, 0, -t29, t28, -t29, -t99 * t45 - t123, -t28, t102 * t48 - t28 * t96 + t29 * t99 + t49 * t87, -t11, t10, t93 * t45 + t123, t10 * t96 + t102 * t34 + t11 * t93 + t35 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -0.2e1 * t103, t95, 0, t73, 0.2e1 * t102 * t99 + 0.2e1 * t62, t95, t73, 0, 0.2e1 * t102 * t93 + 0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, t107, -t148, t6, 0.2e1 * t105 - t148, pkin(4) * t26 + t137, t126, -t5 * pkin(4) + t4 * qJ(5) + t17 * qJD(5), -t107 * t120 - t1, t126 + t191, t120 * t26 - t137, t2 * qJ(5) + t13 * qJD(5) + t1 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, -t29, t28, -t29, pkin(4) * t45 + t136, -t28, -t29 * pkin(4) - t28 * qJ(5) + t49 * qJD(5), -t11, t10, t120 * t45 - t136, t10 * qJ(5) + t35 * qJD(5) + t11 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t103, -t102, 0, t88, -pkin(4) * t102 + t181, -t102, t88, 0, t102 * t120 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t110, 0, t121, 0, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t26, 0, t5, -t107, 0, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, t29, 0, 0, t45, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
