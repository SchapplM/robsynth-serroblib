% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:31
% EndTime: 2019-03-09 03:59:40
% DurationCPUTime: 3.17s
% Computational Cost: add. (4848->279), mult. (9690->491), div. (0->0), fcn. (9442->8), ass. (0->144)
t153 = sin(pkin(10));
t154 = cos(pkin(10));
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t62 = t153 * t82 + t154 * t80;
t63 = -t153 * t80 + t154 * t82;
t73 = t80 * pkin(3) + qJ(2);
t104 = t62 * pkin(4) - pkin(8) * t63 + t73;
t175 = -pkin(1) - pkin(7);
t146 = qJ(4) - t175;
t148 = t82 * qJD(3);
t100 = -t80 * qJD(4) - t146 * t148;
t149 = t80 * qJD(3);
t99 = -t82 * qJD(4) + t146 * t149;
t89 = t154 * t100 + t153 * t99;
t184 = -qJD(5) * t104 - t89;
t57 = t63 * qJD(3);
t58 = t62 * qJD(3);
t67 = pkin(3) * t148 + qJD(2);
t114 = t146 * t153;
t115 = t146 * t154;
t91 = -t82 * t114 - t80 * t115;
t183 = -t57 * pkin(4) - t58 * pkin(8) + qJD(5) * t91 - t67;
t182 = (-t153 * t57 + t154 * t58) * pkin(3);
t79 = sin(qJ(5));
t76 = t79 ^ 2;
t81 = cos(qJ(5));
t77 = t81 ^ 2;
t155 = -t77 - t76;
t173 = cos(qJ(6));
t134 = t173 * t81;
t78 = sin(qJ(6));
t163 = t78 * t79;
t109 = t134 - t163;
t34 = t109 * t62;
t123 = t153 * pkin(3) + pkin(8);
t181 = t123 * t57;
t130 = t173 * qJD(6);
t180 = t173 * qJD(5) + t130;
t151 = qJD(5) * t81;
t161 = t79 * t58;
t113 = t63 * t151 - t161;
t152 = qJD(5) * t79;
t140 = t63 * t152;
t160 = t81 * t58;
t179 = t140 + t160;
t156 = t76 - t77;
t127 = qJD(5) * t156;
t178 = qJD(5) + qJD(6);
t21 = t81 * t104 - t79 * t91;
t22 = t79 * t104 + t81 * t91;
t120 = t21 * t79 - t22 * t81;
t8 = t183 * t79 + t184 * t81;
t9 = -t183 * t81 + t184 * t79;
t177 = t120 * qJD(5) + t79 * t8 - t81 * t9;
t61 = t63 ^ 2;
t176 = 0.2e1 * qJD(2);
t174 = t57 * pkin(5);
t135 = t173 * t79;
t65 = t78 * t81 + t135;
t44 = t178 * t65;
t172 = t44 * t62;
t171 = t44 * t63;
t28 = t153 * t100 - t154 * t99;
t45 = -t80 * t114 + t82 * t115;
t170 = t45 * t28;
t169 = t57 * t65;
t72 = -t154 * pkin(3) - pkin(4);
t168 = t58 * t72;
t167 = t62 * t57;
t46 = t63 * t58;
t166 = t63 * t79;
t165 = t109 * t44;
t43 = t178 * t163 - t180 * t81;
t164 = t65 * t43;
t162 = t79 * t57;
t159 = t81 * t63;
t144 = t63 * t163;
t14 = -t58 * t135 - t78 * t140 - qJD(6) * t144 + (t180 * t63 - t58 * t78) * t81;
t33 = t65 * t63;
t158 = -t65 * t14 + t43 * t33;
t150 = qJD(6) * t78;
t147 = qJ(2) * qJD(3);
t42 = 0.2e1 * t167;
t145 = t79 * t160;
t143 = 0.2e1 * qJD(5) * t72;
t142 = pkin(5) * t152;
t141 = pkin(5) * t150;
t138 = t79 * t151;
t137 = t80 * t148;
t136 = t62 ^ 2 + t61;
t133 = t175 * qJD(3);
t132 = t155 * t57;
t126 = t61 * t138;
t125 = pkin(5) * t130;
t12 = t58 * t134 - t78 * t161 + t171;
t35 = t63 * t134 - t144;
t122 = t109 * t12 + t35 * t44;
t121 = t21 * t81 + t22 * t79;
t119 = -t28 * t63 + t45 * t58;
t118 = t43 * t62 - t169;
t117 = t46 - t167;
t116 = pkin(9) + t123;
t38 = t62 * t151 + t162;
t16 = -pkin(9) * t166 + t22;
t83 = t179 * pkin(9) + t174 + t9;
t87 = t62 * pkin(5) - pkin(9) * t159 + t21;
t85 = t173 * t87;
t88 = -t113 * pkin(9) - t8;
t1 = -qJD(6) * t85 + t16 * t150 - t173 * t88 - t78 * t83;
t111 = t79 * t123;
t110 = t81 * t123;
t108 = qJD(5) * t123;
t107 = t78 * t116;
t106 = 0.2e1 * t117;
t105 = t79 * t107;
t103 = qJD(5) * t107;
t101 = t116 * t173;
t98 = t79 * t101;
t97 = -t168 - t181;
t96 = t123 * t62 - t72 * t63;
t94 = qJD(5) * t101;
t92 = -t121 * qJD(5) - t9 * t79 - t8 * t81;
t86 = t78 * t87;
t84 = t91 * t57 + t89 * t62 + t119;
t2 = -qJD(6) * t86 - t16 * t130 + t173 * t83 - t78 * t88;
t75 = qJ(2) * t176;
t66 = -t81 * pkin(5) + t72;
t60 = t116 * t81;
t40 = t173 * t60 - t105;
t39 = -t78 * t60 - t98;
t37 = t62 * t152 - t81 * t57;
t32 = t65 * t62;
t26 = pkin(5) * t166 + t45;
t24 = t63 * t127 + t145;
t20 = t109 * t57 - t172;
t19 = qJD(6) * t105 + t79 * t103 - t60 * t130 - t81 * t94;
t18 = qJD(6) * t98 + t81 * t103 + t60 * t150 + t79 * t94;
t17 = t113 * pkin(5) + t28;
t13 = -t178 * t34 - t169;
t11 = -t57 * t134 + t78 * t162 + t172;
t7 = t173 * t16 + t86;
t6 = -t78 * t16 + t85;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t75, -0.2e1 * t137, 0.2e1 * (t80 ^ 2 - t82 ^ 2) * qJD(3), 0, 0.2e1 * t137, 0, 0, 0.2e1 * qJD(2) * t80 + 0.2e1 * t82 * t147, 0.2e1 * qJD(2) * t82 - 0.2e1 * t80 * t147, 0, t75, -0.2e1 * t46, -0.2e1 * t57 * t63 + 0.2e1 * t58 * t62, 0, t42, 0, 0, 0.2e1 * t57 * t73 + 0.2e1 * t62 * t67, -0.2e1 * t58 * t73 + 0.2e1 * t63 * t67, -0.2e1 * t84, 0.2e1 * t73 * t67 + 0.2e1 * t91 * t89 + 0.2e1 * t170, -0.2e1 * t77 * t46 - 0.2e1 * t126, 0.2e1 * t61 * t127 + 0.4e1 * t63 * t145, 0.2e1 * t57 * t159 - 0.2e1 * t179 * t62, -0.2e1 * t76 * t46 + 0.2e1 * t126, -0.2e1 * t113 * t62 - 0.2e1 * t63 * t162, t42, 0.2e1 * t113 * t45 + 0.2e1 * t28 * t166 + 0.2e1 * t21 * t57 + 0.2e1 * t62 * t9, 0.2e1 * t28 * t159 - 0.2e1 * t179 * t45 - 0.2e1 * t22 * t57 + 0.2e1 * t62 * t8, 0.2e1 * t121 * t58 + 0.2e1 * t177 * t63, 0.2e1 * t21 * t9 - 0.2e1 * t22 * t8 + 0.2e1 * t170, -0.2e1 * t35 * t12, 0.2e1 * t12 * t33 - 0.2e1 * t14 * t35, -0.2e1 * t12 * t62 + 0.2e1 * t35 * t57, 0.2e1 * t33 * t14, -0.2e1 * t14 * t62 - 0.2e1 * t33 * t57, t42, 0.2e1 * t14 * t26 + 0.2e1 * t17 * t33 + 0.2e1 * t2 * t62 + 0.2e1 * t57 * t6, 0.2e1 * t1 * t62 - 0.2e1 * t12 * t26 + 0.2e1 * t17 * t35 - 0.2e1 * t57 * t7, 0.2e1 * t1 * t33 + 0.2e1 * t12 * t6 - 0.2e1 * t14 * t7 - 0.2e1 * t2 * t35, -0.2e1 * t1 * t7 + 0.2e1 * t17 * t26 + 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t84, 0, 0, 0, 0, 0, 0, t79 * t106 - t136 * t151, t81 * t106 + t136 * t152, 0, -t120 * t57 + t92 * t62 + t119, 0, 0, 0, 0, 0, 0, t13 * t62 - t14 * t63 - t32 * t57 + t33 * t58, t11 * t62 + t12 * t63 - t34 * t57 + t35 * t58, t11 * t33 - t12 * t32 - t13 * t35 - t14 * t34, -t1 * t34 - t11 * t7 + t13 * t6 - t17 * t63 - t2 * t32 + t26 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62 * t132 - 0.2e1 * t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t11 * t34 - 0.2e1 * t13 * t32 - 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, 0, -t148, 0, -t80 * t133, -t82 * t133, 0, 0, 0, 0, -t58, 0, -t57, 0, -t28, -t89, t182 (t89 * t153 - t28 * t154) * pkin(3), -t24, -0.4e1 * t63 * t138 + t156 * t58, t38, t24, -t37, 0, -t28 * t81 + t97 * t79 + (t45 * t79 - t96 * t81) * qJD(5), t28 * t79 + t97 * t81 + (t45 * t81 + t96 * t79) * qJD(5), t92, -t8 * t110 - t9 * t111 + t28 * t72 + (-t21 * t110 - t22 * t111) * qJD(5), -t12 * t65 - t35 * t43, -t122 + t158, -t118, -t109 * t14 + t33 * t44, t20, 0, -t109 * t17 + t14 * t66 + t142 * t33 + t19 * t62 + t26 * t44 + t39 * t57, -t12 * t66 + t142 * t35 + t17 * t65 + t18 * t62 - t26 * t43 - t40 * t57, -t1 * t109 + t12 * t39 - t14 * t40 + t18 * t33 - t19 * t35 - t2 * t65 + t43 * t6 - t44 * t7, -t1 * t40 + t142 * t26 + t17 * t66 - t18 * t7 + t19 * t6 + t2 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t148, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, -t182, 0, 0, 0, 0, 0, 0, -t179, -t113, -t132, -t155 * t181 + t168, 0, 0, 0, 0, 0, 0, -t109 * t58 - t171, t43 * t63 + t58 * t65, -t109 * t11 - t13 * t65 - t32 * t43 - t34 * t44, -pkin(5) * t140 - t11 * t40 + t13 * t39 - t18 * t34 - t19 * t32 + t58 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t138, -0.2e1 * t127, 0, -0.2e1 * t138, 0, 0, t79 * t143, t81 * t143, 0, 0, -0.2e1 * t164, -0.2e1 * t109 * t43 - 0.2e1 * t65 * t44, 0, -0.2e1 * t165, 0, 0, -0.2e1 * t109 * t142 + 0.2e1 * t44 * t66, 0.2e1 * t142 * t65 - 0.2e1 * t43 * t66, -0.2e1 * t109 * t18 - 0.2e1 * t19 * t65 + 0.2e1 * t39 * t43 - 0.2e1 * t40 * t44, 0.2e1 * t142 * t66 - 0.2e1 * t18 * t40 + 0.2e1 * t19 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t58, 0, t67, 0, 0, 0, 0, 0, 0, -t37, -t38, -t155 * t58, -t177, 0, 0, 0, 0, 0, 0, t20, t118, t122 + t158, -t1 * t65 + t109 * t2 - t43 * t7 - t44 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t13 - t11 * t65 + t32 * t44 - t34 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t19 - t18 * t65 - t39 * t44 - t40 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t164 - 0.2e1 * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, 0, -t113, t57, t9, t8, 0, 0, 0, 0, -t12, 0, -t14, t57, -t62 * t141 + t173 * t174 + t2 (-t130 * t62 - t57 * t78) * pkin(5) + t1 (t173 * t12 - t14 * t78 + (-t173 * t33 + t35 * t78) * qJD(6)) * pkin(5) (t173 * t2 - t1 * t78 + (t173 * t7 - t6 * t78) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37, 0, 0, 0, 0, 0, 0, 0, 0, t13, t11, 0 (t173 * t13 - t11 * t78 + (t173 * t34 + t32 * t78) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, 0, -t152, 0, -t81 * t108, t79 * t108, 0, 0, 0, 0, -t43, 0, -t44, 0, t19, t18 (t173 * t43 - t44 * t78 + (t109 * t173 + t65 * t78) * qJD(6)) * pkin(5) (t173 * t19 - t18 * t78 + (t173 * t40 - t39 * t78) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0 (-t173 * t44 - t43 * t78 + (-t109 * t78 + t173 * t65) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t141, -0.2e1 * t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t14, t57, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, -t44, 0, t19, t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
