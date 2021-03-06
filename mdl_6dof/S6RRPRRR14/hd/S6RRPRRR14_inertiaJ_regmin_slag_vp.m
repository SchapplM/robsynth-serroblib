% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:42:46
% EndTime: 2019-05-09 11:43:06
% DurationCPUTime: 2.98s
% Computational Cost: add. (6714->248), mult. (19015->519), div. (0->0), fcn. (22622->16), ass. (0->145)
t100 = sin(pkin(8));
t104 = cos(pkin(8));
t101 = sin(pkin(7));
t105 = cos(pkin(7));
t106 = cos(pkin(6));
t102 = sin(pkin(6));
t113 = cos(qJ(2));
t135 = t102 * t113;
t120 = -t101 * t135 + t106 * t105;
t103 = cos(pkin(14));
t131 = t105 * t113;
t123 = t102 * t131;
t137 = t101 * t106;
t118 = t123 + t137;
t110 = sin(qJ(2));
t136 = t102 * t110;
t99 = sin(pkin(14));
t86 = t99 * t136;
t174 = -t103 * t118 + t86;
t175 = -t120 * t100 + t104 * t174;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t109 = sin(qJ(4));
t164 = cos(qJ(4));
t57 = t99 * t137 + (t103 * t110 + t99 * t131) * t102;
t39 = -t109 * t175 + t57 * t164;
t45 = -t100 * t174 - t120 * t104;
t25 = t39 * t108 + t45 * t112;
t173 = -0.2e1 * t25;
t38 = t57 * t109 + t175 * t164;
t172 = 0.2e1 * t38;
t171 = -0.2e1 * t39;
t132 = t105 * t100;
t134 = t103 * t104;
t59 = t109 * t132 + (t109 * t134 + t164 * t99) * t101;
t138 = t101 * t103;
t75 = t100 * t138 - t104 * t105;
t47 = t108 * t59 + t112 * t75;
t170 = -0.2e1 * t47;
t121 = t104 * t164;
t122 = t100 * t164;
t155 = t101 * t99;
t58 = -t105 * t122 + t109 * t155 - t121 * t138;
t169 = 0.2e1 * t58;
t168 = -0.2e1 * t59;
t167 = -0.2e1 * t108;
t166 = 0.2e1 * t112;
t94 = t101 ^ 2;
t165 = pkin(2) * t94;
t163 = pkin(1) * t110;
t111 = cos(qJ(6));
t162 = pkin(5) * t111;
t161 = pkin(11) * t100;
t160 = pkin(11) * t104;
t107 = sin(qJ(6));
t159 = pkin(12) * t107;
t158 = t101 * pkin(2);
t154 = t105 * t99;
t82 = pkin(10) * t135 + t106 * t163;
t56 = qJ(3) * t118 + t82;
t91 = t106 * t113 * pkin(1);
t61 = t106 * pkin(2) + t91 + (-qJ(3) * t105 - pkin(10)) * t136;
t140 = qJ(3) * t101;
t69 = (-pkin(2) * t113 - t110 * t140 - pkin(1)) * t102;
t36 = t103 * t56 + t61 * t154 + t69 * t155;
t24 = -pkin(11) * t175 + t36;
t133 = t103 * t105;
t35 = t61 * t133 + t69 * t138 - t99 * t56;
t27 = t120 * pkin(3) - t57 * t160 + t35;
t43 = -t101 * t61 + t105 * t69;
t29 = pkin(3) * t174 - t57 * t161 + t43;
t14 = t164 * t24 + (t100 * t29 + t104 * t27) * t109;
t11 = -t45 * pkin(12) + t14;
t15 = -t100 * t27 + t104 * t29;
t12 = t38 * pkin(4) - t39 * pkin(12) + t15;
t5 = -t108 * t11 + t112 * t12;
t3 = -t38 * pkin(5) - t5;
t157 = t3 * t107;
t156 = t3 * t111;
t78 = pkin(2) * t154 + qJ(3) * t138;
t153 = t107 * t25;
t152 = t107 * t47;
t151 = t108 * t38;
t150 = t108 * t58;
t149 = t111 * t25;
t148 = t111 * t47;
t147 = t112 * t38;
t146 = t112 * t58;
t95 = t102 ^ 2;
t145 = t113 * t95;
t90 = pkin(2) * t133;
t60 = t105 * pkin(3) + t90 + (-qJ(3) - t160) * t155;
t68 = (-pkin(3) * t103 - t99 * t161 - pkin(2)) * t101;
t42 = -t100 * t60 + t104 * t68;
t30 = t58 * pkin(4) - t59 * pkin(12) + t42;
t55 = (t101 * t134 + t132) * pkin(11) + t78;
t34 = t164 * t55 + (t100 * t68 + t104 * t60) * t109;
t32 = -t75 * pkin(12) + t34;
t20 = -t108 * t32 + t112 * t30;
t16 = -t58 * pkin(5) - t20;
t144 = t16 * t107;
t143 = t16 * t111;
t26 = -t45 * t108 + t39 * t112;
t19 = t38 * t107 + t26 * t111;
t142 = t19 * t107;
t48 = -t108 * t75 + t112 * t59;
t41 = t107 * t58 + t111 * t48;
t141 = t41 * t107;
t139 = t100 * t109;
t130 = t107 * t108;
t129 = t107 * t111;
t128 = t107 * t112;
t127 = t111 * t108;
t126 = t111 * t112;
t125 = 0.2e1 * t102 * t106;
t124 = t108 * t166;
t6 = t108 * t12 + t112 * t11;
t21 = t108 * t30 + t112 * t32;
t13 = -t109 * t24 + t27 * t121 + t29 * t122;
t33 = -t109 * t55 + t60 * t121 + t68 * t122;
t10 = t45 * pkin(4) - t13;
t31 = t75 * pkin(4) - t33;
t98 = t111 ^ 2;
t97 = t108 ^ 2;
t96 = t107 ^ 2;
t84 = -t112 * pkin(5) - t108 * pkin(13) - pkin(4);
t81 = -pkin(10) * t136 + t91;
t80 = t108 * t104 + t112 * t139;
t79 = -t112 * t104 + t108 * t139;
t77 = -t99 * t140 + t90;
t72 = pkin(12) * t126 + t107 * t84;
t71 = -pkin(12) * t128 + t111 * t84;
t65 = -t107 * t122 + t111 * t80;
t64 = -t107 * t80 - t111 * t122;
t40 = t107 * t48 - t111 * t58;
t22 = t47 * pkin(5) - t48 * pkin(13) + t31;
t18 = t26 * t107 - t38 * t111;
t17 = t58 * pkin(13) + t21;
t9 = t107 * t22 + t111 * t17;
t8 = -t107 * t17 + t111 * t22;
t7 = t25 * pkin(5) - t26 * pkin(13) + t10;
t4 = t38 * pkin(13) + t6;
t2 = t107 * t7 + t111 * t4;
t1 = -t107 * t4 + t111 * t7;
t23 = [1, 0, 0, t95 * t110 ^ 2, 0.2e1 * t110 * t145, t110 * t125, t113 * t125, t106 ^ 2, 0.2e1 * pkin(1) * t145 + 0.2e1 * t81 * t106, -0.2e1 * t82 * t106 - 0.2e1 * t95 * t163, 0.2e1 * t120 * t35 + 0.2e1 * t174 * t43, -0.2e1 * t120 * t36 + 0.2e1 * t43 * t57, -0.2e1 * t174 * t36 - 0.2e1 * t35 * t57, t35 ^ 2 + t36 ^ 2 + t43 ^ 2, t39 ^ 2, t38 * t171, t45 * t171, t45 * t172, t45 ^ 2, -0.2e1 * t13 * t45 + 0.2e1 * t15 * t38, 0.2e1 * t14 * t45 + 0.2e1 * t15 * t39, t26 ^ 2, t26 * t173, t26 * t172, t38 * t173, t38 ^ 2, 0.2e1 * t10 * t25 + 0.2e1 * t5 * t38, 0.2e1 * t10 * t26 - 0.2e1 * t6 * t38, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t25, t18 * t173, t25 ^ 2, 0.2e1 * t1 * t25 + 0.2e1 * t3 * t18, 0.2e1 * t3 * t19 - 0.2e1 * t2 * t25; 0, 0, 0, 0, 0, t136, t135, t106, t81, -t82, t35 * t105 + t120 * t77 - t43 * t138 - t158 * t174, -t78 * t120 - t36 * t105 + (-pkin(2) * t57 + t43 * t99) * t101, -t77 * t57 + t78 * (t103 * t123 - t86) + (-t35 * t99 + (t106 * t78 + t36) * t103) * t101, -t43 * t158 + t35 * t77 + t36 * t78, t39 * t59, -t38 * t59 - t58 * t39, -t75 * t39 - t45 * t59, t75 * t38 + t45 * t58, t45 * t75, -t13 * t75 + t15 * t58 - t33 * t45 + t42 * t38, t14 * t75 + t15 * t59 + t34 * t45 + t42 * t39, t26 * t48, -t48 * t25 - t26 * t47, t26 * t58 + t48 * t38, -t25 * t58 - t47 * t38, t38 * t58, t10 * t47 + t20 * t38 + t31 * t25 + t5 * t58, t10 * t48 - t21 * t38 + t31 * t26 - t6 * t58, t19 * t41, -t41 * t18 - t19 * t40, t19 * t47 + t41 * t25, -t18 * t47 - t40 * t25, t25 * t47, t1 * t47 + t16 * t18 + t8 * t25 + t3 * t40, t16 * t19 - t2 * t47 - t9 * t25 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t103 * t165 + 0.2e1 * t77 * t105, -0.2e1 * t78 * t105 - 0.2e1 * t99 * t165, 0.2e1 * (t103 * t78 - t77 * t99) * t101, t94 * pkin(2) ^ 2 + t77 ^ 2 + t78 ^ 2, t59 ^ 2, t58 * t168, t75 * t168, t75 * t169, t75 ^ 2, -0.2e1 * t33 * t75 + 0.2e1 * t42 * t58, 0.2e1 * t34 * t75 + 0.2e1 * t42 * t59, t48 ^ 2, t48 * t170, t48 * t169, t58 * t170, t58 ^ 2, 0.2e1 * t20 * t58 + 0.2e1 * t31 * t47, -0.2e1 * t21 * t58 + 0.2e1 * t31 * t48, t41 ^ 2, -0.2e1 * t41 * t40, 0.2e1 * t41 * t47, t40 * t170, t47 ^ 2, 0.2e1 * t16 * t40 + 0.2e1 * t8 * t47, 0.2e1 * t16 * t41 - 0.2e1 * t9 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t57, 0, t43, 0, 0, 0, 0, 0, t104 * t38 - t122 * t45, t104 * t39 + t45 * t139, 0, 0, 0, 0, 0, -t122 * t25 - t79 * t38, -t122 * t26 - t80 * t38, 0, 0, 0, 0, 0, t79 * t18 + t64 * t25, t79 * t19 - t65 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t155, 0, -t158, 0, 0, 0, 0, 0, t104 * t58 - t122 * t75, t104 * t59 + t75 * t139, 0, 0, 0, 0, 0, -t122 * t47 - t79 * t58, -t122 * t48 - t80 * t58, 0, 0, 0, 0, 0, t79 * t40 + t64 * t47, t79 * t41 - t65 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, -t45, t13, -t14, t26 * t108, -t108 * t25 + t26 * t112, t151, t147, 0, -pkin(4) * t25 - pkin(12) * t151 - t10 * t112, -pkin(4) * t26 - pkin(12) * t147 + t10 * t108, t19 * t127 (-t111 * t18 - t142) * t108, -t19 * t112 + t127 * t25, t18 * t112 - t130 * t25, -t25 * t112, -t1 * t112 + t71 * t25 + (pkin(12) * t18 + t157) * t108, t2 * t112 - t72 * t25 + (pkin(12) * t19 + t156) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, -t75, t33, -t34, t48 * t108, -t108 * t47 + t48 * t112, t150, t146, 0, -pkin(4) * t47 - pkin(12) * t150 - t31 * t112, -pkin(4) * t48 - pkin(12) * t146 + t31 * t108, t41 * t127 (-t111 * t40 - t141) * t108, -t41 * t112 + t127 * t47, t40 * t112 - t130 * t47, -t47 * t112, -t8 * t112 + t71 * t47 + (pkin(12) * t40 + t144) * t108, t9 * t112 - t72 * t47 + (pkin(12) * t41 + t143) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t139, 0, 0, 0, 0, 0, t112 * t122, -t108 * t122, 0, 0, 0, 0, 0, -t64 * t112 + t130 * t79, t65 * t112 + t127 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t97, t124, 0, 0, 0, pkin(4) * t166, pkin(4) * t167, t98 * t97, -0.2e1 * t97 * t129, t126 * t167, t107 * t124, t112 ^ 2, -0.2e1 * t71 * t112 + 0.2e1 * t97 * t159, 0.2e1 * t97 * pkin(12) * t111 + 0.2e1 * t72 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t38, t5, -t6, t142, -t107 * t18 + t19 * t111, t153, t149, 0, -pkin(5) * t18 - pkin(13) * t153 - t156, -pkin(5) * t19 - pkin(13) * t149 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, t58, t20, -t21, t141, -t107 * t40 + t41 * t111, t152, t148, 0, -pkin(5) * t40 - pkin(13) * t152 - t143, -pkin(5) * t41 - pkin(13) * t148 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t80, 0, 0, 0, 0, 0, -t79 * t111, t79 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t112, 0, -t108 * pkin(12), -t112 * pkin(12), t107 * t127 (-t96 + t98) * t108, -t128, -t126, 0, -pkin(12) * t127 + (-pkin(5) * t108 + pkin(13) * t112) * t107, pkin(13) * t126 + (t159 - t162) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t96, 0.2e1 * t129, 0, 0, 0, 0.2e1 * t162, -0.2e1 * pkin(5) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, t47, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t130, -t112, t71, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t111, 0, -t107 * pkin(13), -t111 * pkin(13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t23;
