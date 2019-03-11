% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:29
% EndTime: 2019-03-09 13:14:34
% DurationCPUTime: 1.60s
% Computational Cost: add. (5244->177), mult. (11333->305), div. (0->0), fcn. (12130->10), ass. (0->123)
t89 = sin(pkin(11));
t90 = cos(pkin(11));
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t69 = -t89 * t94 + t90 * t97;
t70 = t89 * t97 + t90 * t94;
t93 = sin(qJ(4));
t96 = cos(qJ(4));
t47 = t93 * t69 + t96 * t70;
t95 = cos(qJ(6));
t88 = t95 ^ 2;
t91 = sin(qJ(6));
t138 = t91 ^ 2 - t88;
t119 = t138 * qJD(6);
t141 = -qJ(3) - pkin(7);
t77 = t141 * t94;
t78 = t141 * t97;
t48 = t90 * t77 + t89 * t78;
t44 = -t70 * pkin(8) + t48;
t49 = t89 * t77 - t90 * t78;
t45 = t69 * pkin(8) + t49;
t112 = t96 * t44 - t93 * t45;
t120 = qJD(2) * t141;
t65 = t97 * qJD(3) + t94 * t120;
t66 = -t94 * qJD(3) + t97 * t120;
t42 = -t89 * t65 + t90 * t66;
t136 = qJD(2) * t97;
t137 = qJD(2) * t94;
t68 = t90 * t136 - t89 * t137;
t35 = -t68 * pkin(8) + t42;
t43 = t90 * t65 + t89 * t66;
t67 = t70 * qJD(2);
t36 = -t67 * pkin(8) + t43;
t20 = -t112 * qJD(4) - t93 * t35 - t96 * t36;
t152 = pkin(2) * t89;
t150 = cos(qJ(5));
t24 = -t47 * pkin(9) + t112;
t110 = t96 * t69 - t93 * t70;
t111 = -t93 * t44 - t96 * t45;
t25 = t110 * pkin(9) - t111;
t92 = sin(qJ(5));
t17 = t150 * t25 + t92 * t24;
t109 = -qJD(4) * t47 - t96 * t67 - t93 * t68;
t98 = t109 * pkin(9) - t20;
t21 = t111 * qJD(4) + t96 * t35 - t93 * t36;
t32 = t110 * qJD(4) - t93 * t67 + t96 * t68;
t99 = t32 * pkin(9) - t21;
t5 = t17 * qJD(5) + t150 * t99 + t92 * t98;
t3 = t5 * t91;
t16 = -t150 * t24 + t92 * t25;
t86 = qJD(6) * t95;
t151 = t16 * t86 + t3;
t33 = -t150 * t110 + t92 * t47;
t18 = -t33 * qJD(5) + t92 * t109 + t150 * t32;
t34 = t92 * t110 + t150 * t47;
t149 = t34 * t18;
t148 = t34 * t95;
t19 = t34 * qJD(5) - t150 * t109 + t92 * t32;
t147 = t91 * t19;
t82 = t90 * pkin(2) + pkin(3);
t64 = t96 * t152 + t93 * t82;
t146 = t92 * t64;
t144 = t95 * t18;
t143 = t95 * t19;
t125 = t150 * t64;
t107 = t93 * t152 - t96 * t82;
t61 = pkin(4) - t107;
t103 = t92 * t61 + t125;
t102 = t64 * qJD(4);
t59 = t107 * qJD(4);
t122 = t150 * t102 - t92 * t59;
t29 = t103 * qJD(5) + t122;
t40 = -t150 * t61 - pkin(5) + t146;
t140 = t29 * t91 + t40 * t86;
t135 = qJD(5) * t92;
t129 = pkin(4) * t135;
t128 = t150 * pkin(4);
t84 = -t128 - pkin(5);
t139 = t91 * t129 + t84 * t86;
t134 = qJD(6) * t91;
t133 = -0.2e1 * pkin(1) * qJD(2);
t121 = qJD(5) * t150;
t132 = t92 * t102 - t61 * t121 + t150 * t59;
t131 = pkin(5) * t134;
t130 = pkin(5) * t86;
t85 = pkin(2) * t137;
t127 = t91 * t86;
t126 = -t97 * pkin(2) - pkin(1);
t50 = t67 * pkin(3) + t85;
t124 = -0.4e1 * t91 * t148;
t38 = t40 * t134;
t123 = -t29 * t95 + t38;
t118 = pkin(4) * t121;
t55 = -t69 * pkin(3) + t126;
t37 = -t110 * pkin(4) + t55;
t22 = t33 * pkin(5) - t34 * pkin(10) + t37;
t117 = t95 * t17 + t91 * t22;
t116 = t91 * t17 - t95 * t22;
t41 = pkin(10) + t103;
t115 = t33 * t41 - t34 * t40;
t83 = t92 * pkin(4) + pkin(10);
t114 = t33 * t83 - t34 * t84;
t75 = t84 * t134;
t108 = -t95 * t129 + t75;
t106 = t91 * t18 + t34 * t86;
t105 = t34 * t134 - t144;
t10 = t33 * t86 + t147;
t104 = t33 * t134 - t143;
t28 = t64 * t135 + t132;
t101 = t18 * t40 - t19 * t41 + t28 * t33 + t29 * t34;
t26 = -t109 * pkin(4) + t50;
t100 = t18 * t84 - t19 * t83 + (-t150 * t33 + t34 * t92) * qJD(5) * pkin(4);
t80 = 0.2e1 * t127;
t74 = -0.2e1 * t119;
t31 = t34 ^ 2;
t14 = t16 * t134;
t8 = -t34 * t119 + t91 * t144;
t7 = qJD(6) * t124 - t138 * t18;
t6 = t19 * pkin(5) - t18 * pkin(10) + t26;
t4 = -t24 * t121 + t25 * t135 - t150 * t98 + t92 * t99;
t2 = -t117 * qJD(6) + t91 * t4 + t95 * t6;
t1 = t116 * qJD(6) + t95 * t4 - t91 * t6;
t9 = [0, 0, 0, 0.2e1 * t94 * t136, 0.2e1 * (-t94 ^ 2 + t97 ^ 2) * qJD(2), 0, 0, 0, t94 * t133, t97 * t133, -0.2e1 * t42 * t70 + 0.2e1 * t43 * t69 - 0.2e1 * t48 * t68 - 0.2e1 * t49 * t67, 0.2e1 * t126 * t85 + 0.2e1 * t48 * t42 + 0.2e1 * t49 * t43, 0.2e1 * t47 * t32, 0.2e1 * t47 * t109 + 0.2e1 * t32 * t110, 0, 0, 0, -0.2e1 * t55 * t109 - 0.2e1 * t50 * t110, 0.2e1 * t55 * t32 + 0.2e1 * t50 * t47, 0.2e1 * t149, -0.2e1 * t18 * t33 - 0.2e1 * t34 * t19, 0, 0, 0, 0.2e1 * t37 * t19 + 0.2e1 * t26 * t33, 0.2e1 * t37 * t18 + 0.2e1 * t26 * t34, -0.2e1 * t31 * t127 + 0.2e1 * t88 * t149, 0.2e1 * t31 * t119 + t18 * t124, -0.2e1 * t105 * t33 + 0.2e1 * t34 * t143, -0.2e1 * t106 * t33 - 0.2e1 * t34 * t147, 0.2e1 * t33 * t19, 0.2e1 * t106 * t16 - 0.2e1 * t116 * t19 + 0.2e1 * t2 * t33 + 0.2e1 * t34 * t3, 0.2e1 * t1 * t33 - 0.2e1 * t105 * t16 - 0.2e1 * t117 * t19 + 0.2e1 * t148 * t5; 0, 0, 0, 0, 0, t136, -t137, 0, -pkin(7) * t136, pkin(7) * t137 (-t67 * t89 - t68 * t90) * pkin(2) (t42 * t90 + t43 * t89) * pkin(2), 0, 0, t32, t109, 0, t21, t20, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t104, 0, t14 + (-t115 * qJD(6) - t5) * t95 + t101 * t91, t101 * t95 + t115 * t134 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t102, 0.2e1 * t59, 0, 0, 0, 0, 0, -0.2e1 * t29, 0.2e1 * t28, t80, t74, 0, 0, 0, 0.2e1 * t123, 0.2e1 * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, -t109, t32, 0, 0, 0, 0, 0, t19, t18, 0, 0, 0, 0, 0, -t104, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t109, 0, t21, t20, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t104, 0, t14 + (-t114 * qJD(6) - t5) * t95 + t100 * t91, t100 * t95 + t114 * t134 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t59, 0, 0, 0, 0, 0 (-t125 + (-pkin(4) - t61) * t92) * qJD(5) - t122 (-t128 + t146) * qJD(5) + t132, t80, t74, 0, 0, 0, t38 + t75 + (-t29 - t129) * t95, t139 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t129, -0.2e1 * t118, t80, t74, 0, 0, 0, 0.2e1 * t108, 0.2e1 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, -t5, t4, t8, t7, t10, -t104, 0, t14 + (-pkin(5) * t18 - pkin(10) * t19) * t91 + (-t5 + (-pkin(5) * t34 - pkin(10) * t33) * qJD(6)) * t95, pkin(5) * t105 + pkin(10) * t104 + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, t80, t74, 0, 0, 0, t123 - t131, -t130 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t118, t80, t74, 0, 0, 0, t108 - t131, -t130 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t74, 0, 0, 0, -0.2e1 * t131, -0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t106, t19, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t134, 0, t91 * t28 - t41 * t86, t134 * t41 + t95 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t134, 0, -t91 * t118 - t83 * t86, -t118 * t95 + t134 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t134, 0, -pkin(10) * t86, pkin(10) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
