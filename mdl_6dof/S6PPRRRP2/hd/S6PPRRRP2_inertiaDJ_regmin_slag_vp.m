% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:14
% EndTime: 2019-03-08 18:58:18
% DurationCPUTime: 1.41s
% Computational Cost: add. (1460->200), mult. (4457->354), div. (0->0), fcn. (4608->12), ass. (0->120)
t71 = sin(qJ(4));
t153 = -0.4e1 * t71;
t144 = t71 * pkin(10);
t73 = cos(qJ(4));
t99 = -t73 * pkin(4) - t144;
t56 = -pkin(3) + t99;
t146 = pkin(9) * t73;
t72 = cos(qJ(5));
t59 = t72 * t146;
t70 = sin(qJ(5));
t152 = t70 * t56 + t59;
t68 = cos(pkin(7));
t135 = sin(pkin(7));
t143 = sin(qJ(3));
t94 = t135 * t143;
t42 = -t73 * t68 + t71 * t94;
t74 = cos(qJ(3));
t105 = t74 * t135;
t96 = qJD(3) * t105;
t75 = -qJD(4) * t42 + t73 * t96;
t63 = t72 ^ 2;
t137 = t70 ^ 2 - t63;
t104 = t137 * qJD(5);
t145 = pkin(10) * t73;
t98 = pkin(4) * t71 - t145;
t50 = t98 * qJD(4);
t60 = qJD(5) * t72;
t128 = t71 * qJD(4);
t109 = t72 * t128;
t129 = qJD(5) * t73;
t117 = t70 * t129;
t83 = t109 + t117;
t23 = t83 * pkin(9) - t70 * t50 - t56 * t60;
t67 = cos(pkin(12));
t142 = t67 * t68;
t65 = sin(pkin(12));
t66 = sin(pkin(6));
t151 = (-t74 * t142 + t143 * t65) * t66;
t92 = pkin(5) * t70 - qJ(6) * t72;
t85 = pkin(9) + t92;
t39 = t85 * t71;
t40 = t92 * qJD(5) - t70 * qJD(6);
t93 = t72 * pkin(5) + t70 * qJ(6);
t52 = -pkin(4) - t93;
t150 = qJD(4) * (-t52 * t73 + t144) - qJD(5) * t39 - t40 * t71;
t149 = t93 * qJD(5) - t72 * qJD(6);
t113 = t70 * t128;
t24 = pkin(9) * t113 - qJD(5) * t152 + t72 * t50;
t148 = 0.2e1 * qJD(6);
t147 = pkin(9) * t70;
t141 = t70 * t71;
t62 = t71 ^ 2;
t136 = -t73 ^ 2 + t62;
t132 = qJD(4) * t72;
t131 = qJD(5) * t70;
t130 = qJD(5) * t71;
t126 = t73 * qJD(4);
t125 = t73 * qJD(6);
t124 = qJ(6) * qJD(4);
t123 = -0.2e1 * pkin(3) * qJD(4);
t122 = -0.2e1 * pkin(4) * qJD(5);
t120 = pkin(5) * t128;
t119 = pkin(10) * t131;
t118 = pkin(10) * t60;
t116 = t72 * t129;
t69 = cos(pkin(6));
t31 = t69 * t94 + (t143 * t142 + t65 * t74) * t66;
t41 = -t66 * t67 * t135 + t69 * t68;
t17 = t31 * t71 - t41 * t73;
t115 = t17 * t60;
t114 = t42 * t60;
t112 = t70 * t126;
t111 = t70 * t60;
t110 = t71 * t126;
t108 = t72 * t126;
t106 = t71 * t124;
t103 = t136 * qJD(4);
t102 = 0.2e1 * t110;
t101 = t70 * t108;
t100 = t70 * t105;
t97 = qJD(4) * t105;
t18 = t31 * t73 + t41 * t71;
t30 = -t69 * t105 + t151;
t10 = t18 * t70 - t30 * t72;
t11 = t18 * t72 + t30 * t70;
t91 = t10 * t72 - t11 * t70;
t43 = t71 * t68 + t73 * t94;
t33 = t72 * t105 + t70 * t43;
t34 = t72 * t43 - t100;
t90 = t33 * t72 - t34 * t70;
t35 = -t73 * qJ(6) + t152;
t36 = -t72 * t56 + (pkin(5) + t147) * t73;
t89 = -t35 * t70 + t36 * t72;
t26 = qJD(3) * t151 - t69 * t96;
t8 = t18 * qJD(4) - t26 * t71;
t6 = t17 * t131 - t8 * t72;
t84 = qJD(3) * t94;
t32 = t43 * qJD(4) + t71 * t96;
t20 = t42 * t131 - t32 * t72;
t27 = t31 * qJD(3);
t9 = -qJD(4) * t17 - t26 * t73;
t2 = t11 * qJD(5) - t27 * t72 + t9 * t70;
t82 = -t10 * t128 + t17 * t112 + t71 * t115 + t8 * t141 + t2 * t73;
t13 = -qJD(5) * t100 + t43 * t60 + t70 * t75 - t72 * t84;
t81 = t42 * t112 + t71 * t114 - t33 * t128 + t13 * t73 + t32 * t141;
t3 = -t10 * qJD(5) + t27 * t70 + t9 * t72;
t78 = t91 * qJD(5) + t2 * t70 + t3 * t72;
t12 = t33 * qJD(5) - t70 * t84 - t72 * t75;
t77 = t90 * qJD(5) - t12 * t72 + t13 * t70;
t16 = t106 - t23 - t125;
t21 = -t120 - t24;
t76 = t89 * qJD(5) + t16 * t72 + t21 * t70;
t58 = pkin(10) * t116;
t44 = -t70 * t130 + t108;
t22 = t85 * t126 + t149 * t71;
t19 = t32 * t70 + t114;
t5 = t8 * t70 + t115;
t4 = (t42 * t132 - t12) * t73 + (-qJD(4) * t34 - t20) * t71;
t1 = (t17 * t132 + t3) * t73 + (-qJD(4) * t11 - t6) * t71;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t2 + 0.2e1 * t11 * t3 + 0.2e1 * t17 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t13 - t11 * t12 + t17 * t32 + t2 * t33 + t3 * t34 + t8 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t34 * t12 + 0.2e1 * t33 * t13 + 0.2e1 * t42 * t32; 0, 0, 0, -t27, t26, 0, 0, 0, 0, 0, t30 * t128 - t27 * t73, t30 * t126 + t27 * t71, 0, 0, 0, 0, 0, t82, t1, t82, t91 * t126 + (t2 * t72 - t3 * t70 + (-t10 * t70 - t11 * t72) * qJD(5)) * t71, -t1, t10 * t21 + t11 * t16 + t17 * t22 + t2 * t36 + t3 * t35 + t8 * t39; 0, 0, 0, -t84, -t96, 0, 0, 0, 0, 0, -t71 * t97 - t73 * t84, t71 * t84 - t73 * t97, 0, 0, 0, 0, 0, t81, t4, t81, t90 * t126 + (t12 * t70 + t13 * t72 + (-t33 * t70 - t34 * t72) * qJD(5)) * t71, -t4, -t12 * t35 + t13 * t36 + t34 * t16 + t33 * t21 + t42 * t22 + t32 * t39; 0, 0, 0, 0, 0, t102, -0.2e1 * t103, 0, 0, 0, t71 * t123, t73 * t123, 0.2e1 * t63 * t110 - 0.2e1 * t62 * t111, t101 * t153 + 0.2e1 * t62 * t104, 0.2e1 * t71 * t117 + 0.2e1 * t136 * t132, -0.2e1 * t70 * t103 + 0.2e1 * t71 * t116, -0.2e1 * t110, 0.2e1 * t56 * t109 - 0.2e1 * t24 * t73 + 0.2e1 * (t70 * t110 + t62 * t60) * pkin(9), -0.2e1 * t23 * t73 - 0.2e1 * t152 * t128 + 0.2e1 * (t102 * t72 - t62 * t131) * pkin(9), 0.2e1 * (qJD(4) * t39 * t70 + t21) * t73 + 0.2e1 * (-qJD(4) * t36 + t22 * t70 + t39 * t60) * t71, 0.2e1 * t89 * t126 + 0.2e1 * (-t16 * t70 + t21 * t72 + (-t35 * t72 - t36 * t70) * qJD(5)) * t71, 0.2e1 * (-t39 * t132 - t16) * t73 + 0.2e1 * (qJD(4) * t35 + t39 * t131 - t22 * t72) * t71, 0.2e1 * t35 * t16 + 0.2e1 * t36 * t21 + 0.2e1 * t39 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, 0, 0, 0, 0, t6, t5, t6, t78, -t5, pkin(10) * t78 + t17 * t40 + t8 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t75, 0, 0, 0, 0, 0, t20, t19, t20, t77, -t19, pkin(10) * t77 + t32 * t52 + t42 * t40; 0, 0, 0, 0, 0, 0, 0, t126, -t128, 0, -pkin(9) * t126, pkin(9) * t128, -t71 * t104 + t101, t111 * t153 - t137 * t126, t113 - t116, t83, 0, t58 + (-pkin(4) * t72 + t147) * t130 + (t99 * t70 - t59) * qJD(4) (pkin(9) * t71 * t72 + t70 * t98) * qJD(5) + (t70 * t146 + t72 * t99) * qJD(4), t58 + (t52 * t130 - t22) * t72 - t150 * t70, t76 (-t22 + (t52 * t71 + t145) * qJD(5)) * t70 + t150 * t72, pkin(10) * t76 + t22 * t52 + t39 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111, -0.2e1 * t104, 0, 0, 0, t70 * t122, t72 * t122, 0.2e1 * t52 * t131 - 0.2e1 * t40 * t72, 0, -0.2e1 * t40 * t70 - 0.2e1 * t52 * t60, 0.2e1 * t52 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6) + t11 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, -t13, 0, -t12, -t13 * pkin(5) - t12 * qJ(6) + t34 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t71 * t60 - t112, t128, t24, t23, t24 + 0.2e1 * t120 (-pkin(5) * t126 - qJ(6) * t130) * t72 + (-t73 * t124 + (pkin(5) * qJD(5) - qJD(6)) * t71) * t70, 0.2e1 * t106 - t23 - 0.2e1 * t125, -t21 * pkin(5) + t16 * qJ(6) + t35 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t131, 0, -t118, t119, -t118, -t149, -t119, -t149 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, qJ(6) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t44, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
