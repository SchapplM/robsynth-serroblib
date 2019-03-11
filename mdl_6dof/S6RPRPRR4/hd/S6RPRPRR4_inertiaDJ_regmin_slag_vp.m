% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:37
% EndTime: 2019-03-09 03:46:39
% DurationCPUTime: 1.06s
% Computational Cost: add. (904->154), mult. (2060->268), div. (0->0), fcn. (1721->8), ass. (0->112)
t59 = sin(pkin(10)) * pkin(1) + pkin(7);
t128 = pkin(4) + t59;
t71 = sin(qJ(6));
t72 = sin(qJ(5));
t74 = cos(qJ(6));
t75 = cos(qJ(5));
t45 = t71 * t75 + t74 * t72;
t76 = cos(qJ(3));
t30 = t45 * t76;
t77 = -pkin(3) - pkin(8);
t135 = t76 * t77;
t73 = sin(qJ(3));
t117 = t73 * qJ(4);
t60 = -cos(pkin(10)) * pkin(1) - pkin(2);
t80 = t60 - t117;
t27 = t80 + t135;
t41 = t128 * t73;
t31 = t72 * t41;
t121 = t75 * t27 + t31;
t113 = qJD(5) * t76;
t100 = t75 * t113;
t63 = t73 * qJD(3);
t99 = t72 * t63;
t134 = t99 - t100;
t69 = t76 ^ 2;
t90 = qJD(3) * (t73 ^ 2 - t69);
t66 = t72 ^ 2;
t120 = -t75 ^ 2 + t66;
t89 = t120 * qJD(5);
t107 = qJD(5) + qJD(6);
t109 = t76 * qJD(4);
t42 = t128 * t76;
t110 = t42 * qJD(5);
t133 = (t117 - t135) * qJD(3) - t109 - t110;
t132 = 0.2e1 * qJD(4);
t131 = t73 * pkin(5);
t101 = t72 * t113;
t96 = t75 * t63;
t37 = t96 + t101;
t114 = qJD(5) * t75;
t115 = qJD(5) * t72;
t118 = qJ(4) * t76;
t88 = pkin(3) * t63 - t73 * qJD(4);
t21 = (pkin(8) * t73 - t118) * qJD(3) + t88;
t64 = t76 * qJD(3);
t47 = t59 * t64;
t34 = pkin(4) * t64 + t47;
t6 = -t41 * t114 + t27 * t115 - t75 * t21 - t72 * t34;
t5 = t37 * pkin(9) - t6;
t130 = t74 * t5;
t129 = t76 * pkin(3);
t127 = pkin(9) - t77;
t122 = t75 * t76;
t16 = -pkin(9) * t122 + t121;
t126 = t71 * t16;
t125 = t72 * t76;
t124 = t74 * t16;
t123 = t74 * t75;
t116 = qJD(3) * t42;
t112 = qJD(5) * t77;
t111 = qJD(6) * t71;
t108 = qJ(4) * qJD(5);
t106 = t74 * t122;
t105 = 0.2e1 * qJD(3) * t60;
t104 = pkin(5) * t64;
t103 = pkin(5) * t111;
t102 = qJD(6) * t74 * pkin(5);
t98 = t73 * t64;
t97 = t59 * t63;
t95 = t72 * t114;
t91 = -t72 * t21 + t75 * t34;
t92 = pkin(9) * t76 - t27;
t4 = (-pkin(9) * t72 * t73 + pkin(5) * t76) * qJD(3) + (t92 * t75 - t31) * qJD(5) + t91;
t94 = t74 * t4 - t71 * t5;
t49 = t127 * t75;
t32 = t75 * t41;
t13 = t92 * t72 + t131 + t32;
t93 = -t13 - t131;
t87 = t72 * t96;
t86 = t74 * t13 - t126;
t85 = t71 * t13 + t124;
t48 = t127 * t72;
t84 = -t74 * t48 - t71 * t49;
t83 = -t71 * t48 + t74 * t49;
t19 = t107 * t123 - t72 * t111 - t71 * t115;
t14 = t73 * t19 + t45 * t64;
t33 = t128 * t63;
t79 = -t33 + (-t73 * t77 - t118) * qJD(5);
t78 = t109 + (-t117 - t129) * qJD(3);
t61 = t72 * pkin(5) + qJ(4);
t56 = pkin(5) * t114 + qJD(4);
t55 = 0.2e1 * t98;
t46 = -t71 * t72 + t123;
t44 = qJD(5) * t49;
t43 = t127 * t115;
t40 = t80 - t129;
t38 = t73 * t114 + t72 * t64;
t36 = -t73 * t115 + t75 * t64;
t35 = qJ(4) * t64 - t88;
t29 = t71 * t125 - t106;
t24 = pkin(5) * t122 + t42;
t18 = t107 * t45;
t17 = -pkin(5) * t101 + (-pkin(5) * t75 - t128) * t63;
t15 = -t18 * t73 + t46 * t64;
t12 = t107 * t30 - t71 * t99 + t74 * t96;
t11 = -qJD(6) * t106 + (t107 * t125 + t96) * t71 + t134 * t74;
t9 = -t84 * qJD(6) + t74 * t43 + t71 * t44;
t8 = t83 * qJD(6) - t71 * t43 + t74 * t44;
t7 = -t121 * qJD(5) + t91;
t2 = -t85 * qJD(6) + t94;
t1 = -t86 * qJD(6) - t71 * t4 - t130;
t3 = [0, 0, 0, 0, t55, -0.2e1 * t90, 0, 0, 0, t73 * t105, t76 * t105, 0, -0.2e1 * t35 * t76 - 0.2e1 * t40 * t63, 0.2e1 * t35 * t73 - 0.2e1 * t40 * t64, -0.2e1 * t40 * t35, -0.2e1 * t66 * t98 + 0.2e1 * t69 * t95, -0.2e1 * t69 * t89 - 0.4e1 * t76 * t87, -0.2e1 * t73 * t100 + 0.2e1 * t72 * t90, 0.2e1 * t73 * t101 + 0.2e1 * t75 * t90, t55, 0.2e1 * (-t75 * t116 + t7) * t73 + 0.2e1 * ((-t72 * t27 + t32) * qJD(3) - t33 * t75 - t72 * t110) * t76, 0.2e1 * (t72 * t116 + t6) * t73 + 0.2e1 * (-t121 * qJD(3) - t75 * t110 + t33 * t72) * t76, -0.2e1 * t30 * t11, 0.2e1 * t11 * t29 - 0.2e1 * t30 * t12, 0.2e1 * t11 * t73 - 0.2e1 * t30 * t64, 0.2e1 * t12 * t73 + 0.2e1 * t29 * t64, t55, -0.2e1 * t24 * t12 - 0.2e1 * t17 * t29 + 0.2e1 * t2 * t73 + 0.2e1 * t64 * t86, 0.2e1 * t1 * t73 + 0.2e1 * t24 * t11 - 0.2e1 * t17 * t30 - 0.2e1 * t64 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t64, -t63, 0, -t47, t97, t78, t47, -t97, t78 * t59, t76 * t89 + t87, -t120 * t63 + 0.4e1 * t76 * t95, t36, -t38, 0, -t133 * t75 + t79 * t72, t133 * t72 + t79 * t75, t11 * t46 + t30 * t18, -t11 * t45 + t46 * t12 - t18 * t29 + t30 * t19, t15, -t14, 0, -t61 * t12 + t17 * t45 + t24 * t19 - t56 * t29 - t64 * t83 + t9 * t73, t61 * t11 + t17 * t46 - t24 * t18 - t56 * t30 - t64 * t84 + t8 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t64, 0, t63, t64, t35, 0, 0, 0, 0, 0, t38, t36, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, qJ(4) * t132, -0.2e1 * t95, 0.2e1 * t89, 0, 0, 0, 0.2e1 * qJD(4) * t72 + 0.2e1 * t75 * t108, 0.2e1 * qJD(4) * t75 - 0.2e1 * t108 * t72, -0.2e1 * t46 * t18, 0.2e1 * t18 * t45 - 0.2e1 * t46 * t19, 0, 0, 0, 0.2e1 * t61 * t19 + 0.2e1 * t56 * t45, -0.2e1 * t61 * t18 + 0.2e1 * t56 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, t47, 0, 0, 0, 0, 0, t36, -t38, 0, 0, 0, 0, 0, t15, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t37, t64, t7, t6, 0, 0, t11, t12, t64, t74 * t104 + (t71 * t93 - t124) * qJD(6) + t94, -t130 + (-t4 - t104) * t71 + (t74 * t93 + t126) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t134, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t114, 0, -t72 * t112, -t75 * t112, 0, 0, -t18, -t19, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t114, 0, 0, 0, 0, 0, -t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, t64, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
