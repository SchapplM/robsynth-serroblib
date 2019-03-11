% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:03
% EndTime: 2019-03-09 08:48:07
% DurationCPUTime: 1.10s
% Computational Cost: add. (1700->149), mult. (3749->278), div. (0->0), fcn. (3704->8), ass. (0->104)
t76 = cos(qJ(6));
t72 = t76 ^ 2;
t73 = sin(qJ(6));
t117 = t73 ^ 2 - t72;
t132 = qJD(6) * t117;
t134 = 0.2e1 * t132;
t115 = sin(pkin(10));
t116 = cos(pkin(10));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t86 = -t115 * t77 - t116 * t75;
t133 = qJD(2) * t86;
t74 = sin(qJ(5));
t126 = cos(qJ(5));
t103 = t115 * t75;
t64 = t116 * t77;
t94 = t103 - t64;
t88 = t126 * t94;
t30 = -t74 * t86 - t88;
t31 = -t126 * t86 + t74 * t94;
t63 = t115 * pkin(2) + qJ(4);
t65 = -t116 * pkin(2) - pkin(3);
t93 = -pkin(4) + t65;
t85 = t126 * t93;
t36 = t74 * t63 + pkin(5) - t85;
t82 = t126 * t63 + t74 * t93;
t37 = -pkin(9) + t82;
t106 = qJD(5) * t126;
t114 = qJD(5) * t74;
t118 = -qJ(3) - pkin(7);
t59 = t118 * t75;
t60 = t118 * t77;
t34 = -t115 * t60 - t116 * t59;
t26 = pkin(8) * t86 + t34;
t35 = t115 * t59 - t116 * t60;
t27 = t94 * pkin(8) + t35;
t102 = qJD(2) * t118;
t45 = t77 * qJD(3) + t75 * t102;
t87 = -t75 * qJD(3) + t77 * t102;
t25 = t115 * t87 + t116 * t45;
t78 = -pkin(8) * t133 + t25;
t24 = t115 * t45 - t116 * t87;
t62 = qJD(2) * t103;
t89 = qJD(2) * t64 - t62;
t80 = -t89 * pkin(8) + t24;
t5 = t27 * t106 + t26 * t114 - t126 * t80 + t74 * t78;
t131 = -t5 + (t30 * t37 - t31 * t36) * qJD(6);
t130 = 0.2e1 * t133 * t35 - 0.2e1 * t24 * t86 - 0.2e1 * t25 * t94 + 0.2e1 * t34 * t89;
t129 = 2 * qJD(4);
t128 = t5 * t73;
t127 = -t77 * pkin(2) - pkin(1);
t13 = -qJD(5) * t88 - t114 * t86 - t126 * t89 + t133 * t74;
t125 = t31 * t13;
t124 = t31 * t76;
t33 = t74 * qJD(4) + t82 * qJD(5);
t123 = t33 * t73;
t122 = t33 * t76;
t14 = t31 * qJD(5) + t126 * t133 + t74 * t89;
t121 = t73 * t14;
t120 = t76 * t13;
t119 = t76 * t14;
t68 = qJD(6) * t73;
t69 = qJD(6) * t76;
t113 = t75 * qJD(2);
t112 = t77 * qJD(2);
t111 = -0.2e1 * pkin(1) * qJD(2);
t110 = -0.2e1 * pkin(5) * qJD(6);
t66 = pkin(2) * t113;
t109 = t73 * t69;
t108 = 0.4e1 * t73 * t124;
t107 = t34 * t24 + t35 * t25;
t105 = qJD(6) * (pkin(5) + t36);
t104 = qJD(6) * t126;
t28 = t94 * pkin(3) + qJ(4) * t86 + t127;
t19 = -pkin(3) * t133 - t89 * qJ(4) + qJD(4) * t86 + t66;
t99 = pkin(5) * t13 - pkin(9) * t14;
t98 = pkin(5) * t31 + pkin(9) * t30;
t12 = t126 * t27 + t74 * t26;
t21 = -t94 * pkin(4) - t28;
t8 = t30 * pkin(5) - t31 * pkin(9) + t21;
t97 = t76 * t12 + t73 * t8;
t96 = t73 * t12 - t76 * t8;
t92 = -t73 * t13 + t31 * t69;
t91 = -t31 * t68 - t120;
t10 = t30 * t69 + t121;
t90 = t126 * t31 + t30 * t74;
t11 = -t126 * t26 + t74 * t27;
t32 = -t126 * qJD(4) - qJD(5) * t85 + t63 * t114;
t81 = -qJD(6) * t11 - t13 * t36 - t14 * t37 + t30 * t32 + t31 * t33;
t79 = t126 * t13 - t14 * t74 + (-t126 * t30 + t31 * t74) * qJD(5);
t15 = pkin(4) * t133 - t19;
t61 = 0.2e1 * t109;
t56 = -0.2e1 * t132;
t48 = t73 * t104 + t76 * t114;
t47 = -t76 * t104 + t73 * t114;
t29 = t31 ^ 2;
t9 = t30 * t68 - t119;
t7 = t73 * t120 + t132 * t31;
t6 = qJD(6) * t108 - t117 * t13;
t4 = -t26 * t106 + t27 * t114 - t126 * t78 - t74 * t80;
t3 = t14 * pkin(5) + t13 * pkin(9) + t15;
t2 = -t97 * qJD(6) + t76 * t3 + t73 * t4;
t1 = t96 * qJD(6) - t73 * t3 + t76 * t4;
t16 = [0, 0, 0, 0.2e1 * t75 * t112, 0.2e1 * (-t75 ^ 2 + t77 ^ 2) * qJD(2), 0, 0, 0, t75 * t111, t77 * t111, t130, 0.2e1 * t127 * t66 + 0.2e1 * t107, -0.2e1 * t133 * t28 + 0.2e1 * t19 * t94, t130, 0.2e1 * t19 * t86 - 0.2e1 * t28 * t89, 0.2e1 * t28 * t19 + 0.2e1 * t107, -0.2e1 * t125, 0.2e1 * t13 * t30 - 0.2e1 * t31 * t14, 0, 0, 0, 0.2e1 * t21 * t14 + 0.2e1 * t15 * t30, -0.2e1 * t21 * t13 + 0.2e1 * t15 * t31, -0.2e1 * t109 * t29 - 0.2e1 * t125 * t72, t13 * t108 + t29 * t134, 0.2e1 * t119 * t31 + 0.2e1 * t30 * t91, -0.2e1 * t121 * t31 - 0.2e1 * t30 * t92, 0.2e1 * t30 * t14, 0.2e1 * t11 * t92 + 0.2e1 * t31 * t128 - 0.2e1 * t14 * t96 + 0.2e1 * t2 * t30, 0.2e1 * t1 * t30 + 0.2e1 * t11 * t91 + 0.2e1 * t124 * t5 - 0.2e1 * t14 * t97; 0, 0, 0, 0, 0, t112, -t113, 0, -pkin(7) * t112, pkin(7) * t113 (t116 * t62 + (-t116 ^ 2 * t77 + t115 * t86) * qJD(2)) * pkin(2) (t115 * t25 - t116 * t24) * pkin(2), -t24, -qJD(4) * t94 + t133 * t63 + t65 * t89, t25, t35 * qJD(4) + t24 * t65 + t25 * t63, 0, 0, t13, t14, 0, t5, -t4, t7, t6, -t10, t9, 0, -t131 * t76 + t81 * t73, t131 * t73 + t81 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t63 * t129, 0, 0, 0, 0, 0, 0.2e1 * t33, -0.2e1 * t32, t61, t56, 0, 0, 0, -0.2e1 * t36 * t68 + 0.2e1 * t122, -0.2e1 * t36 * t69 - 0.2e1 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t133, 0, -t89, t19, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t90 + t73 * t79, t68 * t90 + t76 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t106, 0, 0, 0, 0, 0, t48, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t5, t4, -t7, -t6, t10, -t9, 0, -t5 * t76 + t99 * t73 + (t11 * t73 - t76 * t98) * qJD(6), t128 + t99 * t76 + (t11 * t76 + t73 * t98) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, -0.2e1 * t109, t134, 0, 0, 0, t105 * t73 - t122, t105 * t76 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t106, 0, 0, 0, 0, 0, -t48, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t56, 0, 0, 0, t73 * t110, t76 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t92, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t68, 0, t73 * t32 - t37 * t69, t76 * t32 + t37 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t73 - t69 * t74, -t106 * t76 + t68 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t68, 0, -pkin(9) * t69, pkin(9) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t16;
