% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:38
% EndTime: 2019-03-09 02:54:42
% DurationCPUTime: 1.13s
% Computational Cost: add. (1388->138), mult. (2966->249), div. (0->0), fcn. (2802->8), ass. (0->90)
t65 = sin(pkin(10));
t66 = cos(pkin(10));
t67 = sin(qJ(6));
t69 = cos(qJ(6));
t124 = -t67 * t65 + t69 * t66;
t122 = t124 * qJD(6);
t102 = cos(pkin(9));
t68 = sin(qJ(3));
t101 = sin(pkin(9));
t70 = cos(qJ(3));
t92 = t101 * t70;
t45 = -t102 * t68 - t92;
t90 = qJD(3) * t101;
t91 = qJD(3) * t102;
t37 = t68 * t90 - t70 * t91;
t47 = t69 * t65 + t67 * t66;
t78 = t47 * t37;
t126 = t122 * t45 + t78;
t38 = -t68 * t91 - t70 * t90;
t93 = t102 * t70;
t44 = -t101 * t68 + t93;
t11 = t122 * t44 + t38 * t47;
t110 = t124 * t37;
t121 = t47 * qJD(6);
t18 = t121 * t45 - t110;
t9 = -t121 * t44 + t124 * t38;
t112 = t44 * t38;
t116 = t37 * t45;
t76 = 0.2e1 * t112 + 0.2e1 * t116;
t123 = (-t101 * t37 + t102 * t38) * pkin(3);
t120 = 0.2e1 * t122;
t119 = 2 * qJD(2);
t55 = t101 * pkin(3) + qJ(5);
t118 = pkin(8) + t55;
t71 = -pkin(1) - pkin(7);
t103 = qJ(4) - t71;
t98 = t70 * qJD(3);
t36 = -t68 * qJD(4) - t103 * t98;
t99 = t68 * qJD(3);
t74 = -t70 * qJD(4) + t103 * t99;
t24 = t101 * t36 - t102 * t74;
t48 = t103 * t68;
t32 = -t101 * t48 + t103 * t93;
t117 = t32 * t24;
t60 = -t102 * pkin(3) - pkin(4);
t113 = t38 * t60;
t111 = t44 * t65;
t109 = t65 * t38;
t108 = t66 * t38;
t50 = pkin(3) * t98 + qJD(2);
t19 = -t37 * pkin(4) - t38 * qJ(5) - t44 * qJD(5) + t50;
t25 = t101 * t74 + t102 * t36;
t6 = t65 * t19 + t66 * t25;
t104 = t68 * pkin(3) + qJ(2);
t30 = -t45 * pkin(4) - t44 * qJ(5) + t104;
t33 = -t102 * t48 - t103 * t92;
t14 = t65 * t30 + t66 * t33;
t105 = t65 ^ 2 + t66 ^ 2;
t100 = qJD(5) * t45;
t97 = qJ(2) * qJD(3);
t94 = t105 * t37;
t5 = t66 * t19 - t65 * t25;
t13 = t66 * t30 - t65 * t33;
t89 = 0.2e1 * t105 * qJD(5);
t88 = t5 * t66 + t6 * t65;
t87 = t5 * t65 - t6 * t66;
t7 = -t66 * t44 * pkin(8) - t45 * pkin(5) + t13;
t8 = -pkin(8) * t111 + t14;
t86 = t67 * t8 - t69 * t7;
t85 = t67 * t7 + t69 * t8;
t84 = t13 * t65 - t14 * t66;
t83 = t24 * t44 + t32 * t38;
t41 = t118 * t65;
t42 = t118 * t66;
t80 = -t69 * t41 - t67 * t42;
t79 = -t67 * t41 + t69 * t42;
t75 = t37 * t55 + t100 + t113;
t72 = t25 * t45 + t33 * t37 + t83;
t49 = -t66 * pkin(5) + t60;
t27 = t124 * t44;
t26 = t47 * t44;
t23 = pkin(5) * t111 + t32;
t21 = -t47 * qJD(5) - t79 * qJD(6);
t20 = -qJD(5) * t124 - t80 * qJD(6);
t15 = pkin(5) * t109 + t24;
t4 = -pkin(8) * t109 + t6;
t3 = -t37 * pkin(5) - pkin(8) * t108 + t5;
t2 = -t85 * qJD(6) + t69 * t3 - t67 * t4;
t1 = t86 * qJD(6) - t67 * t3 - t69 * t4;
t10 = [0, 0, 0, 0, t119, qJ(2) * t119, -0.2e1 * t68 * t98, 0.2e1 * (t68 ^ 2 - t70 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t68 + 0.2e1 * t70 * t97, 0.2e1 * qJD(2) * t70 - 0.2e1 * t68 * t97, 0.2e1 * t72, 0.2e1 * t104 * t50 + 0.2e1 * t33 * t25 + 0.2e1 * t117, -0.2e1 * t13 * t37 - 0.2e1 * t5 * t45 + 0.2e1 * t83 * t65, 0.2e1 * t14 * t37 + 0.2e1 * t6 * t45 + 0.2e1 * t83 * t66, -0.2e1 * t88 * t44 + 0.2e1 * (-t13 * t66 - t14 * t65) * t38, 0.2e1 * t13 * t5 + 0.2e1 * t14 * t6 + 0.2e1 * t117, 0.2e1 * t27 * t9, -0.2e1 * t27 * t11 - 0.2e1 * t9 * t26, -0.2e1 * t27 * t37 - 0.2e1 * t9 * t45, 0.2e1 * t11 * t45 + 0.2e1 * t26 * t37, 0.2e1 * t116, 0.2e1 * t23 * t11 + 0.2e1 * t15 * t26 - 0.2e1 * t2 * t45 + 0.2e1 * t86 * t37, -0.2e1 * t1 * t45 + 0.2e1 * t15 * t27 + 0.2e1 * t23 * t9 + 0.2e1 * t85 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t72, -t65 * t76, -t66 * t76, 0, t84 * t37 + t87 * t45 - t83, 0, 0, 0, 0, 0, -t44 * t11 - t38 * t26 + (-t126 - t78) * t45, -t38 * t27 - t44 * t9 + (t18 - t110) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0.2e1 * t105 * t116 + 0.2e1 * t112, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98, 0, -t71 * t99, -t71 * t98, -t123 (t101 * t25 - t102 * t24) * pkin(3), -t24 * t66 + t75 * t65, t24 * t65 + t75 * t66, -t87, -t84 * qJD(5) + t24 * t60 - t87 * t55, t122 * t27 + t9 * t47, -t47 * t11 - t121 * t27 - t122 * t26 + t124 * t9, -t126, t18, 0, t49 * t11 + t121 * t23 - t124 * t15 - t21 * t45 - t80 * t37, t122 * t23 + t15 * t47 - t20 * t45 + t79 * t37 + t49 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98, 0, t123, t108, -t109, -t94, -t105 * t100 - t55 * t94 - t113, 0, 0, 0, 0, 0, t9, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t55 * t89, t47 * t120, -0.2e1 * t121 * t47 + 0.2e1 * t122 * t124, 0, 0, 0, 0.2e1 * t49 * t121, t49 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t66 * t37, t65 * t37, -t105 * t38, t88, 0, 0, 0, 0, 0, t18, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t108, 0, t24, 0, 0, 0, 0, 0, t11, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t11, -t37, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t121, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
