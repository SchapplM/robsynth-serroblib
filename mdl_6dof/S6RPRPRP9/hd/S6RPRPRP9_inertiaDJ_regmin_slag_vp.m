% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:27
% EndTime: 2019-03-09 03:29:31
% DurationCPUTime: 1.49s
% Computational Cost: add. (1700->211), mult. (3811->360), div. (0->0), fcn. (3340->6), ass. (0->102)
t125 = sin(qJ(5));
t126 = cos(qJ(5));
t71 = sin(pkin(9));
t72 = cos(pkin(9));
t49 = t125 * t71 - t126 * t72;
t100 = qJD(5) * t125;
t101 = qJD(5) * t126;
t43 = t71 * t100 - t72 * t101;
t73 = sin(qJ(3));
t134 = t43 * t73;
t74 = cos(qJ(3));
t133 = t74 * t43;
t40 = t49 * t74;
t108 = t126 * t71;
t50 = t125 * t72 + t108;
t44 = t50 * qJD(5);
t124 = t44 * t73;
t16 = t40 * qJD(3) + t124;
t121 = t74 * t44;
t68 = t73 * qJD(3);
t92 = t125 * t68;
t93 = t126 * t68;
t17 = -t71 * t92 + t72 * t93 + t121;
t39 = t49 * t73;
t132 = (-t39 * t74 + t40 * t73) * qJD(3) - t16 * t73 - t74 * t17;
t75 = -pkin(1) - pkin(7);
t104 = -t71 * t75 + pkin(4);
t128 = t72 * pkin(8);
t41 = -t74 * qJD(4) + qJD(2) + (pkin(3) * t74 + qJ(4) * t73) * qJD(3);
t36 = t72 * t41;
t12 = t36 + (t104 * t74 + t73 * t128) * qJD(3);
t112 = t71 * t68;
t116 = t74 * qJD(3);
t109 = t75 * t116;
t28 = t72 * t109 + t71 * t41;
t20 = pkin(8) * t112 + t28;
t127 = t73 * pkin(3);
t87 = -t74 * qJ(4) + t127;
t51 = qJ(2) + t87;
t47 = t72 * t51;
t25 = t104 * t73 - t74 * t128 + t47;
t123 = t71 * t74;
t122 = t73 * t75;
t61 = t72 * t122;
t33 = t71 * t51 + t61;
t29 = -pkin(8) * t123 + t33;
t78 = t125 * t25 + t126 * t29;
t4 = -t78 * qJD(5) + t126 * t12 - t125 * t20;
t131 = -0.2e1 * t43;
t130 = 0.2e1 * qJD(2);
t129 = 2 * qJD(6);
t120 = pkin(8) + qJ(4);
t119 = t71 ^ 2 + t72 ^ 2;
t118 = qJD(4) * t73;
t117 = t73 * qJD(6);
t115 = qJ(2) * qJD(3);
t114 = t71 * t122;
t113 = pkin(5) * t116;
t111 = t72 * t68;
t65 = t75 * t68;
t110 = t73 * t116;
t67 = -t72 * pkin(4) - pkin(3);
t105 = qJ(6) * t116;
t103 = t119 * t74;
t48 = pkin(4) * t123 - t74 * t75;
t102 = qJD(4) * t125;
t99 = t126 * qJD(4);
t98 = t119 * qJD(4);
t97 = 0.2e1 * t110;
t94 = 0.2e1 * t98;
t89 = t120 * t125;
t27 = -t71 * t109 + t36;
t85 = -t27 * t71 + t28 * t72;
t32 = t47 - t114;
t84 = -t32 * t71 + t33 * t72;
t42 = -pkin(4) * t112 + t65;
t82 = t120 * t108;
t58 = t120 * t72;
t13 = qJD(5) * t82 + t58 * t100 + t71 * t102 - t72 * t99;
t31 = t126 * t58 - t71 * t89;
t81 = -t31 * t116 + t13 * t73;
t14 = t58 * t101 + t72 * t102 + (-qJD(5) * t89 + t99) * t71;
t30 = t125 * t58 + t82;
t80 = -t30 * t116 - t14 * t73;
t79 = t50 * t68 + t133;
t77 = -t125 * t29 + t126 * t25;
t3 = t29 * t100 - t25 * t101 - t125 * t12 - t126 * t20;
t38 = t50 * t74;
t18 = qJD(3) * t38 - t134;
t19 = -t71 * t93 - t72 * t92 - t133;
t37 = t50 * t73;
t76 = -t18 * t73 + t38 * t68 + (-qJD(3) * t37 - t19) * t74;
t24 = t49 * pkin(5) - t50 * qJ(6) + t67;
t23 = t49 * t68 - t121;
t9 = t44 * pkin(5) + t43 * qJ(6) - t50 * qJD(6);
t8 = t38 * pkin(5) + t40 * qJ(6) + t48;
t7 = -t73 * pkin(5) - t77;
t6 = t73 * qJ(6) + t78;
t5 = t19 * pkin(5) + t17 * qJ(6) + t40 * qJD(6) + t42;
t2 = -t113 - t4;
t1 = t105 - t3 + t117;
t10 = [0, 0, 0, 0, t130, qJ(2) * t130, -0.2e1 * t110, 0.2e1 * (t73 ^ 2 - t74 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t73 + 0.2e1 * t74 * t115, 0.2e1 * qJD(2) * t74 - 0.2e1 * t73 * t115, 0.2e1 * t27 * t73 + 0.2e1 * (t32 + 0.2e1 * t114) * t116, -0.2e1 * t28 * t73 + 0.2e1 * (-t33 + 0.2e1 * t61) * t116, 0.2e1 * (-t27 * t72 - t28 * t71) * t74 + 0.2e1 * (t32 * t72 + t33 * t71) * t68, -0.2e1 * t75 ^ 2 * t110 + 0.2e1 * t32 * t27 + 0.2e1 * t33 * t28, 0.2e1 * t40 * t17, 0.2e1 * t17 * t38 + 0.2e1 * t40 * t19, -0.2e1 * t40 * t116 - 0.2e1 * t17 * t73, -0.2e1 * t38 * t116 - 0.2e1 * t19 * t73, t97, 0.2e1 * t77 * t116 + 0.2e1 * t48 * t19 + 0.2e1 * t42 * t38 + 0.2e1 * t4 * t73, -0.2e1 * t78 * t116 - 0.2e1 * t48 * t17 + 0.2e1 * t3 * t73 - 0.2e1 * t42 * t40, -0.2e1 * t7 * t116 + 0.2e1 * t8 * t19 - 0.2e1 * t2 * t73 + 0.2e1 * t5 * t38, -0.2e1 * t1 * t38 - 0.2e1 * t7 * t17 - 0.2e1 * t6 * t19 - 0.2e1 * t2 * t40, 0.2e1 * t1 * t73 + 0.2e1 * t6 * t116 + 0.2e1 * t8 * t17 + 0.2e1 * t5 * t40, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t73 + (t84 - 0.2e1 * t122) * t116, 0, 0, 0, 0, 0, t76, -t132, t76, t16 * t38 - t37 * t17 - t18 * t40 + t39 * t19, t132, -t1 * t39 - t6 * t16 + t7 * t18 + t2 * t37 - t5 * t74 + t8 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t119) * t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t39 * t16 + 0.2e1 * t37 * t18 - 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t116, 0, -t65, -t109, -t71 * t118 + (t87 * t71 - t61) * qJD(3), -t72 * t118 + (t87 * t72 + t114) * qJD(3), t85, -pkin(3) * t65 + t85 * qJ(4) + t84 * qJD(4), -t17 * t50 + t40 * t43, t17 * t49 - t50 * t19 + t43 * t38 + t40 * t44, t50 * t116 - t134, -t49 * t116 - t124, 0, t67 * t19 + t42 * t49 + t48 * t44 + t80, -t67 * t17 + t42 * t50 - t48 * t43 + t81, t24 * t19 + t9 * t38 + t8 * t44 + t5 * t49 + t80, -t1 * t49 + t13 * t38 - t14 * t40 - t30 * t17 - t31 * t19 + t2 * t50 - t7 * t43 - t6 * t44, t24 * t17 + t9 * t40 + t8 * t43 - t5 * t50 - t81, t1 * t31 - t6 * t13 + t7 * t14 + t2 * t30 + t5 * t24 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t116, -t111, t112, qJD(3) * t103, t73 * t98 + (qJ(4) * t103 - t127) * qJD(3), 0, 0, 0, 0, 0, t23, t79, t23, t16 * t49 + t18 * t50 - t37 * t43 + t39 * t44, -t79, t39 * t13 + t37 * t14 - t16 * t31 + t18 * t30 + t24 * t68 - t74 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, qJ(4) * t94, t50 * t131, 0.2e1 * t43 * t49 - 0.2e1 * t50 * t44, 0, 0, 0, 0.2e1 * t67 * t44, t67 * t131, 0.2e1 * t24 * t44 + 0.2e1 * t9 * t49, 0.2e1 * t13 * t49 + 0.2e1 * t14 * t50 - 0.2e1 * t30 * t43 - 0.2e1 * t31 * t44, 0.2e1 * t24 * t43 - 0.2e1 * t9 * t50, -0.2e1 * t31 * t13 + 0.2e1 * t30 * t14 + 0.2e1 * t24 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111, 0, t65, 0, 0, 0, 0, 0, t19, -t17, t19, 0, t17, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, t44, 0, t43, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, t116, t4, t3, t4 + 0.2e1 * t113, pkin(5) * t17 - t19 * qJ(6) - t38 * qJD(6), 0.2e1 * t105 - t3 + 0.2e1 * t117, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t16, -t18, 0, -t16, -t18 * pkin(5) - t16 * qJ(6) - t39 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, -t14, t13, -t14, pkin(5) * t43 - t44 * qJ(6) - t49 * qJD(6), -t13, -t14 * pkin(5) - t13 * qJ(6) + t31 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJ(6) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
