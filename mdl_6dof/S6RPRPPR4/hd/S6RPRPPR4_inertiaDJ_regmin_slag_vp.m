% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:23
% EndTime: 2019-03-09 02:48:26
% DurationCPUTime: 1.04s
% Computational Cost: add. (1455->150), mult. (3490->286), div. (0->0), fcn. (3378->8), ass. (0->97)
t76 = sin(pkin(10));
t72 = t76 ^ 2;
t78 = cos(pkin(10));
t112 = t78 ^ 2 + t72;
t127 = t112 * qJD(4);
t122 = cos(qJ(3));
t115 = pkin(7) + qJ(2);
t77 = sin(pkin(9));
t63 = t115 * t77;
t79 = cos(pkin(9));
t65 = t115 * t79;
t81 = sin(qJ(3));
t38 = t122 * t65 - t81 * t63;
t60 = t122 * t77 + t81 * t79;
t29 = t60 * qJD(2) + t38 * qJD(3);
t126 = t78 * t60 * qJD(5) - t29;
t82 = cos(qJ(6));
t109 = qJD(6) * t82;
t80 = sin(qJ(6));
t110 = qJD(6) * t80;
t48 = t76 * t109 - t78 * t110;
t57 = t80 * t76 + t82 * t78;
t59 = t82 * t76 - t80 * t78;
t71 = -t79 * pkin(2) - pkin(1);
t125 = 0.2e1 * t71;
t124 = pkin(4) + pkin(5);
t123 = pkin(8) * t60;
t101 = t122 * t79;
t85 = -t81 * t77 + t101;
t49 = t85 * qJD(3);
t43 = t76 * t49;
t50 = t60 * qJD(3);
t121 = t76 * t50;
t120 = t78 * t50;
t114 = -pkin(8) + qJ(4);
t23 = t50 * pkin(3) - t49 * qJ(4) - t60 * qJD(4);
t102 = t122 * t63;
t28 = qJD(3) * t102 - qJD(2) * t101 + (qJD(2) * t77 + qJD(3) * t65) * t81;
t11 = t76 * t23 - t78 * t28;
t34 = -pkin(3) * t85 - t60 * qJ(4) + t71;
t19 = t76 * t34 + t78 * t38;
t113 = qJ(4) * t127;
t111 = qJ(5) * t78;
t108 = t72 * qJD(5);
t107 = t76 * qJD(5);
t15 = -qJ(5) * t85 + t19;
t100 = t76 * qJ(5) + pkin(3);
t24 = t76 * t28;
t10 = t78 * t23 + t24;
t35 = t76 * t38;
t18 = t78 * t34 - t35;
t5 = t50 * qJ(5) - qJD(5) * t85 + t11;
t99 = 0.2e1 * (t77 ^ 2 + t79 ^ 2) * qJD(2);
t7 = -t50 * pkin(4) - t10;
t98 = t5 * t76 - t7 * t78;
t37 = t81 * t65 + t102;
t12 = t76 * t123 + t15;
t8 = t35 + (-t34 - t123) * t78 + t124 * t85;
t97 = t82 * t12 + t80 * t8;
t96 = t80 * t12 - t82 * t8;
t94 = pkin(4) * t76 - t111;
t20 = t94 * t60 + t37;
t9 = t94 * t49 - t126;
t95 = t20 * t49 + t60 * t9;
t93 = t10 * t78 + t11 * t76;
t92 = -t10 * t76 + t11 * t78;
t91 = t29 * t60 + t37 * t49;
t47 = t57 * qJD(6);
t90 = -t47 * t85 - t59 * t50;
t62 = t114 * t76;
t64 = t114 * t78;
t89 = t82 * t62 - t80 * t64;
t88 = t80 * t62 + t82 * t64;
t87 = -qJ(4) * t50 + qJD(4) * t85;
t86 = -t124 * t76 + t111;
t84 = -pkin(3) * t49 + t87;
t61 = -t78 * pkin(4) - t100;
t83 = t49 * t61 + t87;
t56 = 0.2e1 * t127;
t51 = t124 * t78 + t100;
t44 = t78 * t49;
t32 = t57 * t60;
t31 = t59 * t60;
t30 = t112 * t49;
t27 = t59 * qJD(4) - t88 * qJD(6);
t26 = -t57 * qJD(4) - t89 * qJD(6);
t21 = -t48 * t85 + t57 * t50;
t17 = t86 * t60 - t37;
t16 = pkin(4) * t85 - t18;
t14 = t60 * t47 - t49 * t59;
t13 = -t48 * t60 - t57 * t49;
t6 = t86 * t49 + t126;
t4 = pkin(8) * t43 + t5;
t3 = -t24 + (-pkin(8) * t49 - t23) * t78 - t124 * t50;
t2 = -t97 * qJD(6) + t82 * t3 - t80 * t4;
t1 = t96 * qJD(6) - t80 * t3 - t82 * t4;
t22 = [0, 0, 0, 0, 0, t99, qJ(2) * t99, 0.2e1 * t60 * t49, 0.2e1 * t49 * t85 - 0.2e1 * t60 * t50, 0, 0, 0, t50 * t125, t49 * t125, -0.2e1 * t10 * t85 + 0.2e1 * t18 * t50 + 0.2e1 * t91 * t76, 0.2e1 * t11 * t85 - 0.2e1 * t19 * t50 + 0.2e1 * t91 * t78, -0.2e1 * t93 * t60 + 0.2e1 * (-t18 * t78 - t19 * t76) * t49, 0.2e1 * t18 * t10 + 0.2e1 * t19 * t11 + 0.2e1 * t37 * t29, -0.2e1 * t16 * t50 + 0.2e1 * t7 * t85 + 0.2e1 * t95 * t76, -0.2e1 * t98 * t60 + 0.2e1 * (-t15 * t76 + t16 * t78) * t49, 0.2e1 * t15 * t50 - 0.2e1 * t5 * t85 - 0.2e1 * t95 * t78, 0.2e1 * t15 * t5 + 0.2e1 * t16 * t7 + 0.2e1 * t20 * t9, -0.2e1 * t32 * t13, -0.2e1 * t13 * t31 - 0.2e1 * t32 * t14, -0.2e1 * t13 * t85 - 0.2e1 * t32 * t50, -0.2e1 * t14 * t85 - 0.2e1 * t31 * t50, -0.2e1 * t85 * t50, 0.2e1 * t17 * t14 + 0.2e1 * t2 * t85 - 0.2e1 * t6 * t31 + 0.2e1 * t96 * t50, 0.2e1 * t1 * t85 - 0.2e1 * t17 * t13 + 0.2e1 * t6 * t32 + 0.2e1 * t97 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, t120, -t121, -t30, t93, t120, -t30, t121, t98, 0, 0, 0, 0, 0, t21, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t50, 0, -t29, t28, -t29 * t78 + t84 * t76, t29 * t76 + t84 * t78, t92, -t29 * pkin(3) + (-t18 * t76 + t19 * t78) * qJD(4) + t92 * qJ(4), -t60 * t108 + t83 * t76 - t9 * t78, t5 * t78 + t7 * t76, -t9 * t76 + (t60 * t107 - t83) * t78, t9 * t61 + (qJ(4) * t5 + qJD(4) * t15) * t78 + (qJ(4) * t7 + qJD(4) * t16 - qJD(5) * t20) * t76, -t13 * t59 - t32 * t47, t13 * t57 - t59 * t14 - t47 * t31 - t32 * t48, t90, t21, 0, -t31 * t107 + t51 * t14 + t17 * t48 + t27 * t85 - t89 * t50 + t6 * t57, t32 * t107 - t51 * t13 - t17 * t47 + t26 * t85 + t88 * t50 + t6 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0.2e1 * t113, 0.2e1 * t78 * t107, t56, 0.2e1 * t108, -0.2e1 * t61 * t107 + 0.2e1 * t113, -0.2e1 * t59 * t47, 0.2e1 * t47 * t57 - 0.2e1 * t59 * t48, 0, 0, 0, 0.2e1 * t57 * t107 + 0.2e1 * t51 * t48, 0.2e1 * t59 * t107 - 0.2e1 * t51 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t44, 0, t29, t43, 0, -t44, t9, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, 0, 0, 0, 0, -t48, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t44, 0, t7, 0, 0, 0, 0, 0, -t110 * t85 - t82 * t50, -t109 * t85 + t80 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 * qJD(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -t50, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, 0, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t22;
