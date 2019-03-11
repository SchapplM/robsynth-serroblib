% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:03
% EndTime: 2019-03-09 03:03:05
% DurationCPUTime: 0.72s
% Computational Cost: add. (1353->145), mult. (2833->260), div. (0->0), fcn. (2584->8), ass. (0->94)
t55 = sin(pkin(10));
t60 = cos(qJ(3));
t92 = cos(pkin(10));
t75 = t92 * t60;
t58 = sin(qJ(3));
t89 = t58 * qJD(3);
t38 = qJD(3) * t75 - t55 * t89;
t57 = sin(qJ(5));
t76 = t92 * t58;
t44 = t55 * t60 + t76;
t59 = cos(qJ(5));
t90 = qJD(5) * t59;
t81 = t44 * t90;
t63 = t57 * t38 + t81;
t109 = t63 * pkin(5);
t102 = t55 * t58;
t49 = sin(pkin(9)) * pkin(1) + pkin(7);
t94 = qJ(4) + t49;
t41 = t94 * t60;
t26 = -t94 * t102 + t92 * t41;
t18 = t59 * t26;
t43 = -t75 + t102;
t51 = -cos(pkin(9)) * pkin(1) - pkin(2);
t64 = -t60 * pkin(3) + t51;
t19 = t43 * pkin(4) - t44 * pkin(8) + t64;
t98 = t57 * t19 + t18;
t37 = t44 * qJD(3);
t65 = -qJ(6) * t38 - qJD(6) * t44;
t73 = qJD(3) * t94;
t31 = t60 * qJD(4) - t58 * t73;
t62 = -t58 * qJD(4) - t60 * t73;
t11 = t92 * t31 + t55 * t62;
t52 = pkin(3) * t89;
t17 = t37 * pkin(4) - t38 * pkin(8) + t52;
t78 = -t57 * t11 + t59 * t17;
t95 = qJ(6) * t44;
t1 = t37 * pkin(5) + t65 * t59 + (-t18 + (-t19 + t95) * t57) * qJD(5) + t78;
t87 = t59 * t11 + t57 * t17 + t19 * t90;
t2 = -qJ(6) * t81 + (-qJD(5) * t26 + t65) * t57 + t87;
t77 = t59 * t19 - t57 * t26;
t5 = t43 * pkin(5) - t59 * t95 + t77;
t6 = -t57 * t95 + t98;
t69 = t5 * t57 - t59 * t6;
t108 = t69 * qJD(5) - t1 * t59 - t2 * t57;
t107 = 0.2e1 * qJD(5);
t106 = t43 * t37;
t105 = t44 * t57;
t104 = t44 * t59;
t53 = t57 ^ 2;
t103 = t53 * t38;
t54 = t59 ^ 2;
t34 = t54 * t38;
t101 = t57 * t37;
t100 = t59 * t37;
t99 = t59 * t38;
t97 = t44 * t100 + t43 * t99;
t96 = t53 - t54;
t48 = t55 * pkin(3) + pkin(8);
t93 = qJ(6) + t48;
t91 = qJD(5) * t57;
t88 = t60 * qJD(3);
t50 = -t92 * pkin(3) - pkin(4);
t86 = t50 * t107;
t85 = 0.2e1 * t88;
t84 = pkin(5) * t91;
t83 = t43 * t91;
t82 = t44 * t91;
t80 = t57 * t90;
t79 = -0.4e1 * t57 * t104;
t10 = t55 * t31 - t92 * t62;
t25 = t55 * t41 + t94 * t76;
t74 = t96 * qJD(5);
t72 = qJD(5) * t93;
t70 = -t5 * t59 - t57 * t6;
t68 = -t37 * t48 + t38 * t50;
t39 = t93 * t57;
t40 = t93 * t59;
t67 = t39 * t57 + t40 * t59;
t66 = t43 * t48 - t44 * t50;
t23 = t43 * t90 + t101;
t22 = t82 - t99;
t32 = t59 * qJD(6) - t57 * t72;
t33 = -t57 * qJD(6) - t59 * t72;
t61 = t32 * t59 - t33 * t57 + (t39 * t59 - t40 * t57) * qJD(5);
t45 = -t59 * pkin(5) + t50;
t42 = t44 ^ 2;
t27 = t44 * t34;
t21 = t83 - t100;
t20 = -t34 - t103;
t9 = pkin(5) * t105 + t25;
t7 = t10 + t109;
t4 = -t98 * qJD(5) + t78;
t3 = t26 * t91 - t87;
t8 = [0, 0, 0, 0, t58 * t85, 0.2e1 * (-t58 ^ 2 + t60 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t51 * t89, t51 * t85, 0.2e1 * t10 * t44 - 0.2e1 * t11 * t43 + 0.2e1 * t25 * t38 - 0.2e1 * t26 * t37, 0.2e1 * t25 * t10 + 0.2e1 * t26 * t11 + 0.2e1 * t64 * t52, -0.2e1 * t42 * t80 + 0.2e1 * t27, t96 * t42 * t107 + t38 * t79, -0.2e1 * t43 * t82 + 0.2e1 * t97, -0.2e1 * t44 * t101 - 0.2e1 * t63 * t43, 0.2e1 * t106, 0.2e1 * t10 * t105 + 0.2e1 * t63 * t25 + 0.2e1 * t77 * t37 + 0.2e1 * t4 * t43, 0.2e1 * t10 * t104 - 0.2e1 * t22 * t25 + 0.2e1 * t3 * t43 - 0.2e1 * t98 * t37, 0.2e1 * t108 * t44 + 0.2e1 * t70 * t38, 0.2e1 * t5 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * t9 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t43 + t11 * t44 + t25 * t37 + t26 * t38, 0, 0, 0, 0, 0, 0 (-t44 * t37 - t38 * t43) * t59 + t97, 0, t9 * t37 + t7 * t43 - t69 * t38 + (t70 * qJD(5) - t1 * t57 + t2 * t59) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t44 * t38 + 0.2e1 * t106, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t44 * t103 + 0.2e1 * t106 + 0.2e1 * t27; 0, 0, 0, 0, 0, 0, t88, -t89, 0, -t49 * t88, t49 * t89 (-t37 * t55 - t92 * t38) * pkin(3) (-t92 * t10 + t11 * t55) * pkin(3), -t44 * t74 + t57 * t99, qJD(5) * t79 - t103 + t34, t23, -t21, 0, -t10 * t59 + t68 * t57 + (t25 * t57 - t66 * t59) * qJD(5), t10 * t57 + t68 * t59 + (t25 * t59 + t66 * t57) * qJD(5) (-t33 * t44 + t38 * t39 + t2 + (-t40 * t44 - t5) * qJD(5)) * t59 + (-t32 * t44 - t38 * t40 - t1 + (-t39 * t44 - t6) * qJD(5)) * t57, -t1 * t39 + t2 * t40 + t6 * t32 + t5 * t33 + t7 * t45 + t9 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0 (-t92 * t37 + t38 * t55) * pkin(3), 0, 0, 0, 0, 0, t21, t23, -t20, pkin(5) * t83 + t37 * t45 + t67 * t38 + t61 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80, -0.2e1 * t74, 0, 0, 0, t57 * t86, t59 * t86, 0.2e1 * t61, 0.2e1 * t40 * t32 - 0.2e1 * t39 * t33 + 0.2e1 * t45 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, -t21, -t23, t20, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * qJD(5) + t32 * t57 + t33 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t63, t37, t4, t3, t22 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t22, 0, -t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t91, 0, -t48 * t90, t48 * t91, -pkin(5) * t90, t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90, 0, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
