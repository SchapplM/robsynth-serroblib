% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:43
% EndTime: 2019-03-09 01:58:45
% DurationCPUTime: 0.65s
% Computational Cost: add. (1193->131), mult. (2556->224), div. (0->0), fcn. (2413->8), ass. (0->86)
t56 = cos(pkin(10));
t101 = cos(qJ(4));
t74 = qJD(4) * t101;
t55 = sin(pkin(10));
t59 = sin(qJ(4));
t95 = t59 * t55;
t36 = qJD(4) * t95 - t56 * t74;
t58 = sin(qJ(5));
t42 = t101 * t55 + t56 * t59;
t60 = cos(qJ(5));
t86 = qJD(5) * t60;
t79 = t42 * t86;
t62 = -t36 * t58 + t79;
t105 = t62 * pkin(5);
t48 = sin(pkin(9)) * pkin(1) + qJ(3);
t102 = pkin(7) + t48;
t38 = t102 * t55;
t39 = t102 * t56;
t25 = t101 * t39 - t38 * t59;
t15 = t60 * t25;
t77 = t101 * t56;
t41 = -t77 + t95;
t43 = -cos(pkin(9)) * pkin(1) - t56 * pkin(3) - pkin(2);
t16 = pkin(4) * t41 - pkin(8) * t42 + t43;
t91 = t16 * t58 + t15;
t53 = t58 ^ 2;
t54 = t60 ^ 2;
t72 = (t53 - t54) * qJD(5);
t37 = t42 * qJD(4);
t63 = qJ(6) * t36 - qJD(6) * t42;
t26 = pkin(4) * t37 + pkin(8) * t36;
t9 = t38 * t74 - qJD(3) * t77 + (qJD(3) * t55 + qJD(4) * t39) * t59;
t76 = t26 * t60 + t58 * t9;
t88 = qJ(6) * t42;
t1 = t37 * pkin(5) + t63 * t60 + (-t15 + (-t16 + t88) * t58) * qJD(5) + t76;
t84 = t16 * t86 + t26 * t58 - t60 * t9;
t2 = -qJ(6) * t79 + (-qJD(5) * t25 + t63) * t58 + t84;
t75 = t16 * t60 - t25 * t58;
t5 = pkin(5) * t41 - t60 * t88 + t75;
t6 = -t58 * t88 + t91;
t66 = t5 * t58 - t6 * t60;
t104 = qJD(5) * t66 - t1 * t60 - t2 * t58;
t103 = -0.2e1 * t36;
t10 = qJD(3) * t42 + qJD(4) * t25;
t100 = t10 * t60;
t99 = t41 * t37;
t98 = t42 * t58;
t97 = t53 * t36;
t31 = t54 * t36;
t96 = t58 * t37;
t94 = t60 * t36;
t93 = t60 * t37;
t92 = -qJ(6) - pkin(8);
t90 = -t41 * t94 + t42 * t93;
t87 = qJD(5) * t58;
t85 = -0.2e1 * pkin(4) * qJD(5);
t83 = t58 * t94;
t82 = pkin(5) * t87;
t81 = t41 * t87;
t80 = t42 * t87;
t78 = t58 * t86;
t73 = qJD(5) * t92;
t71 = 0.2e1 * (t55 ^ 2 + t56 ^ 2) * qJD(3);
t70 = pkin(4) * t36 - pkin(8) * t37;
t69 = pkin(4) * t42 + pkin(8) * t41;
t67 = t5 * t60 + t58 * t6;
t24 = t101 * t38 + t39 * t59;
t65 = t36 * t41 - t37 * t42;
t44 = t92 * t58;
t45 = t92 * t60;
t64 = t44 * t58 + t45 * t60;
t21 = t80 + t94;
t22 = t41 * t86 + t96;
t34 = t60 * qJD(6) + t58 * t73;
t35 = -t58 * qJD(6) + t60 * t73;
t61 = t34 * t60 - t35 * t58 + (-t44 * t60 + t45 * t58) * qJD(5);
t50 = -pkin(5) * t60 - pkin(4);
t40 = t42 ^ 2;
t27 = t42 * t31;
t20 = t81 - t93;
t17 = t31 + t97;
t11 = pkin(5) * t98 + t24;
t7 = t10 + t105;
t4 = -qJD(5) * t91 + t76;
t3 = t25 * t87 - t84;
t8 = [0, 0, 0, 0, 0, 0, t71, t48 * t71, t42 * t103, 0.2e1 * t65, 0, 0, 0, 0.2e1 * t43 * t37, t43 * t103, -0.2e1 * t40 * t78 - 0.2e1 * t27, 0.2e1 * t40 * t72 + 0.4e1 * t42 * t83, -0.2e1 * t41 * t80 + 0.2e1 * t90, -0.2e1 * t41 * t62 - 0.2e1 * t42 * t96, 0.2e1 * t99, 0.2e1 * t10 * t98 + 0.2e1 * t24 * t62 + 0.2e1 * t37 * t75 + 0.2e1 * t4 * t41, 0.2e1 * t100 * t42 - 0.2e1 * t21 * t24 + 0.2e1 * t3 * t41 - 0.2e1 * t37 * t91, 0.2e1 * t104 * t42 + 0.2e1 * t36 * t67, 0.2e1 * t1 * t5 + 0.2e1 * t11 * t7 + 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t65 + t90, 0, t11 * t37 + t7 * t41 + t66 * t36 + (-qJD(5) * t67 - t1 * t58 + t2 * t60) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t42 * t97 - 0.2e1 * t27 + 0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, 0, 0, 0, 0, -t20, -t22, t17, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, -t10, t9, -t42 * t72 - t83, -0.4e1 * t42 * t78 - t31 + t97, t22, -t20, 0, -t100 + t70 * t58 + (t24 * t58 - t60 * t69) * qJD(5), t10 * t58 + t70 * t60 + (t24 * t60 + t58 * t69) * qJD(5) (-t35 * t42 + t36 * t44 + t2 + (t42 * t45 - t5) * qJD(5)) * t60 + (-t34 * t42 - t36 * t45 - t1 + (t42 * t44 - t6) * qJD(5)) * t58, t1 * t44 + t11 * t82 - t2 * t45 + t34 * t6 + t35 * t5 + t50 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, 0, 0, 0, 0, t20, t22, -t17, pkin(5) * t81 + t36 * t64 + t37 * t50 + t42 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t64 + t58 * t34 + t60 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t78, -0.2e1 * t72, 0, 0, 0, t58 * t85, t60 * t85, 0.2e1 * t61, -0.2e1 * t34 * t45 + 0.2e1 * t35 * t44 + 0.2e1 * t50 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t62, t37, t4, t3, t21 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t21, 0, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t87, 0, -pkin(8) * t86, pkin(8) * t87, -pkin(5) * t86, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
