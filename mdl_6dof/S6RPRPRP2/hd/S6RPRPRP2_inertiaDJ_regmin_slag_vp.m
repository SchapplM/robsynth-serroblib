% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP2
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
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:13
% EndTime: 2019-03-09 03:06:15
% DurationCPUTime: 0.98s
% Computational Cost: add. (1686->151), mult. (3509->269), div. (0->0), fcn. (3217->8), ass. (0->95)
t58 = cos(qJ(5));
t56 = sin(qJ(5));
t71 = t58 * pkin(5) + t56 * qJ(6);
t114 = t71 * qJD(5) - t58 * qJD(6);
t54 = sin(pkin(10));
t59 = cos(qJ(3));
t93 = cos(pkin(10));
t77 = t93 * t59;
t57 = sin(qJ(3));
t89 = t57 * qJD(3);
t37 = qJD(3) * t77 - t54 * t89;
t78 = t93 * t57;
t42 = t54 * t59 + t78;
t70 = pkin(5) * t56 - qJ(6) * t58;
t117 = -t114 * t42 - t70 * t37;
t103 = t54 * t57;
t41 = -t77 + t103;
t48 = -cos(pkin(9)) * pkin(1) - pkin(2);
t65 = -t59 * pkin(3) + t48;
t18 = t41 * pkin(4) - t42 * pkin(8) + t65;
t46 = sin(pkin(9)) * pkin(1) + pkin(7);
t94 = qJ(4) + t46;
t39 = t94 * t59;
t25 = -t94 * t103 + t93 * t39;
t115 = t56 * t18 + t58 * t25;
t75 = qJD(3) * t94;
t30 = t59 * qJD(4) - t57 * t75;
t63 = -t57 * qJD(4) - t59 * t75;
t11 = t93 * t30 + t54 * t63;
t36 = t42 * qJD(3);
t50 = pkin(3) * t89;
t16 = t36 * pkin(4) - t37 * pkin(8) + t50;
t4 = -qJD(5) * t115 - t56 * t11 + t58 * t16;
t51 = qJD(5) * t58;
t92 = qJD(5) * t56;
t3 = -t58 * t11 - t56 * t16 - t18 * t51 + t25 * t92;
t91 = t41 * qJD(6);
t95 = t36 * qJ(6);
t1 = -t3 + t91 + t95;
t110 = t36 * pkin(5);
t2 = -t110 - t4;
t6 = t41 * qJ(6) + t115;
t69 = t58 * t18 - t56 * t25;
t7 = -t41 * pkin(5) - t69;
t73 = t56 * t7 + t58 * t6;
t113 = t73 * qJD(5) + t1 * t56 - t2 * t58;
t112 = 0.2e1 * qJD(5);
t111 = 0.2e1 * qJD(6);
t10 = t54 * t30 - t93 * t63;
t109 = t10 * t56;
t45 = t54 * pkin(3) + pkin(8);
t108 = t36 * t45;
t107 = t41 * t36;
t106 = t41 * t45;
t105 = t42 * t58;
t52 = t56 ^ 2;
t104 = t52 * t37;
t53 = t58 ^ 2;
t31 = t53 * t37;
t102 = t56 * t36;
t101 = t56 * t37;
t100 = t58 * t36;
t99 = t58 * t37;
t97 = t42 * t100 + t41 * t99;
t96 = t52 - t53;
t90 = t56 * qJD(6);
t87 = t59 * qJD(3);
t47 = -t93 * pkin(3) - pkin(4);
t86 = t47 * t112;
t85 = 0.2e1 * t87;
t84 = t42 * t92;
t83 = t45 * t92;
t82 = t45 * t51;
t81 = t56 * t51;
t80 = -0.4e1 * t56 * t105;
t24 = t54 * t39 + t94 * t78;
t76 = t96 * qJD(5);
t72 = -t56 * t6 + t58 * t7;
t67 = t37 * t47 - t108;
t66 = -t42 * t47 + t106;
t22 = t41 * t51 + t102;
t20 = t41 * t92 - t100;
t23 = t42 * t51 + t101;
t21 = t84 - t99;
t38 = -t71 + t47;
t5 = t10 - t117;
t64 = -t5 + (t38 * t42 - t106) * qJD(5);
t35 = -pkin(5) * t92 + qJ(6) * t51 + t90;
t8 = t70 * t42 + t24;
t62 = qJD(5) * t8 - t35 * t42 + t37 * t38 - t108;
t60 = t72 * qJD(5) + t1 * t58 + t2 * t56;
t40 = t42 ^ 2;
t26 = t42 * t31;
t19 = -t31 - t104;
t9 = [0, 0, 0, 0, t57 * t85, 0.2e1 * (-t57 ^ 2 + t59 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t48 * t89, t48 * t85, 0.2e1 * t10 * t42 - 0.2e1 * t11 * t41 + 0.2e1 * t24 * t37 - 0.2e1 * t25 * t36, 0.2e1 * t24 * t10 + 0.2e1 * t25 * t11 + 0.2e1 * t65 * t50, -0.2e1 * t40 * t81 + 0.2e1 * t26, t96 * t40 * t112 + t37 * t80, -0.2e1 * t41 * t84 + 0.2e1 * t97, -0.2e1 * t42 * t102 - 0.2e1 * t23 * t41, 0.2e1 * t107, 0.2e1 * t42 * t109 + 0.2e1 * t23 * t24 + 0.2e1 * t69 * t36 + 0.2e1 * t4 * t41, 0.2e1 * t10 * t105 - 0.2e1 * t115 * t36 - 0.2e1 * t21 * t24 + 0.2e1 * t3 * t41, 0.2e1 * t8 * t101 - 0.2e1 * t2 * t41 - 0.2e1 * t7 * t36 + 0.2e1 * (t5 * t56 + t8 * t51) * t42, -0.2e1 * t113 * t42 + 0.2e1 * t72 * t37, -0.2e1 * t8 * t99 + 0.2e1 * t1 * t41 + 0.2e1 * t6 * t36 + 0.2e1 * (-t5 * t58 + t8 * t92) * t42, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t41 + t11 * t42 + t24 * t36 + t25 * t37, 0, 0, 0, 0, 0, 0 (-t42 * t36 - t37 * t41) * t58 + t97, 0, 0, 0, t8 * t36 + t73 * t37 + t5 * t41 + t60 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 * t37 + 0.2e1 * t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 * t104 + 0.2e1 * t107 + 0.2e1 * t26; 0, 0, 0, 0, 0, 0, t87, -t89, 0, -t46 * t87, t46 * t89 (-t36 * t54 - t93 * t37) * pkin(3) (-t93 * t10 + t11 * t54) * pkin(3), -t42 * t76 + t56 * t99, qJD(5) * t80 - t104 + t31, t22, -t20, 0, -t10 * t58 + t67 * t56 + (t24 * t56 - t66 * t58) * qJD(5), t109 + t67 * t58 + (t24 * t58 + t66 * t56) * qJD(5), t62 * t56 + t64 * t58, t60, t64 * t56 - t62 * t58, -t8 * t35 + t5 * t38 + t60 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t87, 0 (-t93 * t36 + t37 * t54) * pkin(3), 0, 0, 0, 0, 0, t20, t22, t20, -t19, -t22, -t41 * t35 + t36 * t38 + (t52 + t53) * t45 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t81, -0.2e1 * t76, 0, 0, 0, t56 * t86, t58 * t86, 0.2e1 * t35 * t58 + 0.2e1 * t38 * t92, 0, 0.2e1 * t35 * t56 - 0.2e1 * t38 * t51, -0.2e1 * t38 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, -t20, -t22, -t20, t19, t22, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, t36, t4, t3, t4 + 0.2e1 * t110, -t71 * t37 + (t70 * qJD(5) - t90) * t42, -t3 + 0.2e1 * t91 + 0.2e1 * t95, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, -t23, 0, -t21, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t92, 0, -t82, t83, -t82, -t114, -t83, -t114 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t51, -t92, 0, t51, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, qJ(6) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
