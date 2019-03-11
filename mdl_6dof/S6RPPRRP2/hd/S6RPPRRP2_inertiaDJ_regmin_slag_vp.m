% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP2
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
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:17
% EndTime: 2019-03-09 02:01:19
% DurationCPUTime: 0.87s
% Computational Cost: add. (1527->137), mult. (3227->232), div. (0->0), fcn. (3065->8), ass. (0->88)
t59 = cos(qJ(5));
t57 = sin(qJ(5));
t68 = t59 * pkin(5) + t57 * qJ(6);
t111 = t68 * qJD(5) - t59 * qJD(6);
t55 = cos(pkin(10));
t102 = cos(qJ(4));
t76 = qJD(4) * t102;
t54 = sin(pkin(10));
t58 = sin(qJ(4));
t95 = t58 * t54;
t35 = qJD(4) * t95 - t55 * t76;
t41 = t102 * t54 + t58 * t55;
t67 = pkin(5) * t57 - qJ(6) * t59;
t114 = -t111 * t41 + t67 * t35;
t78 = t102 * t55;
t40 = -t78 + t95;
t42 = -cos(pkin(9)) * pkin(1) - t55 * pkin(3) - pkin(2);
t15 = t40 * pkin(4) - t41 * pkin(8) + t42;
t46 = sin(pkin(9)) * pkin(1) + qJ(3);
t103 = pkin(7) + t46;
t37 = t103 * t54;
t38 = t103 * t55;
t24 = t102 * t38 - t58 * t37;
t112 = t57 * t15 + t59 * t24;
t52 = t57 ^ 2;
t53 = t59 ^ 2;
t75 = (t52 - t53) * qJD(5);
t10 = t37 * t76 - qJD(3) * t78 + (qJD(3) * t54 + qJD(4) * t38) * t58;
t105 = t35 * pkin(8);
t36 = t41 * qJD(4);
t25 = t36 * pkin(4) + t105;
t4 = -qJD(5) * t112 + t57 * t10 + t59 * t25;
t49 = qJD(5) * t59;
t88 = qJD(5) * t57;
t3 = t59 * t10 - t15 * t49 + t24 * t88 - t57 * t25;
t87 = t40 * qJD(6);
t89 = t36 * qJ(6);
t1 = -t3 + t87 + t89;
t104 = t36 * pkin(5);
t2 = -t104 - t4;
t6 = t40 * qJ(6) + t112;
t66 = t59 * t15 - t57 * t24;
t7 = -t40 * pkin(5) - t66;
t70 = t57 * t7 + t59 * t6;
t110 = t70 * qJD(5) + t1 * t57 - t2 * t59;
t109 = -0.2e1 * t35;
t108 = 0.2e1 * qJD(6);
t107 = pkin(8) * t36;
t106 = pkin(8) * t40;
t11 = t41 * qJD(3) + t24 * qJD(4);
t101 = t11 * t57;
t100 = t11 * t59;
t99 = t40 * t36;
t98 = t52 * t35;
t30 = t53 * t35;
t97 = t57 * t35;
t96 = t57 * t36;
t94 = t59 * t35;
t93 = t59 * t36;
t91 = -t40 * t94 + t41 * t93;
t86 = t57 * qJD(6);
t84 = -0.2e1 * pkin(4) * qJD(5);
t83 = t57 * t94;
t82 = pkin(8) * t88;
t81 = pkin(8) * t49;
t80 = t41 * t88;
t79 = t57 * t49;
t74 = 0.2e1 * (t54 ^ 2 + t55 ^ 2) * qJD(3);
t73 = pkin(4) * t35 - t107;
t72 = pkin(4) * t41 + t106;
t69 = t57 * t6 - t59 * t7;
t23 = t102 * t37 + t58 * t38;
t64 = t35 * t40 - t41 * t36;
t22 = t41 * t49 - t97;
t20 = t80 + t94;
t21 = t40 * t49 + t96;
t19 = t40 * t88 - t93;
t43 = -pkin(4) - t68;
t5 = t11 - t114;
t63 = -t5 + (t41 * t43 - t106) * qJD(5);
t34 = -pkin(5) * t88 + qJ(6) * t49 + t86;
t8 = t67 * t41 + t23;
t62 = -qJD(5) * t8 + t34 * t41 + t35 * t43 + t107;
t60 = -t69 * qJD(5) + t1 * t59 + t2 * t57;
t39 = t41 ^ 2;
t26 = t41 * t30;
t16 = t30 + t98;
t9 = [0, 0, 0, 0, 0, 0, t74, t46 * t74, t41 * t109, 0.2e1 * t64, 0, 0, 0, 0.2e1 * t42 * t36, t42 * t109, -0.2e1 * t39 * t79 - 0.2e1 * t26, 0.2e1 * t39 * t75 + 0.4e1 * t41 * t83, -0.2e1 * t40 * t80 + 0.2e1 * t91, -0.2e1 * t22 * t40 - 0.2e1 * t41 * t96, 0.2e1 * t99, 0.2e1 * t41 * t101 + 0.2e1 * t22 * t23 + 0.2e1 * t66 * t36 + 0.2e1 * t4 * t40, 0.2e1 * t41 * t100 - 0.2e1 * t112 * t36 - 0.2e1 * t20 * t23 + 0.2e1 * t3 * t40, -0.2e1 * t8 * t97 - 0.2e1 * t2 * t40 - 0.2e1 * t7 * t36 + 0.2e1 * (t8 * t49 + t5 * t57) * t41, -0.2e1 * t110 * t41 + 0.2e1 * t69 * t35, 0.2e1 * t8 * t94 + 0.2e1 * t1 * t40 + 0.2e1 * t6 * t36 + 0.2e1 * (-t5 * t59 + t8 * t88) * t41, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t59 + t91, 0, 0, 0, -t70 * t35 + t8 * t36 + t5 * t40 + t60 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t41 * t98 - 0.2e1 * t26 + 0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, 0, 0, 0, 0, -t19, -t21, -t19, t16, t21, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t11, t10, -t41 * t75 - t83, -0.4e1 * t41 * t79 - t30 + t98, t21, -t19, 0, -t100 + t73 * t57 + (t23 * t57 - t72 * t59) * qJD(5), t101 + t73 * t59 + (t23 * t59 + t72 * t57) * qJD(5), -t62 * t57 + t63 * t59, t60, t63 * t57 + t62 * t59, t60 * pkin(8) - t8 * t34 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35, 0, 0, 0, 0, 0, t19, t21, t19, -t16, -t21, -t40 * t34 + t36 * t43 + (-t52 - t53) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79, -0.2e1 * t75, 0, 0, 0, t57 * t84, t59 * t84, 0.2e1 * t34 * t59 + 0.2e1 * t43 * t88, 0, 0.2e1 * t34 * t57 - 0.2e1 * t43 * t49, -0.2e1 * t43 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22, t36, t4, t3, t4 + 0.2e1 * t104, t68 * t35 + (t67 * qJD(5) - t86) * t41, -t3 + 0.2e1 * t87 + 0.2e1 * t89, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t20, -t22, 0, -t20, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t49, -t88, 0, t49, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t88, 0, -t81, t82, -t81, -t111, -t82, -t111 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, qJ(6) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t20, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
