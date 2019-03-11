% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:53
% EndTime: 2019-03-08 20:16:54
% DurationCPUTime: 0.88s
% Computational Cost: add. (682->157), mult. (1729->282), div. (0->0), fcn. (1466->8), ass. (0->100)
t47 = cos(qJ(4));
t103 = t47 * pkin(9);
t44 = sin(qJ(4));
t60 = t44 * pkin(4) - t103;
t29 = qJ(3) + t60;
t43 = sin(qJ(5));
t25 = t43 * t29;
t46 = cos(qJ(5));
t49 = -pkin(2) - pkin(8);
t34 = t46 * t44 * t49;
t105 = -t25 - t34;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t97 = t37 - t39;
t66 = t97 * qJD(5);
t104 = 2 * qJD(3);
t41 = sin(pkin(6));
t45 = sin(qJ(2));
t102 = t41 * t45;
t48 = cos(qJ(2));
t101 = t41 * t48;
t100 = t43 * t49;
t99 = t47 * t49;
t98 = -qJ(6) - pkin(9);
t96 = t37 + t39;
t38 = t44 ^ 2;
t40 = t47 ^ 2;
t95 = t38 - t40;
t94 = t38 + t40;
t93 = qJ(6) * t47;
t92 = qJD(2) * t48;
t42 = cos(pkin(6));
t18 = t47 * t101 + t42 * t44;
t91 = qJD(4) * t18;
t90 = qJD(4) * t46;
t89 = qJD(5) * t43;
t88 = qJD(5) * t46;
t87 = qJD(5) * t47;
t86 = qJD(5) * t49;
t85 = t44 * qJD(4);
t84 = t46 * qJD(6);
t83 = t47 * qJD(4);
t82 = qJ(3) * qJD(4);
t81 = -0.2e1 * pkin(4) * qJD(5);
t61 = pkin(4) * t47 + pkin(9) * t44;
t27 = t61 * qJD(4) + qJD(3);
t69 = t49 * t83;
t80 = t43 * t27 + t29 * t88 + t46 * t69;
t79 = pkin(5) * t89;
t78 = t43 * t87;
t77 = t43 * t86;
t76 = t46 * t87;
t75 = t18 * t89;
t35 = qJD(2) * t102;
t74 = t41 * t92;
t73 = t43 * t88;
t72 = t49 * t85;
t71 = t46 * t85;
t70 = t44 * t83;
t68 = pkin(5) - t100;
t67 = qJD(5) * t98;
t65 = t95 * qJD(4);
t64 = 0.2e1 * t70;
t63 = t43 * t69;
t62 = t43 * t71;
t26 = t46 * t29;
t7 = t68 * t44 - t46 * t93 + t26;
t8 = -t43 * t93 - t105;
t59 = t43 * t8 + t46 * t7;
t58 = t43 * t7 - t46 * t8;
t19 = -t44 * t101 + t42 * t47;
t11 = t46 * t102 - t19 * t43;
t12 = t43 * t102 + t19 * t46;
t57 = t11 * t46 + t12 * t43;
t56 = t11 * t43 - t12 * t46;
t10 = t19 * qJD(4) - t47 * t35;
t55 = t10 * t43 + t18 * t88;
t54 = -t10 * t46 + t75;
t53 = t71 + t78;
t22 = t43 * t83 + t44 * t88;
t23 = t43 * t85 - t76;
t52 = t105 * qJD(5) + t46 * t27;
t9 = -t44 * t35 + t91;
t3 = -t12 * qJD(5) + t9 * t43 + t46 * t74;
t4 = t11 * qJD(5) + t43 * t74 - t9 * t46;
t51 = -t57 * qJD(5) - t3 * t43 + t4 * t46;
t14 = t43 * t67 + t84;
t15 = -t43 * qJD(6) + t46 * t67;
t32 = t98 * t43;
t33 = t98 * t46;
t50 = t14 * t46 - t15 * t43 + (-t32 * t46 + t33 * t43) * qJD(5);
t36 = -t46 * pkin(5) - pkin(4);
t28 = (pkin(5) * t43 - t49) * t47;
t20 = t44 * t89 - t46 * t83;
t13 = -t23 * pkin(5) + t72;
t6 = t52 - t63;
t5 = t44 * t77 - t80;
t2 = -qJ(6) * t76 + (-qJD(6) * t47 + (qJ(6) * qJD(4) - t86) * t44) * t43 + t80;
t1 = qJ(6) * t71 + (qJ(6) * t89 + t68 * qJD(4) - t84) * t47 + t52;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t18 * t10 + 0.2e1 * t11 * t3 + 0.2e1 * t12 * t4; 0, 0, -t35, -t74, t35, t74 (qJD(3) * t45 + (-pkin(2) * t45 + qJ(3) * t48) * qJD(2)) * t41, 0, 0, 0, 0, 0 (t44 * t92 + t45 * t83) * t41 (-t45 * t85 + t47 * t92) * t41, 0, 0, 0, 0, 0 (-t43 * t91 + t3) * t44 + (qJD(4) * t11 + t55) * t47 (-t18 * t90 - t4) * t44 + (-qJD(4) * t12 - t54) * t47, t57 * t85 + (t56 * qJD(5) - t3 * t46 - t4 * t43) * t47, t11 * t1 + t10 * t28 + t12 * t2 + t18 * t13 + t3 * t7 + t4 * t8; 0, 0, 0, 0, 0, t104, qJ(3) * t104, -0.2e1 * t70, 0.2e1 * t65, 0, 0, 0, 0.2e1 * qJD(3) * t44 + 0.2e1 * t47 * t82, 0.2e1 * qJD(3) * t47 - 0.2e1 * t44 * t82, -0.2e1 * t39 * t70 - 0.2e1 * t40 * t73, 0.2e1 * t40 * t66 + 0.4e1 * t47 * t62, -0.2e1 * t44 * t78 - 0.2e1 * t95 * t90, 0.2e1 * t43 * t65 - 0.2e1 * t44 * t76, t64, -0.2e1 * t40 * t46 * t86 + 0.2e1 * t26 * t83 + 0.2e1 * (t6 + t63) * t44, 0.2e1 * t40 * t77 + 0.2e1 * t5 * t44 + 0.2e1 * (-t25 + t34) * t83, 0.2e1 * t59 * t85 + 0.2e1 * (t58 * qJD(5) - t1 * t46 - t2 * t43) * t47, 0.2e1 * t7 * t1 + 0.2e1 * t28 * t13 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t56 * qJD(4) - t10) * t47 + (t51 + t91) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * t88, t94 * t89, 0 (-t58 * qJD(4) - t13) * t47 + (qJD(4) * t28 - t59 * qJD(5) - t1 * t43 + t2 * t46) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t96) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, t54, t55, t51, pkin(5) * t75 + t10 * t36 + t11 * t15 + t12 * t14 + t3 * t32 - t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t83, 0, -t72, -t69, -t47 * t66 - t62, -0.4e1 * t47 * t73 + t97 * t85, t22, -t20, 0 (-t43 * t99 - t61 * t46) * qJD(5) + (t60 * t43 - t34) * qJD(4) (t61 * t43 - t46 * t99) * qJD(5) + (-t46 * t103 + (pkin(4) * t46 + t100) * t44) * qJD(4) (t32 * t85 - t15 * t47 + t2 + (t33 * t47 - t7) * qJD(5)) * t46 + (-t33 * t85 - t14 * t47 - t1 + (t32 * t47 - t8) * qJD(5)) * t43, t1 * t32 + t13 * t36 + t8 * t14 + t7 * t15 - t2 * t33 + t28 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t83, 0, 0, 0, 0, 0, -t53, t23, t96 * t83 (-t79 + (-t32 * t43 - t33 * t46) * qJD(4)) * t47 + (qJD(4) * t36 + t50) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t73, -0.2e1 * t66, 0, 0, 0, t43 * t81, t46 * t81, 0.2e1 * t50, -0.2e1 * t33 * t14 + 0.2e1 * t32 * t15 + 0.2e1 * t36 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t23, t83, t6, t5, t53 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t20, 0, -t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t89, 0, -pkin(9) * t88, pkin(9) * t89, -pkin(5) * t88, t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t16;
