% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:22
% EndTime: 2019-03-09 02:11:25
% DurationCPUTime: 0.86s
% Computational Cost: add. (669->129), mult. (1339->221), div. (0->0), fcn. (946->4), ass. (0->85)
t34 = -pkin(7) + qJ(2);
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t81 = t39 * qJD(4);
t99 = t37 * qJD(2) + t34 * t81;
t36 = sin(qJ(5));
t30 = t36 ^ 2;
t38 = cos(qJ(5));
t32 = t38 ^ 2;
t91 = t30 - t32;
t63 = t91 * qJD(5);
t95 = pkin(8) * t37;
t58 = pkin(4) * t39 + t95;
t20 = t58 * qJD(4) + qJD(3);
t35 = pkin(1) + qJ(3);
t94 = t39 * pkin(8);
t57 = t37 * pkin(4) - t94;
t21 = t57 + t35;
t29 = qJD(5) * t38;
t73 = t37 * t29;
t4 = -(qJD(5) * t21 + t99) * t36 + t38 * t20 - t34 * t73;
t52 = pkin(5) * t36 - qJ(6) * t38;
t49 = -t34 + t52;
t8 = t49 * t39;
t53 = t38 * pkin(5) + t36 * qJ(6);
t22 = -pkin(4) - t53;
t92 = t22 * t37;
t9 = t52 * qJD(5) - t36 * qJD(6);
t93 = t39 * t9;
t98 = (t92 + t94) * qJD(4) - qJD(5) * t8 - t93;
t83 = t37 * qJD(4);
t47 = t39 * qJD(2) - t34 * t83;
t97 = t58 * qJD(5) - t47;
t44 = -t53 * qJD(5) + t38 * qJD(6);
t96 = 0.2e1 * qJD(6);
t90 = t30 + t32;
t31 = t37 ^ 2;
t33 = t39 ^ 2;
t89 = t31 - t33;
t88 = t31 + t33;
t87 = qJD(4) * t8;
t86 = qJD(4) * t38;
t28 = qJD(5) * t36;
t85 = qJD(5) * t39;
t80 = qJ(2) * qJD(2);
t79 = qJ(6) * qJD(4);
t78 = -0.2e1 * pkin(4) * qJD(5);
t24 = t38 * t37 * t34;
t77 = pkin(5) * t81;
t76 = pkin(8) * t28;
t75 = pkin(8) * t29;
t74 = t36 * t85;
t72 = t38 * t85;
t71 = t34 * t28;
t70 = t36 * t29;
t69 = t38 * t83;
t68 = t38 * t81;
t67 = t37 * t81;
t65 = t39 * t79;
t64 = t90 * t39;
t62 = t89 * qJD(4);
t61 = 0.2e1 * t67;
t60 = t36 * t20 + t21 * t29 + t99 * t38;
t59 = t36 * t69;
t19 = t36 * t21;
t6 = t37 * qJ(6) + t19 + t24;
t7 = -t38 * t21 + (t34 * t36 - pkin(5)) * t37;
t56 = t36 * t7 + t38 * t6;
t55 = t36 * t6 - t38 * t7;
t5 = -t49 * t83 + (-qJD(2) - t44) * t39;
t46 = -t5 + (t22 * t39 - t95) * qJD(5);
t43 = t57 * qJD(4) - t34 * t85;
t1 = t65 + (qJD(6) - t71) * t37 + t60;
t2 = -t4 - t77;
t42 = -t56 * qJD(5) - t1 * t36 + t2 * t38;
t41 = -t55 * qJD(5) + t1 * t38 + t2 * t36;
t40 = 0.2e1 * qJD(2);
t18 = -t36 * t83 + t72;
t17 = t36 * t81 + t73;
t16 = t88 * t29;
t15 = -t69 - t74;
t14 = t37 * t28 - t68;
t13 = t88 * t28;
t3 = t37 * t71 - t60;
t10 = [0, 0, 0, 0, t40, 0.2e1 * t80, t40, 0.2e1 * qJD(3), 0.2e1 * t35 * qJD(3) + 0.2e1 * t80, -0.2e1 * t67, 0.2e1 * t62, 0, 0, 0, 0.2e1 * qJD(3) * t37 + 0.2e1 * t35 * t81, 0.2e1 * qJD(3) * t39 - 0.2e1 * t35 * t83, -0.2e1 * t32 * t67 - 0.2e1 * t33 * t70, 0.2e1 * t33 * t63 + 0.4e1 * t39 * t59, -0.2e1 * t37 * t74 - 0.2e1 * t89 * t86, 0.2e1 * t36 * t62 - 0.2e1 * t37 * t72, t61, 0.2e1 * t21 * t68 - 0.2e1 * t33 * qJD(2) * t36 + 0.2e1 * t4 * t37 + 0.2e1 * (-t33 * t29 + t36 * t67) * t34, 0.2e1 * t3 * t37 + 0.2e1 * (-qJD(2) * t38 + t71) * t33 + 0.2e1 * (-t19 + t24) * t81, 0.2e1 * (-t36 * t87 - t2) * t37 + 0.2e1 * (-qJD(4) * t7 + t8 * t29 + t5 * t36) * t39, 0.2e1 * t42 * t39 + 0.2e1 * t55 * t83, 0.2e1 * (t8 * t86 + t1) * t37 + 0.2e1 * (qJD(4) * t6 + t8 * t28 - t5 * t38) * t39, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t81, t83, 0, 0, 0, 0, 0, t14, t17, t14, -t90 * t83, -t17, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t13, -t16, 0, -t13 (t56 * qJD(4) - t5) * t39 + (t41 + t87) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t90) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t81, 0, t47, -t99, -t39 * t63 - t59, -0.4e1 * t39 * t70 + t91 * t83, t17, -t14, 0, t43 * t36 - t97 * t38, t97 * t36 + t43 * t38, -t98 * t36 + t46 * t38, t41, t46 * t36 + t98 * t38, t41 * pkin(8) + t5 * t22 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t81, 0, 0, 0, 0, 0, t15, -t18, t15, qJD(4) * t64, t18, -t93 + (pkin(8) * t64 + t92) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t70, -0.2e1 * t63, 0, 0, 0, t36 * t78, t38 * t78, 0.2e1 * t22 * t28 - 0.2e1 * t9 * t38, 0, -0.2e1 * t22 * t29 - 0.2e1 * t9 * t36, 0.2e1 * t22 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t18, t81, t4, t3, t4 + 0.2e1 * t77 (pkin(5) * t83 - qJ(6) * t85) * t38 + (t37 * t79 + (pkin(5) * qJD(5) - qJD(6)) * t39) * t36, 0.2e1 * t65 + (t96 - t71) * t37 + t60, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, t28, 0, -t29, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t14, -t17, 0, -t14, t44 * t37 - t52 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, -t75, t76, -t75, t44, -t76, t44 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, qJ(6) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
