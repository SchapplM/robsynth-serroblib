% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:48
% EndTime: 2019-03-09 02:03:51
% DurationCPUTime: 0.84s
% Computational Cost: add. (774->141), mult. (1533->233), div. (0->0), fcn. (1118->6), ass. (0->90)
t29 = sin(pkin(9)) * pkin(1) + qJ(3);
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t95 = t40 * pkin(8);
t55 = t38 * pkin(4) - t95;
t21 = t29 + t55;
t37 = sin(qJ(5));
t11 = t37 * t21;
t28 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t39 = cos(qJ(5));
t24 = t39 * t38 * t28;
t48 = t24 + t11;
t32 = t37 ^ 2;
t34 = t39 ^ 2;
t88 = t32 + t34;
t89 = t32 - t34;
t59 = t89 * qJD(5);
t78 = qJ(6) * qJD(4);
t102 = (pkin(5) * qJD(5) - qJD(6)) * t40 + t38 * t78;
t50 = pkin(5) * t37 - qJ(6) * t39;
t47 = -t28 + t50;
t81 = t38 * qJD(4);
t51 = t39 * pkin(5) + t37 * qJ(6);
t99 = t51 * qJD(5) - t39 * qJD(6);
t5 = t99 * t40 - t47 * t81;
t6 = t38 * qJ(6) + t48;
t92 = t28 * t37;
t62 = -pkin(5) + t92;
t7 = -t39 * t21 + t62 * t38;
t54 = t37 * t7 + t39 * t6;
t101 = t54 * qJD(4) - t5;
t9 = t47 * t40;
t12 = t50 * qJD(5) - t37 * qJD(6);
t91 = t40 * t12;
t25 = -pkin(4) - t51;
t94 = t25 * t38;
t100 = (t94 + t95) * qJD(4) - qJD(5) * t9 - t91;
t96 = pkin(8) * t38;
t56 = pkin(4) * t40 + t96;
t23 = t56 * qJD(4) + qJD(3);
t44 = -t48 * qJD(5) + t39 * t23;
t98 = 0.2e1 * qJD(3);
t97 = 0.2e1 * qJD(6);
t93 = t25 * t40;
t90 = t40 * t28;
t33 = t38 ^ 2;
t35 = t40 ^ 2;
t87 = t33 - t35;
t86 = t33 + t35;
t85 = qJD(4) * t37;
t84 = qJD(4) * t39;
t83 = qJD(5) * t37;
t30 = qJD(5) * t39;
t82 = qJD(5) * t40;
t79 = t40 * qJD(4);
t77 = -0.2e1 * pkin(4) * qJD(5);
t65 = t28 * t79;
t76 = t21 * t30 + t37 * t23 + t39 * t65;
t75 = pkin(8) * t83;
t74 = pkin(8) * t30;
t73 = t28 * t83;
t72 = t38 * t83;
t71 = t37 * t82;
t70 = t39 * t82;
t69 = t37 * t30;
t68 = t39 * t81;
t67 = t39 * t79;
t66 = t38 * t79;
t63 = t40 * t78;
t61 = t88 * t38;
t60 = t88 * t40;
t57 = t37 * t68;
t53 = t37 * t6 - t39 * t7;
t46 = pkin(5) * t81 - qJ(6) * t82;
t45 = -t5 + (t93 - t96) * qJD(5);
t1 = t63 + (qJD(6) - t73) * t38 + t76;
t2 = t62 * t79 - t44;
t42 = -t53 * qJD(5) + t1 * t39 + t2 * t37;
t41 = qJD(4) * t9 + t42;
t31 = qJD(4) * t33;
t20 = -t37 * t81 + t70;
t19 = t38 * t30 + t37 * t79;
t18 = t86 * t30;
t17 = t68 + t71;
t16 = -t67 + t72;
t15 = t86 * t83;
t8 = (-0.1e1 + t88) * t66;
t4 = -t37 * t65 + t44;
t3 = t28 * t72 - t76;
t10 = [0, 0, 0, 0, 0, t98, t29 * t98, -0.2e1 * t66, -0.2e1 * t35 * qJD(4) + 0.2e1 * t31, 0, 0, 0, 0.2e1 * qJD(3) * t38 + 0.2e1 * t29 * t79, 0.2e1 * qJD(3) * t40 - 0.2e1 * t29 * t81, -0.2e1 * t34 * t66 - 0.2e1 * t35 * t69, 0.2e1 * t35 * t59 + 0.4e1 * t40 * t57, -0.2e1 * t38 * t71 - 0.2e1 * t87 * t84, -0.2e1 * t38 * t70 + 0.2e1 * t87 * t85, 0.2e1 * t66, 0.2e1 * t21 * t67 + 0.2e1 * t4 * t38 + 0.2e1 * (-t35 * t30 + t37 * t66) * t28, 0.2e1 * t35 * t73 + 0.2e1 * t3 * t38 + 0.2e1 * (-t11 + t24) * t79, 0.2e1 * (-t9 * t85 - t2) * t38 + 0.2e1 * (-qJD(4) * t7 + t9 * t30 + t5 * t37) * t40, 0.2e1 * t53 * t81 + 0.2e1 * (-t54 * qJD(5) - t1 * t37 + t2 * t39) * t40, 0.2e1 * (t9 * t84 + t1) * t38 + 0.2e1 * (qJD(4) * t6 - t5 * t39 + t9 * t83) * t40, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t38 + t41 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t15, -t18, 0, -t15, t101 * t40 + t41 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + (-t88 * t87 - t35) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t79, 0, -t28 * t81, -t65, -t40 * t59 - t57, -0.4e1 * t40 * t69 + t89 * t81, t19, -t16, 0 (-t37 * t90 - t56 * t39) * qJD(5) + (t55 * t37 - t24) * qJD(4) (t56 * t37 - t39 * t90) * qJD(5) + (-t39 * t95 + (pkin(4) * t39 + t92) * t38) * qJD(4), -t100 * t37 + t45 * t39, t42, t100 * t39 + t45 * t37, t42 * pkin(8) + t9 * t12 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t81, 0, 0, 0, 0, 0, t16, t19, t16, -qJD(4) * t61, -t19, t38 * t12 + (-pkin(8) * t61 + t93) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t79, 0, 0, 0, 0, 0, -t17, -t20, -t17, qJD(4) * t60, t20, -t91 + (pkin(8) * t60 + t94) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t59, 0, 0, 0, t37 * t77, t39 * t77, -0.2e1 * t12 * t39 + 0.2e1 * t25 * t83, 0, -0.2e1 * t12 * t37 - 0.2e1 * t25 * t30, 0.2e1 * t25 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t20, t79, t4, t3 (0.2e1 * pkin(5) - t92) * t79 + t44, t102 * t37 + t46 * t39, 0.2e1 * t63 + (t97 - t73) * t38 + t76, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t17, -t20, 0, -t17, -t102 * t39 + t46 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t16, -t19, 0, -t16, -t38 * t99 - t50 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t83, 0, -t74, t75, -t74, -t99, -t75, -t99 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(6) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
