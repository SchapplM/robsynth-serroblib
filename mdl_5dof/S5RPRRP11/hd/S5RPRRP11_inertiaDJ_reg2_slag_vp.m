% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:36
% DurationCPUTime: 1.08s
% Computational Cost: add. (1519->133), mult. (3515->250), div. (0->0), fcn. (3263->6), ass. (0->90)
t63 = sin(qJ(4));
t59 = t63 ^ 2;
t65 = cos(qJ(4));
t60 = t65 ^ 2;
t98 = t59 - t60;
t48 = t98 * qJD(4);
t109 = cos(qJ(3));
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t64 = sin(qJ(3));
t43 = -t109 * t62 + t64 * t61;
t44 = t109 * t61 + t64 * t62;
t53 = -t62 * pkin(2) - pkin(1);
t30 = t43 * pkin(3) - t44 * pkin(7) + t53;
t99 = pkin(6) + qJ(2);
t46 = t99 * t61;
t47 = t99 * t62;
t34 = t109 * t47 - t64 * t46;
t118 = t63 * t30 + t65 * t34;
t33 = t109 * t46 + t64 * t47;
t19 = t43 * qJD(2) + t33 * qJD(3);
t39 = t43 * qJD(3);
t40 = t44 * qJD(3);
t83 = t40 * pkin(3) + t39 * pkin(7);
t4 = -qJD(4) * t118 + t63 * t19 + t65 * t83;
t11 = t65 * t30 - t63 * t34;
t55 = qJD(4) * t63;
t56 = qJD(4) * t65;
t3 = t65 * t19 - t30 * t56 + t34 * t55 - t63 * t83;
t117 = t3 * t63 - t4 * t65 + (t11 * t63 - t118 * t65) * qJD(4);
t76 = t65 * pkin(4) + t63 * qJ(5);
t116 = t76 * qJD(4) - t65 * qJD(5);
t95 = t43 * qJD(5);
t97 = t40 * qJ(5);
t1 = -t3 + t95 + t97;
t110 = t40 * pkin(4);
t2 = -t110 - t4;
t8 = t43 * qJ(5) + t118;
t9 = -t43 * pkin(4) - t11;
t115 = t1 * t63 - t2 * t65 + (t63 * t9 + t65 * t8) * qJD(4);
t114 = 0.2e1 * t40;
t113 = 0.2e1 * qJD(5);
t112 = pkin(7) * t40;
t111 = pkin(7) * t43;
t20 = t44 * qJD(2) + t34 * qJD(3);
t108 = t33 * t20;
t107 = t44 * t39;
t106 = t44 * t63;
t105 = t44 * t65;
t104 = t63 * t40;
t101 = t65 * t39;
t100 = t65 * t40;
t94 = t63 * qJD(5);
t32 = t43 * t114;
t92 = -0.2e1 * pkin(3) * qJD(4);
t91 = t63 * t101;
t90 = pkin(7) * t55;
t89 = pkin(7) * t56;
t88 = t63 * t56;
t41 = t44 ^ 2;
t85 = t41 * t88;
t84 = 0.2e1 * (t61 ^ 2 + t62 ^ 2) * qJD(2);
t82 = pkin(3) * t39 - t112;
t81 = pkin(3) * t44 + t111;
t77 = t63 * t8 - t65 * t9;
t75 = pkin(4) * t63 - qJ(5) * t65;
t74 = t11 * t65 + t118 * t63;
t29 = -t63 * t39 + t44 * t56;
t71 = t44 * t55 + t101;
t28 = t43 * t56 + t104;
t26 = t43 * t55 - t100;
t45 = -pkin(3) - t76;
t5 = t116 * t44 - t75 * t39 + t20;
t70 = -t5 + (t44 * t45 - t111) * qJD(4);
t13 = t75 * t44 + t33;
t38 = -pkin(4) * t55 + qJ(5) * t56 + t94;
t69 = -qJD(4) * t13 + t38 * t44 + t39 * t45 + t112;
t67 = -t77 * qJD(4) + t1 * t65 + t2 * t63;
t66 = -t74 * qJD(4) - t3 * t65 - t4 * t63;
t50 = -0.2e1 * t88;
t49 = 0.2e1 * t88;
t22 = (t59 + t60) * t39;
t17 = -0.2e1 * t60 * t107 - 0.2e1 * t85;
t16 = -0.2e1 * t59 * t107 + 0.2e1 * t85;
t15 = t44 * t48 + t91;
t14 = -t98 * t39 + 0.4e1 * t44 * t88;
t10 = t41 * t48 + 0.2e1 * t44 * t91;
t7 = t44 * t104 + t29 * t43;
t6 = 0.2e1 * t44 * t100 - 0.2e1 * t71 * t43;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, qJ(2) * t84, -0.2e1 * t107, 0.2e1 * t39 * t43 - 0.2e1 * t44 * t40, 0, t32, 0, 0, t53 * t114, -0.2e1 * t53 * t39, 0.2e1 * t19 * t43 + 0.2e1 * t20 * t44 - 0.2e1 * t33 * t39 - 0.2e1 * t34 * t40, -0.2e1 * t34 * t19 + 0.2e1 * t108, t17, 0.2e1 * t10, t6, t16, -0.2e1 * t7, t32, 0.2e1 * t20 * t106 + 0.2e1 * t11 * t40 + 0.2e1 * t29 * t33 + 0.2e1 * t4 * t43, 0.2e1 * t20 * t105 - 0.2e1 * t118 * t40 + 0.2e1 * t3 * t43 - 0.2e1 * t71 * t33, 0.2e1 * t117 * t44 + 0.2e1 * t74 * t39, 0.2e1 * t11 * t4 - 0.2e1 * t118 * t3 + 0.2e1 * t108, t17, t6, -0.2e1 * t10, t32, 0.2e1 * t7, t16, 0.2e1 * t5 * t106 + 0.2e1 * t29 * t13 - 0.2e1 * t2 * t43 - 0.2e1 * t9 * t40, -0.2e1 * t115 * t44 + 0.2e1 * t77 * t39, 0.2e1 * t1 * t43 - 0.2e1 * t5 * t105 + 0.2e1 * t71 * t13 + 0.2e1 * t8 * t40, 0.2e1 * t8 * t1 + 0.2e1 * t13 * t5 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t28, t22, -t117, 0, 0, 0, 0, 0, 0, -t26, t22, t28, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t40, 0, -t20, t19, 0, 0, -t15, -t14, t28, t15, -t26, 0, -t20 * t65 + t82 * t63 + (t33 * t63 - t81 * t65) * qJD(4), t20 * t63 + t82 * t65 + (t33 * t65 + t81 * t63) * qJD(4), t66, -t20 * pkin(3) + t66 * pkin(7), -t15, t28, t14, 0, t26, t15, -t69 * t63 + t70 * t65, t67, t70 * t63 + t69 * t65, t67 * pkin(7) - t13 * t38 + t5 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -0.2e1 * t48, 0, t50, 0, 0, t63 * t92, t65 * t92, 0, 0, t49, 0, 0.2e1 * t48, 0, 0, t50, 0.2e1 * t38 * t65 + 0.2e1 * t45 * t55, 0, 0.2e1 * t38 * t63 - 0.2e1 * t45 * t56, -0.2e1 * t45 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, -t29, t40, t4, t3, 0, 0, 0, -t71, 0, t40, t29, 0, t4 + 0.2e1 * t110, t76 * t39 + (t75 * qJD(4) - t94) * t44, -t3 + 0.2e1 * t95 + 0.2e1 * t97, -t2 * pkin(4) + t1 * qJ(5) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, t56, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t55, 0, -t89, t90, 0, 0, 0, t56, 0, 0, t55, 0, -t89, -t116, -t90, -t116 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, qJ(5) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t71, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
