% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:05
% DurationCPUTime: 1.22s
% Computational Cost: add. (1681->154), mult. (4301->287), div. (0->0), fcn. (3981->8), ass. (0->81)
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t101 = -pkin(7) - pkin(6);
t60 = sin(qJ(3));
t83 = t101 * t60;
t72 = t62 * t83;
t63 = cos(qJ(3));
t84 = t101 * t63;
t30 = t59 * t84 + t72;
t44 = t59 * t63 + t60 * t62;
t61 = sin(qJ(2));
t34 = t44 * t61;
t102 = qJD(3) + qJD(4);
t100 = cos(qJ(5));
t99 = t62 * t63;
t45 = t59 * t83;
t31 = -t62 * t84 + t45;
t56 = t60 ^ 2;
t57 = t63 ^ 2;
t98 = t56 + t57;
t97 = qJD(4) * t59;
t96 = qJD(4) * t62;
t58 = sin(qJ(5));
t95 = qJD(5) * t58;
t94 = t60 * qJD(3);
t93 = t61 * qJD(2);
t92 = t63 * qJD(3);
t64 = cos(qJ(2));
t91 = t64 * qJD(2);
t90 = t58 * t59 * pkin(3);
t89 = -0.2e1 * pkin(2) * qJD(3);
t88 = pkin(3) * t97;
t87 = pkin(3) * t96;
t86 = pkin(4) * t95;
t85 = pkin(3) * t94;
t50 = t59 * t92;
t82 = t60 * t92;
t81 = t61 * t91;
t80 = t60 * t91;
t79 = t64 * t94;
t78 = t63 * t91;
t55 = -pkin(3) * t63 - pkin(2);
t77 = t100 * t59;
t76 = t98 * t64;
t75 = qJD(5) * t100;
t74 = t60 * t96 + t62 * t94 + t63 * t97 + t50;
t71 = pkin(4) * t75;
t43 = t59 * t60 - t99;
t35 = t43 * t61;
t19 = -t100 * t35 - t58 * t34;
t28 = t100 * t44 - t58 * t43;
t54 = pkin(3) * t62 + pkin(4);
t22 = -t54 * t75 - t100 * t87 + (qJD(4) + qJD(5)) * t90;
t15 = -qJD(3) * t72 - qJD(4) * t30 - t101 * t50;
t70 = -t44 * pkin(8) + t30;
t69 = t58 * t70;
t68 = -pkin(8) * t74 - t15;
t67 = t100 * t70;
t29 = t102 * t43;
t66 = (-t59 * t75 + (-t58 * t62 - t77) * qJD(4)) * pkin(3);
t16 = -t31 * qJD(4) + (t101 * t99 - t45) * qJD(3);
t65 = t29 * pkin(8) + t16;
t38 = pkin(3) * t77 + t58 * t54;
t37 = t100 * t54 - t90;
t33 = pkin(4) * t43 + t55;
t27 = t100 * t43 + t58 * t44;
t23 = -t54 * t95 + t66;
t21 = pkin(4) * t74 + t85;
t20 = -pkin(8) * t43 + t31;
t18 = -t100 * t34 + t58 * t35;
t14 = t29 * t61 - t44 * t91;
t13 = t102 * t34 + t59 * t80 - t62 * t78;
t10 = t100 * t20 + t69;
t9 = -t58 * t20 + t67;
t6 = qJD(5) * t28 + t100 * t74 - t58 * t29;
t5 = t100 * t29 + t43 * t75 + t44 * t95 + t58 * t74;
t4 = -qJD(5) * t19 + t100 * t14 + t58 * t13;
t3 = t100 * t13 - t58 * t14 + t34 * t75 - t35 * t95;
t2 = -qJD(5) * t69 + t100 * t65 - t20 * t75 - t58 * t68;
t1 = -qJD(5) * t67 - t100 * t68 + t20 * t95 - t58 * t65;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t98) * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t35 - 0.2e1 * t14 * t34 - 0.2e1 * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t18 * t4 - 0.2e1 * t19 * t3 - 0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t91, 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t93 - t79, t60 * t93 - t64 * t92, qJD(2) * t76, (-pkin(2) * t61 + pkin(6) * t76) * qJD(2), 0, 0, 0, 0, 0, 0, t43 * t93 - t64 * t74, t29 * t64 + t44 * t93, t13 * t43 - t14 * t44 - t34 * t29 + t35 * t74, -pkin(3) * t79 - t13 * t31 + t14 * t30 + t15 * t35 - t16 * t34 + t55 * t93, 0, 0, 0, 0, 0, 0, t27 * t93 - t6 * t64, t28 * t93 + t5 * t64, t18 * t5 - t19 * t6 + t27 * t3 - t28 * t4, -t1 * t19 - t10 * t3 + t18 * t2 - t21 * t64 + t33 * t93 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t82, 0.2e1 * (-t56 + t57) * qJD(3), 0, -0.2e1 * t82, 0, 0, t60 * t89, t63 * t89, 0, 0, -0.2e1 * t44 * t29, 0.2e1 * t29 * t43 - 0.2e1 * t44 * t74, 0, 0.2e1 * t43 * t74, 0, 0, 0.2e1 * t43 * t85 + 0.2e1 * t55 * t74, -0.2e1 * t29 * t55 + 0.2e1 * t44 * t85, 0.2e1 * t15 * t43 - 0.2e1 * t16 * t44 + 0.2e1 * t30 * t29 - 0.2e1 * t31 * t74, -0.2e1 * t15 * t31 + 0.2e1 * t16 * t30 + 0.2e1 * t55 * t85, -0.2e1 * t28 * t5, 0.2e1 * t27 * t5 - 0.2e1 * t28 * t6, 0, 0.2e1 * t27 * t6, 0, 0, 0.2e1 * t21 * t27 + 0.2e1 * t33 * t6, 0.2e1 * t21 * t28 - 0.2e1 * t33 * t5, 0.2e1 * t1 * t27 - 0.2e1 * t10 * t6 - 0.2e1 * t2 * t28 + 0.2e1 * t5 * t9, -0.2e1 * t1 * t10 + 0.2e1 * t2 * t9 + 0.2e1 * t21 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t92 - t80, t61 * t94 - t78, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, (-t13 * t59 + t14 * t62 + (t34 * t59 - t35 * t62) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t4, t3, 0, t18 * t23 - t19 * t22 - t3 * t38 + t37 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, -t94, 0, -pkin(6) * t92, pkin(6) * t94, 0, 0, 0, 0, -t29, 0, -t74, 0, t16, t15, (-t59 * t74 + t62 * t29 + (-t43 * t62 + t44 * t59) * qJD(4)) * pkin(3), (-t15 * t59 + t16 * t62 + (-t30 * t59 + t31 * t62) * qJD(4)) * pkin(3), 0, 0, -t5, 0, -t6, 0, t2, t1, t22 * t27 - t23 * t28 + t37 * t5 - t38 * t6, -t1 * t38 - t10 * t22 + t2 * t37 + t23 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, -0.2e1 * t87, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t22, 0, -0.2e1 * t22 * t38 + 0.2e1 * t23 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, (t100 * t4 - t3 * t58 + (t100 * t19 - t18 * t58) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t74, 0, t16, t15, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, (t100 * t5 - t58 * t6 + (-t100 * t27 + t28 * t58) * qJD(5)) * pkin(4), (t100 * t2 - t1 * t58 + (t10 * t100 - t58 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t87, 0, 0, 0, 0, 0, 0, 0, 0, (-pkin(4) - t54) * t95 + t66, -t71 + t22, 0, (t100 * t23 - t22 * t58 + (t100 * t38 - t37 * t58) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
