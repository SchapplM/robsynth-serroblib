% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:14
% EndTime: 2019-12-05 18:12:19
% DurationCPUTime: 1.11s
% Computational Cost: add. (2968->124), mult. (6451->234), div. (0->0), fcn. (6620->8), ass. (0->70)
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t51 = t66 * t60 + t64 * t61;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t75 = t64 * t60 - t66 * t61;
t99 = -t65 * t51 + t63 * t75;
t46 = t51 * qJD(3);
t98 = 0.2e1 * t46;
t91 = pkin(6) + qJ(2);
t53 = t91 * t60;
t54 = t91 * t61;
t38 = -t64 * t53 + t66 * t54;
t31 = -pkin(7) * t75 + t38;
t97 = -t66 * t53 - t64 * t54;
t73 = -t51 * pkin(7) + t97;
t15 = t65 * t31 + t63 * t73;
t26 = t75 * qJD(2) - qJD(3) * t97;
t45 = t75 * qJD(3);
t96 = -0.2e1 * t45;
t95 = t46 * pkin(3);
t94 = cos(qJ(5));
t90 = qJD(4) * t63;
t89 = qJD(4) * t65;
t62 = sin(qJ(5));
t88 = qJD(5) * t62;
t87 = t62 * t63 * pkin(3);
t86 = pkin(3) * t90;
t85 = pkin(3) * t89;
t84 = pkin(4) * t88;
t56 = -t61 * pkin(2) - pkin(1);
t83 = t94 * t63;
t29 = t65 * t73;
t14 = -t63 * t31 + t29;
t82 = t94 * qJD(5);
t81 = pkin(4) * t82;
t80 = 0.2e1 * (t60 ^ 2 + t61 ^ 2) * qJD(2);
t79 = t65 * t45 + t63 * t46;
t78 = -t63 * t51 - t65 * t75;
t74 = -t63 * t45 + t65 * t46 + t51 * t89 - t75 * t90;
t39 = pkin(3) * t75 + t56;
t72 = -pkin(8) * t99 - t14;
t27 = -t51 * qJD(2) - qJD(3) * t38;
t9 = t31 * t90 - t63 * (t45 * pkin(7) + t27) - t65 * (-t46 * pkin(7) - t26) - qJD(4) * t29;
t71 = t94 * t78;
t57 = t65 * pkin(3) + pkin(4);
t34 = -t57 * t82 - t94 * t85 + (qJD(4) + qJD(5)) * t87;
t70 = t94 * t72;
t21 = t62 * t78 - t94 * t99;
t12 = t78 * pkin(8) + t15;
t6 = t94 * t12 - t62 * t72;
t69 = t74 * pkin(8) + t9;
t68 = (-t63 * t82 + (-t62 * t65 - t83) * qJD(4)) * pkin(3);
t10 = t79 * pkin(7) - t15 * qJD(4) + (-t38 * t65 - t63 * t97) * qJD(3) + t99 * qJD(2);
t19 = -t78 * qJD(4) + t79;
t67 = t19 * pkin(8) + t10;
t44 = pkin(3) * t83 + t62 * t57;
t43 = t94 * t57 - t87;
t35 = -t57 * t88 + t68;
t24 = -t78 * pkin(4) + t39;
t20 = -t62 * t99 - t71;
t13 = t74 * pkin(4) + t95;
t8 = qJD(5) * t21 - t62 * t19 + t94 * t74;
t7 = -qJD(5) * t71 + t94 * t19 + t62 * t74 - t88 * t99;
t5 = -t62 * t12 - t70;
t2 = -t6 * qJD(5) + t62 * t69 + t94 * t67;
t1 = qJD(5) * t70 + t12 * t88 - t62 * t67 + t94 * t69;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, qJ(2) * t80, t51 * t96, 0.2e1 * t45 * t75 - 0.2e1 * t51 * t46, 0, t75 * t98, 0, 0, t56 * t98, t56 * t96, 0.2e1 * t26 * t75 - 0.2e1 * t27 * t51 - 0.2e1 * t38 * t46 + 0.2e1 * t45 * t97, -0.2e1 * t38 * t26 + 0.2e1 * t27 * t97, 0.2e1 * t99 * t19, -0.2e1 * t19 * t78 + 0.2e1 * t74 * t99, 0, -0.2e1 * t78 * t74, 0, 0, 0.2e1 * t39 * t74 - 0.2e1 * t78 * t95, -0.2e1 * t39 * t19 - 0.2e1 * t95 * t99, 0.2e1 * t10 * t99 + 0.2e1 * t14 * t19 - 0.2e1 * t15 * t74 - 0.2e1 * t9 * t78, 0.2e1 * t14 * t10 - 0.2e1 * t15 * t9 + 0.2e1 * t39 * t95, -0.2e1 * t21 * t7, 0.2e1 * t7 * t20 - 0.2e1 * t21 * t8, 0, 0.2e1 * t20 * t8, 0, 0, 0.2e1 * t13 * t20 + 0.2e1 * t24 * t8, 0.2e1 * t13 * t21 - 0.2e1 * t24 * t7, 0.2e1 * t1 * t20 - 0.2e1 * t2 * t21 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t24 * t13 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t19, 0, t95, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t46, 0, t27, t26, 0, 0, 0, 0, -t19, 0, -t74, 0, t10, t9, (-t63 * t74 + t65 * t19 + (-t63 * t99 + t65 * t78) * qJD(4)) * pkin(3), (t10 * t65 - t63 * t9 + (-t14 * t63 + t15 * t65) * qJD(4)) * pkin(3), 0, 0, -t7, 0, -t8, 0, t2, t1, t34 * t20 - t35 * t21 + t43 * t7 - t44 * t8, -t1 * t44 + t2 * t43 - t6 * t34 + t5 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86, -0.2e1 * t85, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t35, 0.2e1 * t34, 0, -0.2e1 * t44 * t34 + 0.2e1 * t43 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t74, 0, t10, t9, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, (t94 * t7 - t62 * t8 + (-t94 * t20 + t21 * t62) * qJD(5)) * pkin(4), (t94 * t2 - t1 * t62 + (-t5 * t62 + t94 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, 0, 0, 0, 0, 0, 0, 0, (-pkin(4) - t57) * t88 + t68, -t81 + t34, 0, (t94 * t35 - t34 * t62 + (-t43 * t62 + t94 * t44) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t84, -0.2e1 * t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t81, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
