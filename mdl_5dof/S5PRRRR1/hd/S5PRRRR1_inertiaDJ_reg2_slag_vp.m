% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaDJ_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:23
% EndTime: 2019-12-05 17:03:27
% DurationCPUTime: 1.02s
% Computational Cost: add. (645->103), mult. (2131->226), div. (0->0), fcn. (2068->8), ass. (0->91)
t41 = sin(qJ(5));
t37 = t41 ^ 2;
t45 = cos(qJ(5));
t39 = t45 ^ 2;
t86 = t37 - t39;
t103 = t86 * qJD(5);
t106 = t37 + t39;
t99 = cos(qJ(4));
t65 = t99 * qJD(4);
t105 = t99 * qJD(3) + t65;
t44 = sin(qJ(2));
t46 = cos(qJ(3));
t42 = sin(qJ(4));
t43 = sin(qJ(3));
t90 = t42 * t43;
t55 = t99 * t46 - t90;
t19 = t55 * t44;
t47 = cos(qJ(2));
t17 = t45 * t19 - t47 * t41;
t59 = t41 * t19 + t47 * t45;
t104 = t17 * t45 + t41 * t59;
t102 = qJD(3) + qJD(4);
t101 = 2 * pkin(2);
t48 = pkin(2) ^ 2;
t100 = 2 * t48;
t68 = t99 * t43;
t80 = t43 * qJD(3);
t71 = t44 * t80;
t77 = t47 * qJD(2);
t84 = qJD(4) * t42;
t10 = t68 * t77 - t42 * t71 - t43 * t44 * t84 + (t105 * t44 + t42 * t77) * t46;
t98 = t10 * t41;
t97 = t10 * t45;
t14 = t102 * t90 - t105 * t46;
t96 = t14 * t46;
t27 = t42 * t46 + t68;
t18 = t27 * t44;
t93 = t18 * t10;
t92 = t27 * t14;
t15 = t102 * t27;
t91 = t41 * t15;
t89 = t45 * t14;
t88 = t45 * t15;
t38 = t43 ^ 2;
t40 = t46 ^ 2;
t85 = t38 + t40;
t82 = qJD(5) * t41;
t36 = qJD(5) * t45;
t81 = qJD(5) * t46;
t79 = t44 * qJD(2);
t78 = t46 * qJD(3);
t76 = -0.2e1 * t55 * t15;
t75 = t41 * t89;
t74 = pkin(2) * t84;
t73 = t41 * t36;
t72 = t43 * t78;
t70 = t44 * t77;
t69 = t99 * t10;
t67 = qJD(5) * t99;
t64 = -0.2e1 * t72;
t24 = t27 ^ 2;
t63 = t24 * t73;
t62 = pkin(2) * t65;
t61 = t106 * t80;
t60 = t17 * t41 - t45 * t59;
t58 = -t41 * t14 + t27 * t36;
t57 = -t27 * t82 - t89;
t56 = t99 * t27 - t42 * t55;
t54 = -t46 * t79 - t47 * t80;
t53 = t41 * t81 + t45 * t80;
t52 = -t41 * t80 + t45 * t81;
t51 = t106 * t99;
t9 = t15 * t44 - t55 * t77;
t2 = t59 * qJD(5) - t41 * t79 + t45 * t9;
t3 = -t17 * qJD(5) + t41 * t9 + t45 * t79;
t50 = -t104 * qJD(5) + t2 * t41 - t3 * t45;
t1 = -t60 * qJD(5) - t2 * t45 - t3 * t41;
t49 = t99 * t14 - t15 * t42 + (t27 * t42 + t55 * t99) * qJD(4);
t31 = -0.2e1 * t73;
t30 = 0.2e1 * t73;
t25 = -0.2e1 * t103;
t23 = (-t41 * t67 - t45 * t84) * pkin(2);
t22 = (t41 * t84 - t45 * t67) * pkin(2);
t21 = t51 * qJD(4) * pkin(2);
t12 = -t36 * t55 + t91;
t11 = t55 * t82 + t88;
t7 = t27 * t103 + t75;
t6 = t18 * t82 - t97;
t5 = t18 * t36 + t98;
t4 = t86 * t14 - 0.4e1 * t27 * t73;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t85) * t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t19 * t9 - 0.2e1 * t70 + 0.2e1 * t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t2 - 0.2e1 * t3 * t59 + 0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t77, 0, 0, 0, 0, 0, 0, 0, 0, t54, t43 * t79 - t47 * t78, t85 * t77, 0, 0, 0, 0, 0, 0, 0, -t47 * t15 - t55 * t79, t47 * t14 + t27 * t79, t10 * t27 - t18 * t14 - t19 * t15 - t55 * t9, t54 * pkin(2), 0, 0, 0, 0, 0, 0, -t15 * t59 + t58 * t18 + t27 * t98 - t3 * t55, -t17 * t15 + t57 * t18 - t2 * t55 + t27 * t97, t60 * t14 + t50 * t27, (t50 * t46 + t60 * t80) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72, 0.2e1 * (-t38 + t40) * qJD(3), 0, t64, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t14 * t55 - 0.2e1 * t27 * t15, 0, t76, 0, 0, (-t15 * t46 - t55 * t80) * t101, (t27 * t80 + t96) * t101, 0, t48 * t64, -0.2e1 * t39 * t92 - 0.2e1 * t63, 0.2e1 * t24 * t103 + 0.4e1 * t27 * t75, 0.2e1 * t27 * t88 - 0.2e1 * t55 * t57, -0.2e1 * t37 * t92 + 0.2e1 * t63, -0.2e1 * t27 * t91 + 0.2e1 * t55 * t58, t76, (-t46 * t88 - t53 * t55) * t101, (t46 * t91 - t52 * t55) * t101, (-t106 * t96 - t27 * t61) * t101, -t46 * t61 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t77 - t44 * t78, -t46 * t77 + t71, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, (-t69 - t42 * t9 + (t18 * t42 + t99 * t19) * qJD(4)) * pkin(2), 0, 0, 0, 0, 0, 0, t6, t5, t1, (-t69 + (qJD(4) * t18 + t1) * t42 + t104 * t65) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, -t80, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t15, 0, 0, 0, t49 * pkin(2), 0, -t7, t4, t12, t7, t11, 0, (-t56 * t36 + t49 * t41) * pkin(2), (t49 * t45 + t56 * t82) * pkin(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t74, -0.2e1 * t62, 0, 0, t30, t25, 0, t31, 0, 0, 0.2e1 * t23, 0.2e1 * t22, 0.2e1 * t21, (-t99 + t51) * t84 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, t1, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t15, 0, 0, 0, 0, 0, -t7, t4, t12, t7, t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t62, 0, 0, t30, t25, 0, t31, 0, 0, t23, t22, t21, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t25, 0, t31, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t58, t15, t53 * pkin(2), t52 * pkin(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t82, 0, (-t42 * t36 - t41 * t65) * pkin(2), (t42 * t82 - t45 * t65) * pkin(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t82, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
