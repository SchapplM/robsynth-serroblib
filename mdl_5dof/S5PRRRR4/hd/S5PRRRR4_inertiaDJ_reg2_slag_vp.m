% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR4
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:08:01
% DurationCPUTime: 0.71s
% Computational Cost: add. (766->101), mult. (1860->175), div. (0->0), fcn. (1491->6), ass. (0->84)
t64 = sin(qJ(3));
t113 = t64 * pkin(2);
t87 = pkin(7) + t113;
t115 = qJD(4) + qJD(5);
t114 = pkin(8) + pkin(7);
t66 = cos(qJ(3));
t112 = t66 * pkin(2);
t111 = cos(qJ(5));
t62 = sin(qJ(5));
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t44 = t111 * t63 + t62 * t65;
t27 = t115 * t44;
t108 = t62 * t63;
t86 = t111 * t65;
t43 = -t86 + t108;
t110 = t43 * t27;
t82 = t111 * qJD(5);
t26 = -qJD(4) * t86 + t115 * t108 - t65 * t82;
t109 = t44 * t26;
t98 = t63 * qJD(4);
t93 = pkin(4) * t98;
t100 = pkin(2) * qJD(3);
t95 = t64 * t100;
t45 = t93 + t95;
t57 = -t65 * pkin(4) - pkin(3);
t48 = t57 - t112;
t106 = t48 * t27 + t45 * t43;
t105 = -t48 * t26 + t45 * t44;
t104 = t57 * t27 + t43 * t93;
t103 = -t57 * t26 + t44 * t93;
t56 = -pkin(3) - t112;
t58 = t65 * qJD(4);
t91 = t63 * t100;
t102 = t56 * t58 + t64 * t91;
t60 = t63 ^ 2;
t61 = t65 ^ 2;
t101 = t60 + t61;
t99 = qJD(5) * t62;
t97 = pkin(3) * t98;
t96 = pkin(3) * t58;
t94 = t66 * t100;
t92 = pkin(4) * t99;
t90 = t65 * t100;
t89 = t62 * t114;
t88 = t63 * t58;
t85 = t101 * t66;
t59 = t65 * pkin(8);
t40 = t65 * t87 + t59;
t79 = -pkin(8) - t87;
t74 = t79 * t63;
t69 = t111 * t74;
t24 = -t62 * t40 + t69;
t70 = t62 * t74;
t25 = t111 * t40 + t70;
t71 = qJD(4) * t79;
t80 = t66 * t90;
t67 = t63 * t71 + t80;
t81 = t66 * t91;
t68 = t65 * t71 - t81;
t5 = -qJD(5) * t69 - t111 * t67 + t40 * t99 - t62 * t68;
t6 = -qJD(5) * t70 + t111 * t68 - t40 * t82 - t62 * t67;
t84 = t24 * t26 - t25 * t27 + t5 * t43 - t6 * t44;
t49 = t65 * pkin(7) + t59;
t76 = t114 * t111;
t73 = t63 * t76;
t12 = t115 * t73 + t49 * t99 + t89 * t58;
t78 = t63 * t89;
t31 = t111 * t49 - t78;
t13 = -t31 * qJD(5) + (-t65 * t76 + t78) * qJD(4);
t30 = -t62 * t49 - t73;
t83 = t12 * t43 - t13 * t44 + t30 * t26 - t31 * t27;
t77 = pkin(4) * t82;
t75 = t87 * qJD(4);
t72 = t56 * t98 - t64 * t90;
t52 = -0.2e1 * t88;
t51 = 0.2e1 * t88;
t42 = 0.2e1 * (-t60 + t61) * qJD(4);
t37 = t85 * t100;
t17 = -0.2e1 * t109;
t16 = 0.2e1 * t110;
t7 = 0.2e1 * t43 * t26 - 0.2e1 * t44 * t27;
t4 = (t111 * t26 - t27 * t62 + (-t111 * t43 + t44 * t62) * qJD(5)) * pkin(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t109 + 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t24 - t26 * t25 - t43 * t6 - t44 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t94, 0, 0, t51, t42, 0, t52, 0, 0, 0.2e1 * t72, 0.2e1 * t102, 0.2e1 * t37, 0.2e1 * (t87 * t112 * t101 + t56 * t113) * qJD(3), t17, t7, 0, t16, 0, 0, 0.2e1 * t106, 0.2e1 * t105, 0.2e1 * t84, 0.2e1 * t24 * t6 - 0.2e1 * t25 * t5 + 0.2e1 * t48 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t12 - t43 * t13 - t26 * t31 - t27 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94, 0, 0, t51, t42, 0, t52, 0, 0, t72 - t97, -t96 + t102, t37, (-pkin(3) * t64 + pkin(7) * t85) * t100, t17, t7, 0, t16, 0, 0, t104 + t106, t103 + t105, t83 + t84, -t25 * t12 + t24 * t13 + t6 * t30 - t5 * t31 + t45 * t57 + t48 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t42, 0, t52, 0, 0, -0.2e1 * t97, -0.2e1 * t96, 0, 0, t17, t7, 0, t16, 0, 0, 0.2e1 * t104, 0.2e1 * t103, 0.2e1 * t83, -0.2e1 * t31 * t12 + 0.2e1 * t30 * t13 + 0.2e1 * t57 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, (-t111 * t27 - t26 * t62 + (t111 * t44 + t43 * t62) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t98, 0, -t65 * t75 - t81, t63 * t75 - t80, 0, 0, 0, 0, -t26, 0, -t27, 0, t6, t5, t4, (t111 * t6 - t5 * t62 + (t111 * t25 - t24 * t62) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t98, 0, -pkin(7) * t58, pkin(7) * t98, 0, 0, 0, 0, -t26, 0, -t27, 0, t13, t12, t4, (t111 * t13 - t12 * t62 + (t111 * t31 - t30 * t62) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, 0, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, 0, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
