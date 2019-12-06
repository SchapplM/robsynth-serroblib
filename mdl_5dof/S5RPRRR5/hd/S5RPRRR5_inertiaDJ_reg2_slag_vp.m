% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:29
% EndTime: 2019-12-05 18:16:33
% DurationCPUTime: 0.77s
% Computational Cost: add. (1095->106), mult. (2427->177), div. (0->0), fcn. (1990->8), ass. (0->82)
t117 = sin(pkin(9)) * pkin(1);
t62 = cos(pkin(9)) * pkin(1) + pkin(2);
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t102 = t74 * t117 + t72 * t62;
t88 = -t72 * t117 + t74 * t62;
t119 = qJD(4) + qJD(5);
t42 = t88 * qJD(3);
t71 = sin(qJ(4));
t67 = t71 ^ 2;
t73 = cos(qJ(4));
t68 = t73 ^ 2;
t22 = (-t68 - t67) * t42;
t118 = pkin(8) + pkin(7);
t116 = t73 * pkin(4);
t115 = cos(qJ(5));
t70 = sin(qJ(5));
t52 = t115 * t71 + t70 * t73;
t33 = t119 * t52;
t112 = t70 * t71;
t91 = t115 * t73;
t51 = -t91 + t112;
t114 = t51 * t33;
t87 = t115 * qJD(5);
t32 = -qJD(4) * t91 + t119 * t112 - t73 * t87;
t113 = t52 * t32;
t111 = t71 * t42;
t110 = t73 * t42;
t43 = t102 * qJD(3);
t100 = t71 * qJD(4);
t96 = pkin(4) * t100;
t35 = t43 + t96;
t44 = -pkin(3) - t88;
t41 = t44 - t116;
t107 = t41 * t33 + t35 * t51;
t106 = -t41 * t32 + t35 * t52;
t65 = t73 * qJD(4);
t105 = t43 * t71 + t44 * t65;
t64 = -pkin(3) - t116;
t104 = t64 * t33 + t51 * t96;
t103 = -t64 * t32 + t52 * t96;
t101 = qJD(5) * t70;
t99 = pkin(7) + t102;
t98 = pkin(3) * t100;
t97 = pkin(3) * t65;
t95 = pkin(4) * t101;
t94 = t70 * t118;
t93 = t71 * t65;
t66 = t73 * pkin(8);
t34 = t73 * t99 + t66;
t90 = -pkin(8) - t99;
t81 = t90 * t71;
t77 = t115 * t81;
t12 = -t70 * t34 + t77;
t78 = t70 * t81;
t13 = t115 * t34 + t78;
t80 = qJD(4) * t90;
t75 = -t71 * t80 - t110;
t76 = t73 * t80 - t111;
t3 = -qJD(5) * t77 + t34 * t101 + t115 * t75 - t70 * t76;
t4 = -qJD(5) * t78 + t115 * t76 - t34 * t87 + t70 * t75;
t92 = t12 * t32 - t13 * t33 + t3 * t51 - t4 * t52;
t89 = t44 * t100 - t43 * t73;
t55 = t73 * pkin(7) + t66;
t83 = t118 * t115;
t79 = t71 * t83;
t14 = t55 * t101 + t119 * t79 + t94 * t65;
t85 = t71 * t94;
t37 = t115 * t55 - t85;
t15 = -t37 * qJD(5) + (-t73 * t83 + t85) * qJD(4);
t36 = -t70 * t55 - t79;
t86 = t14 * t51 - t15 * t52 + t36 * t32 - t37 * t33;
t84 = pkin(4) * t87;
t82 = t99 * qJD(4);
t57 = -0.2e1 * t93;
t56 = 0.2e1 * t93;
t50 = 0.2e1 * (-t67 + t68) * qJD(4);
t24 = -0.2e1 * t113;
t23 = 0.2e1 * t114;
t9 = 0.2e1 * t51 * t32 - 0.2e1 * t52 * t33;
t6 = (t115 * t32 - t33 * t70 + (-t115 * t51 + t52 * t70) * qJD(5)) * pkin(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t43, -0.2e1 * t42, 0, 0.2e1 * t102 * t42 - 0.2e1 * t88 * t43, t56, t50, 0, t57, 0, 0, 0.2e1 * t89, 0.2e1 * t105, -0.2e1 * t22, -0.2e1 * t99 * t22 + 0.2e1 * t44 * t43, t24, t9, 0, t23, 0, 0, 0.2e1 * t107, 0.2e1 * t106, 0.2e1 * t92, 0.2e1 * t12 * t4 - 0.2e1 * t13 * t3 + 0.2e1 * t41 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t33 - t13 * t32 - t3 * t52 - t4 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113 + 0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, 0, 0, t56, t50, 0, t57, 0, 0, t89 - t98, -t97 + t105, -t22, -t43 * pkin(3) - pkin(7) * t22, t24, t9, 0, t23, 0, 0, t104 + t107, t103 + t106, t86 + t92, t12 * t15 - t13 * t14 - t3 * t37 + t35 * t64 + t4 * t36 + t41 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t14 - t51 * t15 - t32 * t37 - t33 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t50, 0, t57, 0, 0, -0.2e1 * t98, -0.2e1 * t97, 0, 0, t24, t9, 0, t23, 0, 0, 0.2e1 * t104, 0.2e1 * t103, 0.2e1 * t86, -0.2e1 * t37 * t14 + 0.2e1 * t36 * t15 + 0.2e1 * t64 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t100, 0, -t73 * t82 - t111, t71 * t82 - t110, 0, 0, 0, 0, -t32, 0, -t33, 0, t4, t3, t6, (t115 * t4 - t3 * t70 + (t115 * t13 - t12 * t70) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t65, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, (-t115 * t33 - t32 * t70 + (t115 * t52 + t51 * t70) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t100, 0, -pkin(7) * t65, pkin(7) * t100, 0, 0, 0, 0, -t32, 0, -t33, 0, t15, t14, t6, (t115 * t15 - t14 * t70 + (t115 * t37 - t36 * t70) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t33, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t33, 0, t15, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
