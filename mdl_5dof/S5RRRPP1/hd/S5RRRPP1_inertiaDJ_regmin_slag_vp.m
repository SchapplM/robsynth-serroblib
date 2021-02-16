% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:48
% EndTime: 2021-01-15 22:14:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (859->96), mult. (1947->154), div. (0->0), fcn. (1533->6), ass. (0->80)
t105 = cos(pkin(8));
t80 = cos(qJ(3));
t74 = t80 * qJD(4);
t78 = sin(qJ(3));
t79 = sin(qJ(2));
t100 = t79 * pkin(1) + pkin(7);
t89 = -qJ(4) - t100;
t85 = qJD(3) * t89;
t106 = pkin(1) * qJD(2);
t81 = cos(qJ(2));
t101 = t81 * t106;
t90 = t80 * t101;
t38 = t78 * t85 + t74 + t90;
t77 = sin(pkin(8));
t82 = (-qJD(4) - t101) * t78 + t80 * t85;
t12 = -t105 * t82 + t77 * t38;
t13 = t105 * t38 + t77 * t82;
t76 = t80 * qJ(4);
t54 = t80 * t100 + t76;
t86 = t89 * t78;
t34 = -t105 * t86 + t77 * t54;
t35 = t105 * t54 + t77 * t86;
t95 = t105 * t78;
t56 = t77 * t80 + t95;
t52 = t56 * qJD(3);
t104 = t78 * qJD(3);
t94 = t105 * t80;
t53 = qJD(3) * t94 - t77 * t104;
t115 = t77 * t78;
t55 = -t94 + t115;
t126 = t12 * t56 - t13 * t55 + t34 * t53 - t35 * t52;
t130 = 0.2e1 * t126;
t114 = -qJ(4) - pkin(7);
t93 = qJD(3) * t114;
t51 = t78 * t93 + t74;
t83 = -t78 * qJD(4) + t80 * t93;
t31 = -t105 * t83 + t77 * t51;
t32 = t105 * t51 + t77 * t83;
t63 = t80 * pkin(7) + t76;
t39 = -t114 * t95 + t77 * t63;
t40 = t105 * t63 + t114 * t115;
t127 = t31 * t56 - t32 * t55 + t39 * t53 - t40 * t52;
t129 = 0.2e1 * t127;
t128 = t127 + t126;
t125 = 2 * qJD(5);
t124 = t81 * pkin(1);
t72 = pkin(3) * t104;
t16 = t52 * pkin(4) - t53 * qJ(5) - t56 * qJD(5) + t72;
t73 = t79 * t106;
t14 = t16 + t73;
t71 = -t80 * pkin(3) - pkin(2);
t36 = t55 * pkin(4) - t56 * qJ(5) + t71;
t33 = t36 - t124;
t122 = t14 * t55 + t33 * t52;
t121 = -t14 * t56 - t33 * t53;
t120 = t16 * t55 + t36 * t52;
t113 = -t16 * t56 - t36 * t53;
t59 = t73 + t72;
t62 = t71 - t124;
t111 = t62 * t52 + t59 * t55;
t110 = t62 * t53 + t59 * t56;
t109 = t71 * t52 + t55 * t72;
t108 = t71 * t53 + t56 * t72;
t70 = -pkin(2) - t124;
t75 = t80 * qJD(3);
t107 = t70 * t75 + t78 * t73;
t103 = pkin(2) * t104;
t102 = pkin(2) * t75;
t99 = t34 * t12 + t35 * t13;
t96 = t39 * t31 + t40 * t32;
t88 = t100 * qJD(3);
t87 = t70 * t104 - t80 * t73;
t84 = t12 * t39 + t13 * t40 + t34 * t31 + t35 * t32;
t69 = -t105 * pkin(3) - pkin(4);
t67 = t77 * pkin(3) + qJ(5);
t65 = 0.2e1 * t78 * t75;
t58 = 0.2e1 * (-t78 ^ 2 + t80 ^ 2) * qJD(3);
t30 = (-t105 * t53 - t52 * t77) * pkin(3);
t15 = -qJD(5) * t55 - t67 * t52 + t69 * t53;
t1 = [0, 0, 0, 0, -0.2e1 * t73, -0.2e1 * t101, t65, t58, 0, 0, 0, 0.2e1 * t87, 0.2e1 * t107, 0.2e1 * t111, 0.2e1 * t110, t130, 0.2e1 * t62 * t59 + 0.2e1 * t99, 0.2e1 * t122, t130, 0.2e1 * t121, 0.2e1 * t33 * t14 + 0.2e1 * t99; 0, 0, 0, 0, -t73, -t101, t65, t58, 0, 0, 0, t87 - t103, -t102 + t107, t109 + t111, t108 + t110, t128, t59 * t71 + t62 * t72 + t84, t120 + t122, t128, t113 + t121, t14 * t36 + t33 * t16 + t84; 0, 0, 0, 0, 0, 0, t65, t58, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102, 0.2e1 * t109, 0.2e1 * t108, t129, 0.2e1 * t71 * t72 + 0.2e1 * t96, 0.2e1 * t120, t129, 0.2e1 * t113, 0.2e1 * t36 * t16 + 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, t75, -t104, 0, -t78 * t101 - t80 * t88, t78 * t88 - t90, -t12, -t13, t30, (-t105 * t12 + t13 * t77) * pkin(3), -t12, t15, t13, t35 * qJD(5) + t12 * t69 + t13 * t67; 0, 0, 0, 0, 0, 0, 0, 0, t75, -t104, 0, -pkin(7) * t75, pkin(7) * t104, -t31, -t32, t30, (-t105 * t31 + t32 * t77) * pkin(3), -t31, t15, t32, t40 * qJD(5) + t31 * t69 + t32 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t67 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t53, 0, t59, t52, 0, -t53, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t53, 0, t72, t52, 0, -t53, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
