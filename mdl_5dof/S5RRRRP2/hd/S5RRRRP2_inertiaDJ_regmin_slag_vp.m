% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:28
% DurationCPUTime: 0.55s
% Computational Cost: add. (1081->111), mult. (2454->168), div. (0->0), fcn. (2046->6), ass. (0->90)
t129 = qJD(3) + qJD(4);
t128 = -pkin(7) - pkin(8);
t88 = cos(qJ(2));
t127 = t88 * pkin(1);
t86 = sin(qJ(2));
t76 = t86 * pkin(1) + pkin(7);
t126 = -pkin(8) - t76;
t125 = cos(qJ(4));
t84 = sin(qJ(4));
t85 = sin(qJ(3));
t124 = t84 * t85;
t87 = cos(qJ(3));
t62 = t125 * t85 + t84 * t87;
t41 = t129 * t62;
t110 = t85 * qJD(3);
t80 = pkin(3) * t110;
t35 = t41 * pkin(4) + t80;
t113 = pkin(1) * qJD(2);
t81 = t86 * t113;
t30 = t35 + t81;
t61 = -t125 * t87 + t124;
t79 = -t87 * pkin(3) - pkin(2);
t48 = t61 * pkin(4) + t79;
t47 = t48 - t127;
t123 = t30 * t61 + t47 * t41;
t100 = t125 * qJD(4);
t101 = t125 * qJD(3);
t40 = t129 * t124 + (-t101 - t100) * t87;
t122 = t30 * t62 - t47 * t40;
t121 = t35 * t61 + t48 * t41;
t120 = t35 * t62 - t48 * t40;
t119 = t41 * qJ(5) + t61 * qJD(5);
t63 = t81 + t80;
t66 = t79 - t127;
t118 = t66 * t41 + t63 * t61;
t117 = -t66 * t40 + t63 * t62;
t116 = t79 * t41 + t61 * t80;
t115 = -t79 * t40 + t62 * t80;
t78 = -pkin(2) - t127;
t82 = t87 * qJD(3);
t114 = t78 * t82 + t85 * t81;
t112 = t62 * qJ(5);
t111 = qJD(4) * t84;
t109 = pkin(2) * t110;
t108 = pkin(2) * t82;
t107 = t88 * t113;
t106 = pkin(3) * t111;
t105 = t125 * pkin(3);
t57 = t126 * t85;
t83 = t87 * pkin(8);
t58 = t87 * t76 + t83;
t27 = t125 * t57 - t84 * t58 - t112;
t56 = t61 * qJ(5);
t92 = -t125 * t58 - t84 * t57;
t28 = -t56 - t92;
t102 = qJD(3) * t126;
t98 = t87 * t107;
t89 = t85 * t102 + t98;
t99 = t85 * t107;
t90 = t87 * t102 - t99;
t10 = -t57 * t100 + t58 * t111 - t125 * t89 - t84 * t90;
t3 = t10 + t119;
t11 = t92 * qJD(4) + t125 * t90 - t84 * t89;
t95 = t40 * qJ(5) - t62 * qJD(5);
t4 = t11 + t95;
t104 = t27 * t40 - t28 * t41 + t3 * t61 - t4 * t62;
t67 = t128 * t85;
t68 = t87 * pkin(7) + t83;
t33 = t125 * t67 - t84 * t68 - t112;
t91 = -t125 * t68 - t84 * t67;
t34 = -t56 - t91;
t93 = t128 * t101;
t97 = qJD(3) * t84 * t128;
t17 = -t67 * t100 + t68 * t111 - t85 * t93 - t87 * t97;
t8 = t17 + t119;
t18 = t91 * qJD(4) - t85 * t97 + t87 * t93;
t9 = t18 + t95;
t103 = t33 * t40 - t34 * t41 + t8 * t61 - t9 * t62;
t96 = pkin(3) * t100;
t94 = t78 * t110 - t87 * t81;
t77 = t105 + pkin(4);
t74 = -0.2e1 * t96;
t73 = -0.2e1 * t106;
t70 = 0.2e1 * t85 * t82;
t60 = 0.2e1 * (-t85 ^ 2 + t87 ^ 2) * qJD(3);
t39 = t40 * pkin(4);
t29 = -0.2e1 * t62 * t40;
t12 = 0.2e1 * t40 * t61 - 0.2e1 * t62 * t41;
t7 = t77 * t40 + (-t41 * t84 + (-t125 * t61 + t62 * t84) * qJD(4)) * pkin(3);
t1 = [0, 0, 0, 0, -0.2e1 * t81, -0.2e1 * t107, t70, t60, 0, 0, 0, 0.2e1 * t94, 0.2e1 * t114, t29, t12, 0, 0, 0, 0.2e1 * t118, 0.2e1 * t117, 0.2e1 * t123, 0.2e1 * t122, 0.2e1 * t104, 0.2e1 * t27 * t4 - 0.2e1 * t28 * t3 + 0.2e1 * t47 * t30; 0, 0, 0, 0, -t81, -t107, t70, t60, 0, 0, 0, t94 - t109, -t108 + t114, t29, t12, 0, 0, 0, t116 + t118, t115 + t117, t121 + t123, t120 + t122, t103 + t104, t27 * t9 - t28 * t8 - t3 * t34 + t30 * t48 + t4 * t33 + t47 * t35; 0, 0, 0, 0, 0, 0, t70, t60, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t108, t29, t12, 0, 0, 0, 0.2e1 * t116, 0.2e1 * t115, 0.2e1 * t121, 0.2e1 * t120, 0.2e1 * t103, 0.2e1 * t33 * t9 - 0.2e1 * t34 * t8 + 0.2e1 * t48 * t35; 0, 0, 0, 0, 0, 0, 0, 0, t82, -t110, 0, -t76 * t82 - t99, t76 * t110 - t98, 0, 0, -t40, -t41, 0, t11, t10, t4, t3, t7, t4 * t77 + (-t3 * t84 + (t125 * t28 - t27 * t84) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, t82, -t110, 0, -pkin(7) * t82, pkin(7) * t110, 0, 0, -t40, -t41, 0, t18, t17, t9, t8, t7, t9 * t77 + (-t8 * t84 + (t125 * t34 - t33 * t84) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t74, t73, t74, 0, 0.2e1 * (t105 - t77) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41, 0, t11, t10, t4, t3, t39, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41, 0, t18, t17, t9, t8, t39, t9 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t96, -t106, -t96, 0, -pkin(4) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
