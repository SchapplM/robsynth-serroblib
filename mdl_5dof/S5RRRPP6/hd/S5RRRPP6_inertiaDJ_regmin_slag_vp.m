% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP6
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
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:03
% DurationCPUTime: 1.02s
% Computational Cost: add. (1358->168), mult. (3391->315), div. (0->0), fcn. (2684->6), ass. (0->96)
t81 = cos(qJ(3));
t114 = qJD(3) * t81;
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t110 = t82 * qJD(2);
t79 = sin(qJ(3));
t99 = t79 * t110;
t130 = t80 * t114 + t99;
t129 = -0.4e1 * t80;
t88 = -t82 * pkin(2) - t80 * pkin(7);
t56 = -pkin(1) + t88;
t124 = t81 * t82;
t64 = pkin(6) * t124;
t120 = t79 * t56 + t64;
t115 = qJD(3) * t79;
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t41 = t78 * t114 - t77 * t115;
t74 = t80 ^ 2;
t91 = (-t82 ^ 2 + t74) * qJD(2);
t75 = t81 ^ 2;
t119 = t79 ^ 2 - t75;
t92 = t119 * qJD(3);
t128 = 2 * qJD(5);
t87 = pkin(2) * t80 - pkin(7) * t82;
t50 = t87 * qJD(2);
t122 = -t56 * t114 - t79 * t50;
t125 = t80 * t81;
t12 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t125 + (-qJD(4) * t80 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t82) * t79 - t122;
t111 = t81 * qJD(4);
t116 = qJ(4) * t81;
t117 = qJ(4) * t80;
t112 = t80 * qJD(2);
t100 = t79 * t112;
t121 = pkin(6) * t100 + t81 * t50;
t8 = -t80 * t111 + (pkin(3) * t80 - t82 * t116) * qJD(2) + (-t64 + (-t56 + t117) * t79) * qJD(3) + t121;
t4 = t78 * t12 + t77 * t8;
t127 = pkin(6) * t79;
t126 = t79 * t80;
t123 = -qJ(4) - pkin(7);
t48 = t81 * t56;
t30 = -t80 * t116 + t48 + (-pkin(3) - t127) * t82;
t32 = -t79 * t117 + t120;
t16 = t77 * t30 + t78 * t32;
t51 = pkin(3) * t126 + t80 * pkin(6);
t113 = qJD(3) * t82;
t109 = t82 * qJD(5);
t108 = qJ(5) * t112 + t4;
t107 = -0.2e1 * pkin(1) * qJD(2);
t106 = -0.2e1 * pkin(2) * qJD(3);
t70 = pkin(6) * t110;
t35 = t130 * pkin(3) + t70;
t71 = pkin(3) * t115;
t105 = t79 * t113;
t103 = t81 * t113;
t98 = t79 * t114;
t97 = t80 * t110;
t96 = t81 * t110;
t69 = -t81 * pkin(3) - pkin(2);
t3 = -t77 * t12 + t78 * t8;
t95 = t123 * t79;
t93 = qJD(3) * t123;
t39 = t79 * t93 + t111;
t84 = -t79 * qJD(4) + t81 * t93;
t26 = t77 * t39 - t78 * t84;
t27 = t78 * t39 + t77 * t84;
t57 = t123 * t81;
t33 = -t77 * t57 - t78 * t95;
t34 = -t78 * t57 + t77 * t95;
t94 = t33 * t26 + t34 * t27;
t90 = 0.2e1 * t97;
t89 = t79 * t96;
t15 = t78 * t30 - t77 * t32;
t46 = t77 * t81 + t78 * t79;
t40 = t46 * qJD(3);
t24 = -t41 * t80 - t77 * t96 - t78 * t99;
t25 = t80 * t40 + t77 * t99 - t78 * t96;
t37 = t46 * t80;
t38 = t78 * t125 - t77 * t126;
t86 = t34 * t24 - t33 * t25 + t26 * t38 - t27 * t37;
t85 = t81 * t112 + t105;
t45 = t77 * t79 - t78 * t81;
t83 = 0.2e1 * t26 * t46 - 0.2e1 * t27 * t45 + 0.2e1 * t33 * t41 - 0.2e1 * t34 * t40;
t67 = -t78 * pkin(3) - pkin(4);
t65 = t77 * pkin(3) + qJ(5);
t28 = t45 * pkin(4) - t46 * qJ(5) + t69;
t22 = -t120 * qJD(3) + t121;
t21 = t85 * pkin(6) + t122;
t20 = t37 * pkin(4) - t38 * qJ(5) + t51;
t18 = t40 * pkin(4) - t41 * qJ(5) - t46 * qJD(5) + t71;
t11 = t82 * pkin(4) - t15;
t10 = -t82 * qJ(5) + t16;
t5 = -t24 * pkin(4) + t25 * qJ(5) - t38 * qJD(5) + t35;
t2 = -pkin(4) * t112 - t3;
t1 = t108 - t109;
t6 = [0, 0, 0, t90, -0.2e1 * t91, 0, 0, 0, t80 * t107, t82 * t107, -0.2e1 * t74 * t98 + 0.2e1 * t75 * t97, t89 * t129 + 0.2e1 * t74 * t92, 0.2e1 * t80 * t105 + 0.2e1 * t81 * t91, 0.2e1 * t80 * t103 - 0.2e1 * t79 * t91, -0.2e1 * t97, 0.2e1 * t48 * t112 - 0.2e1 * t22 * t82 + 0.2e1 * (t74 * t114 + t79 * t97) * pkin(6), -0.2e1 * t21 * t82 - 0.2e1 * t120 * t112 + 0.2e1 * (-t74 * t115 + t81 * t90) * pkin(6), 0.2e1 * t15 * t25 + 0.2e1 * t16 * t24 - 0.2e1 * t3 * t38 - 0.2e1 * t4 * t37, 0.2e1 * t15 * t3 + 0.2e1 * t16 * t4 + 0.2e1 * t51 * t35, -0.2e1 * t11 * t112 + 0.2e1 * t2 * t82 - 0.2e1 * t20 * t24 + 0.2e1 * t5 * t37, -0.2e1 * t1 * t37 + 0.2e1 * t10 * t24 - 0.2e1 * t11 * t25 + 0.2e1 * t2 * t38, -0.2e1 * t1 * t82 + 0.2e1 * t10 * t112 + 0.2e1 * t20 * t25 - 0.2e1 * t5 * t38, 0.2e1 * t10 * t1 + 0.2e1 * t11 * t2 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, t110, -t112, 0, -t70, pkin(6) * t112, -t80 * t92 + t89, -t119 * t110 + t98 * t129, t100 - t103, t85, 0, (pkin(7) * t124 + (-pkin(2) * t81 + t127) * t80) * qJD(3) + (t88 * t79 - t64) * qJD(2), (pkin(6) * t125 + t87 * t79) * qJD(3) + (t82 * t127 + t88 * t81) * qJD(2), -t15 * t41 - t16 * t40 - t3 * t46 - t4 * t45 + t86, -t15 * t26 + t16 * t27 - t3 * t33 + t4 * t34 + t35 * t69 + t51 * t71, -t33 * t112 + t18 * t37 + t20 * t40 - t28 * t24 + t26 * t82 + t5 * t45, -t1 * t45 - t10 * t40 + t11 * t41 + t2 * t46 + t86, t34 * t112 - t18 * t38 - t20 * t41 + t28 * t25 - t27 * t82 - t5 * t46, t1 * t34 + t10 * t27 + t11 * t26 + t20 * t18 + t2 * t33 + t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98, -0.2e1 * t92, 0, 0, 0, t79 * t106, t81 * t106, t83, 0.2e1 * t69 * t71 + 0.2e1 * t94, 0.2e1 * t18 * t45 + 0.2e1 * t28 * t40, t83, -0.2e1 * t18 * t46 - 0.2e1 * t28 * t41, 0.2e1 * t28 * t18 + 0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t115 + t96, -t130, t112, t22, t21, (t24 * t77 + t25 * t78) * pkin(3), (t3 * t78 + t4 * t77) * pkin(3), (pkin(4) - t67) * t112 + t3, -qJD(5) * t37 + t65 * t24 - t67 * t25, t65 * t112 + t108 - 0.2e1 * t109, t10 * qJD(5) + t1 * t65 + t2 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t115, 0, -pkin(7) * t114, pkin(7) * t115, (-t40 * t77 - t41 * t78) * pkin(3), (-t26 * t78 + t27 * t77) * pkin(3), -t26, -qJD(5) * t45 - t65 * t40 + t67 * t41, t27, t34 * qJD(5) + t26 * t67 + t27 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t65 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t24, 0, t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t40, 0, -t41, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
