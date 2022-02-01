% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:15
% EndTime: 2022-01-20 11:43:17
% DurationCPUTime: 0.50s
% Computational Cost: add. (942->109), mult. (2199->181), div. (0->0), fcn. (1935->8), ass. (0->96)
t99 = sin(pkin(9));
t133 = pkin(3) * t99;
t100 = cos(pkin(9));
t102 = sin(qJ(3));
t116 = t102 * qJD(3);
t105 = cos(qJ(3));
t97 = t105 * qJD(3);
t70 = t100 * t97 - t99 * t116;
t132 = t70 * pkin(8);
t74 = -t100 * t105 + t99 * t102;
t131 = t74 * pkin(4);
t75 = t100 * t102 + t99 * t105;
t130 = t75 * pkin(8);
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t41 = -t101 * t74 + t104 * t75;
t69 = t75 * qJD(3);
t16 = t41 * qJD(5) + t101 * t70 + t104 * t69;
t40 = t101 * t75 + t104 * t74;
t94 = pkin(3) * t116;
t53 = t69 * pkin(4) + t94;
t103 = sin(qJ(2));
t118 = pkin(1) * qJD(2);
t95 = t103 * t118;
t45 = t53 + t95;
t106 = cos(qJ(2));
t126 = t106 * pkin(1);
t93 = -t105 * pkin(3) - pkin(2);
t83 = t93 - t126;
t52 = t83 + t131;
t129 = t52 * t16 + t45 * t40;
t15 = -t40 * qJD(5) - t101 * t69 + t104 * t70;
t128 = t52 * t15 + t45 * t41;
t56 = t93 + t131;
t127 = t56 * t15 + t53 * t41;
t125 = -qJ(4) - pkin(7);
t124 = t56 * t16 + t53 * t40;
t80 = t95 + t94;
t123 = t83 * t69 + t80 * t74;
t122 = t83 * t70 + t80 * t75;
t121 = t93 * t69 + t74 * t94;
t120 = t93 * t70 + t75 * t94;
t91 = t103 * pkin(1) + pkin(7);
t117 = -qJ(4) - t91;
t72 = t117 * t102;
t98 = t105 * qJ(4);
t73 = t105 * t91 + t98;
t39 = t100 * t73 + t99 * t72;
t84 = t125 * t102;
t85 = t105 * pkin(7) + t98;
t47 = t100 * t85 + t99 * t84;
t92 = -pkin(2) - t126;
t119 = t102 * t95 + t92 * t97;
t115 = pkin(2) * t116;
t114 = pkin(2) * t97;
t113 = t106 * t118;
t108 = t105 * t113;
t109 = qJD(3) * t117;
t96 = t105 * qJD(4);
t43 = t102 * t109 + t108 + t96;
t44 = (-qJD(4) - t113) * t102 + t105 * t109;
t17 = t100 * t44 - t99 * t43;
t112 = qJD(3) * t125;
t67 = t102 * t112 + t96;
t68 = -t102 * qJD(4) + t105 * t112;
t34 = t100 * t68 - t99 * t67;
t38 = t100 * t72 - t99 * t73;
t46 = t100 * t84 - t99 * t85;
t18 = t100 * t43 + t99 * t44;
t111 = -t17 * t75 - t18 * t74 - t38 * t70 - t39 * t69;
t35 = t100 * t67 + t99 * t68;
t110 = -t34 * t75 - t35 * t74 - t46 * t70 - t47 * t69;
t107 = -t105 * t95 + t92 * t116;
t90 = t100 * pkin(3) + pkin(4);
t87 = 0.2e1 * t102 * t97;
t79 = 0.2e1 * (-t102 ^ 2 + t105 ^ 2) * qJD(3);
t71 = t74 * pkin(8);
t66 = t69 * pkin(8);
t59 = (-t101 * t90 - t104 * t133) * qJD(5);
t58 = (t101 * t133 - t104 * t90) * qJD(5);
t37 = -t71 + t47;
t36 = t46 - t130;
t33 = (-t100 * t70 - t69 * t99) * pkin(3);
t30 = -t71 + t39;
t29 = t38 - t130;
t22 = t35 - t66;
t21 = t34 - t132;
t12 = t18 - t66;
t11 = t17 - t132;
t6 = 0.2e1 * t41 * t15;
t5 = -t101 * t22 + t104 * t21 + (-t101 * t36 - t104 * t37) * qJD(5);
t4 = -t101 * t21 - t104 * t22 + (t101 * t37 - t104 * t36) * qJD(5);
t3 = -0.2e1 * t15 * t40 - 0.2e1 * t41 * t16;
t2 = -t101 * t12 + t104 * t11 + (-t101 * t29 - t104 * t30) * qJD(5);
t1 = -t101 * t11 - t104 * t12 + (t101 * t30 - t104 * t29) * qJD(5);
t7 = [0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t113, t87, t79, 0, 0, 0, 0.2e1 * t107, 0.2e1 * t119, 0.2e1 * t123, 0.2e1 * t122, 0.2e1 * t111, 0.2e1 * t38 * t17 + 0.2e1 * t39 * t18 + 0.2e1 * t83 * t80, t6, t3, 0, 0, 0, 0.2e1 * t129, 0.2e1 * t128; 0, 0, 0, 0, -t95, -t113, t87, t79, 0, 0, 0, t107 - t115, -t114 + t119, t121 + t123, t120 + t122, t110 + t111, t17 * t46 + t18 * t47 + t38 * t34 + t39 * t35 + t80 * t93 + t83 * t94, t6, t3, 0, 0, 0, t124 + t129, t127 + t128; 0, 0, 0, 0, 0, 0, t87, t79, 0, 0, 0, -0.2e1 * t115, -0.2e1 * t114, 0.2e1 * t121, 0.2e1 * t120, 0.2e1 * t110, 0.2e1 * t46 * t34 + 0.2e1 * t47 * t35 + 0.2e1 * t93 * t94, t6, t3, 0, 0, 0, 0.2e1 * t124, 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, t97, -t116, 0, -t102 * t113 - t91 * t97, t91 * t116 - t108, t17, -t18, t33, (t100 * t17 + t18 * t99) * pkin(3), 0, 0, t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t97, -t116, 0, -pkin(7) * t97, pkin(7) * t116, t34, -t35, t33, (t100 * t34 + t35 * t99) * pkin(3), 0, 0, t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t59, 0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t70, 0, t80, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t70, 0, t94, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
