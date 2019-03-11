% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:01
% EndTime: 2019-03-09 05:21:04
% DurationCPUTime: 1.15s
% Computational Cost: add. (2228->138), mult. (4548->242), div. (0->0), fcn. (4388->8), ass. (0->100)
t113 = cos(pkin(10));
t68 = sin(qJ(4));
t69 = sin(qJ(3));
t71 = cos(qJ(4));
t72 = cos(qJ(3));
t116 = t68 * t72 + t71 * t69;
t50 = -t68 * t69 + t71 * t72;
t66 = sin(pkin(10));
t126 = t113 * t116 + t66 * t50;
t127 = qJD(3) + qJD(4);
t36 = t127 * t116;
t37 = t127 * t50;
t81 = t113 * t37 - t66 * t36;
t141 = t126 * t81;
t21 = t113 * t36 + t37 * t66;
t32 = -t113 * t50 + t66 * t116;
t146 = t21 * t32;
t139 = t141 + t146;
t29 = t32 ^ 2;
t150 = (t113 * t21 - t66 * t81) * pkin(4);
t119 = t66 * t68;
t60 = pkin(3) * t71 + pkin(4);
t45 = -pkin(3) * t119 + t113 * t60;
t98 = t113 * t68;
t46 = pkin(3) * t98 + t66 * t60;
t114 = pkin(3) * qJD(4);
t43 = (t66 * t71 + t98) * t114;
t44 = (t113 * t71 - t119) * t114;
t88 = -t126 * t44 - t32 * t43;
t149 = t21 * t45 - t46 * t81 + t88;
t73 = -pkin(1) - pkin(7);
t124 = pkin(8) - t73;
t102 = t124 * t72;
t128 = t124 * t69;
t80 = t68 * t102 + t128 * t71;
t27 = -t116 * qJ(5) - t80;
t79 = -t71 * t102 + t128 * t68;
t74 = -t50 * qJ(5) + t79;
t17 = -t113 * t74 + t27 * t66;
t18 = t113 * t27 + t66 * t74;
t24 = t127 * t79;
t10 = -t37 * qJ(5) - t116 * qJD(5) + t24;
t25 = t127 * t80;
t143 = t36 * qJ(5) - t50 * qJD(5) + t25;
t4 = t10 * t66 - t113 * t143;
t5 = t113 * t10 + t143 * t66;
t148 = t126 * t5 + t17 * t21 + t18 * t81 + t32 * t4;
t70 = cos(qJ(6));
t62 = qJD(6) * t70;
t67 = sin(qJ(6));
t13 = t21 * t67 + t32 * t62;
t112 = qJD(6) * t67;
t147 = t32 * t112 - t21 * t70;
t12 = t126 * t62 + t67 * t81;
t137 = t112 * t126 - t70 * t81;
t65 = t70 ^ 2;
t115 = t67 ^ 2 - t65;
t97 = t115 * qJD(6);
t125 = 2 * qJD(2);
t123 = t17 * t62 + t4 * t67;
t122 = t32 * t67;
t121 = t32 * t70;
t118 = t67 * t70;
t41 = -pkin(5) - t45;
t117 = t41 * t62 + t43 * t67;
t57 = t69 * pkin(3) + qJ(2);
t111 = t69 * qJD(3);
t110 = t72 * qJD(3);
t53 = pkin(3) * t110 + qJD(2);
t109 = qJ(2) * qJD(3);
t108 = t68 * t114;
t107 = t71 * t114;
t56 = -t113 * pkin(4) - pkin(5);
t106 = t56 * t112;
t105 = t56 * t62;
t104 = t67 * t62;
t100 = 0.4e1 * t32 * t118;
t99 = t41 * t112 - t43 * t70;
t28 = pkin(4) * t37 + t53;
t84 = t116 * pkin(4) + t57;
t19 = pkin(5) * t126 + pkin(9) * t32 + t84;
t93 = t18 * t70 + t19 * t67;
t92 = t18 * t67 - t19 * t70;
t55 = pkin(4) * t66 + pkin(9);
t91 = -t21 * t56 - t55 * t81;
t90 = -t126 ^ 2 - t29;
t42 = pkin(9) + t46;
t89 = t126 * t42 + t32 * t41;
t87 = t126 * t55 + t32 * t56;
t77 = -0.2e1 * t139;
t76 = -t21 * t41 - t42 * t81 + t88;
t52 = 0.2e1 * t104;
t49 = -0.2e1 * t97;
t15 = t17 * t112;
t8 = -t118 * t21 + t32 * t97;
t7 = pkin(5) * t81 + pkin(9) * t21 + t28;
t6 = qJD(6) * t100 + t115 * t21;
t2 = -t93 * qJD(6) - t67 * t5 + t70 * t7;
t1 = t92 * qJD(6) - t70 * t5 - t67 * t7;
t3 = [0, 0, 0, 0, t125, qJ(2) * t125, -0.2e1 * t69 * t110, 0.2e1 * (t69 ^ 2 - t72 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t69 + 0.2e1 * t72 * t109, 0.2e1 * qJD(2) * t72 - 0.2e1 * t69 * t109, -0.2e1 * t50 * t36, 0.2e1 * t36 * t116 - 0.2e1 * t50 * t37, 0, 0, 0, 0.2e1 * t53 * t116 + 0.2e1 * t57 * t37, -0.2e1 * t36 * t57 + 0.2e1 * t50 * t53, -0.2e1 * t148, 0.2e1 * t17 * t4 + 0.2e1 * t18 * t5 + 0.2e1 * t84 * t28, -0.2e1 * t29 * t104 + 0.2e1 * t65 * t146, -t100 * t21 + 0.2e1 * t29 * t97, -0.2e1 * t121 * t81 + 0.2e1 * t126 * t147, 0.2e1 * t122 * t81 + 0.2e1 * t126 * t13, 0.2e1 * t141, -0.2e1 * t4 * t122 + 0.2e1 * t126 * t2 - 0.2e1 * t13 * t17 - 0.2e1 * t81 * t92, 0.2e1 * t1 * t126 - 0.2e1 * t4 * t121 + 0.2e1 * t147 * t17 - 0.2e1 * t81 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t148, 0, 0, 0, 0, 0, t90 * t62 + t77 * t67, -t90 * t112 + t77 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t139, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, 0, -t73 * t111, -t73 * t110, 0, 0, -t36, -t37, 0, t25, -t24, t149, t17 * t43 + t18 * t44 - t4 * t45 + t46 * t5, t8, t6, t12, -t137, 0, t15 + (-t89 * qJD(6) - t4) * t70 + t76 * t67, t89 * t112 + t76 * t70 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, 0, 0, 0, 0, 0, -t36, -t37, 0, -t149, 0, 0, 0, 0, 0, t147, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108, -0.2e1 * t107, 0, -0.2e1 * t43 * t45 + 0.2e1 * t44 * t46, t52, t49, 0, 0, 0, 0.2e1 * t99, 0.2e1 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t25, -t24, t150 (-t113 * t4 + t5 * t66) * pkin(4), t8, t6, t12, -t137, 0, t15 + t91 * t67 + (-t87 * qJD(6) - t4) * t70, t87 * t112 + t91 * t70 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, -t150, 0, 0, 0, 0, 0, t147, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t107, 0 (-t113 * t43 + t44 * t66) * pkin(4), t52, t49, 0, 0, 0, t99 + t106, t105 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t49, 0, 0, 0, 0.2e1 * t106, 0.2e1 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, -t137, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t13, t81, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t112, 0, -t42 * t62 - t44 * t67, t42 * t112 - t44 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t112, 0, -t55 * t62, t55 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
