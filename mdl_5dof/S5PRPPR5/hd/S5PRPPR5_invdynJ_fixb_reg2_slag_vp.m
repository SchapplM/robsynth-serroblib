% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:41
% EndTime: 2019-12-31 17:38:43
% DurationCPUTime: 0.78s
% Computational Cost: add. (964->176), mult. (1657->228), div. (0->0), fcn. (1120->8), ass. (0->103)
t66 = sin(pkin(7));
t65 = sin(pkin(8));
t67 = cos(pkin(8));
t70 = sin(qJ(2));
t72 = cos(qJ(2));
t88 = t72 * t65 - t70 * t67;
t18 = t88 * t66;
t68 = cos(pkin(7));
t20 = t88 * t68;
t114 = t72 * qJD(1);
t115 = t70 * qJD(1);
t24 = t65 * t114 - t67 * t115;
t28 = t70 * t65 + t72 * t67;
t141 = -g(1) * t20 - g(2) * t18 - g(3) * t28 + (t65 * qJD(3) - t24) * qJD(2);
t25 = t28 * qJD(1);
t98 = t67 * qJD(3) - t25;
t140 = t98 * qJD(2);
t57 = t72 * qJDD(1);
t104 = qJD(2) * t115 + qJDD(3) - t57;
t137 = pkin(2) + pkin(3);
t16 = -t137 * qJDD(2) + t104;
t108 = qJDD(2) * qJ(3);
t54 = t70 * qJDD(1);
t21 = t108 + t54 + (qJD(3) + t114) * qJD(2);
t5 = t67 * t16 - t65 * t21;
t3 = qJDD(2) * pkin(4) - t5;
t34 = -t65 * qJ(3) - t137 * t67;
t31 = pkin(4) - t34;
t35 = t67 * qJ(3) - t137 * t65;
t32 = -pkin(6) + t35;
t74 = qJD(5) ^ 2;
t139 = qJDD(2) * t31 - t32 * t74 + t141 + t3;
t103 = g(1) * t66 - g(2) * t68;
t94 = g(1) * t68 + g(2) * t66;
t138 = g(3) * t70 + t94 * t72 - t54;
t96 = qJD(3) - t114;
t30 = -t137 * qJD(2) + t96;
t45 = qJD(2) * qJ(3) + t115;
t12 = t65 * t30 + t67 * t45;
t10 = -qJD(2) * pkin(6) + t12;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t7 = t71 * qJD(4) - t69 * t10;
t121 = t7 * qJD(5);
t6 = t65 * t16 + t67 * t21;
t4 = -qJDD(2) * pkin(6) + t6;
t1 = t69 * qJDD(4) + t71 * t4 + t121;
t8 = t69 * qJD(4) + t71 * t10;
t120 = t8 * qJD(5);
t56 = t71 * qJDD(4);
t2 = -t69 * t4 - t120 + t56;
t79 = -(t69 * t8 + t7 * t71) * qJD(5) + t1 * t71 - t2 * t69;
t133 = t66 * t70;
t75 = qJD(2) ^ 2;
t132 = t67 * t75;
t131 = t68 * t70;
t126 = t72 * pkin(2) + t70 * qJ(3);
t63 = t69 ^ 2;
t64 = t71 ^ 2;
t125 = t63 - t64;
t124 = t63 + t64;
t123 = t74 + t75;
t122 = qJ(3) * t72;
t11 = t67 * t30 - t65 * t45;
t9 = qJD(2) * pkin(4) - t11;
t119 = t9 * qJD(2);
t27 = qJD(2) * t28;
t118 = qJDD(2) * pkin(2);
t117 = t27 * qJD(2);
t113 = qJDD(5) * t69;
t112 = qJDD(5) * t71;
t111 = t69 * qJDD(2);
t110 = t71 * qJDD(2);
t109 = qJD(2) * qJD(5);
t107 = t69 * t75 * t71;
t106 = g(1) * t131 + g(2) * t133 - g(3) * t72;
t105 = t72 * pkin(3) + t126;
t101 = t71 * t109;
t100 = t57 + t106;
t97 = qJDD(2) * t124;
t95 = t69 * t101;
t91 = t7 * t69 - t8 * t71;
t90 = t11 * t65 - t12 * t67;
t89 = (-qJD(2) * pkin(2) + t96) * t70 + t45 * t72;
t26 = t88 * qJD(2);
t86 = t26 * qJD(2) + t28 * qJDD(2);
t22 = t104 - t118;
t17 = t28 * t66;
t19 = t28 * t68;
t85 = -g(1) * t19 - g(2) * t17 + g(3) * t88;
t83 = -qJD(5) * t10 + t103;
t82 = t74 * t88 + t86;
t81 = -0.2e1 * t27 * qJD(5) + qJDD(5) * t88;
t80 = t94 * t137 * t70;
t78 = -qJDD(5) * t32 + (-qJD(2) * t31 - t9 - t98) * qJD(5);
t76 = -qJD(5) * qJD(4) + t119 - t4 - t85;
t49 = t68 * t122;
t48 = t66 * t122;
t39 = t72 * qJDD(2) - t75 * t70;
t38 = qJDD(2) * t70 + t75 * t72;
t37 = -t74 * t69 + t112;
t36 = -t74 * t71 - t113;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, t39, -t38, 0, -g(3) + (t70 ^ 2 + t72 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t39, 0, t38, t89 * qJD(2) + t21 * t70 - t22 * t72 - g(3), 0, 0, 0, 0, 0, 0, t86, -qJDD(2) * t88 + t117, 0, -t11 * t26 + t12 * t27 - t5 * t28 - t6 * t88 - g(3), 0, 0, 0, 0, 0, 0, t81 * t69 + t82 * t71, -t82 * t69 + t81 * t71, -t124 * t117 + t88 * t97, t9 * t26 - t91 * t27 + t3 * t28 - t79 * t88 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t100, t138, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -qJDD(3) + t100 + 0.2e1 * t118, 0, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t108 - t138, t21 * qJ(3) + t45 * qJD(3) - t22 * pkin(2) - g(1) * (-pkin(2) * t131 + t49) - g(2) * (-pkin(2) * t133 + t48) - g(3) * t126 - t89 * qJD(1), 0, 0, 0, 0, 0, qJDD(2), -t34 * qJDD(2) + t141 - t5, t35 * qJDD(2) + t140 + t6 + t85, 0, -g(1) * t49 - g(2) * t48 - g(3) * t105 - t90 * qJD(3) + t11 * t24 - t12 * t25 + t5 * t34 + t6 * t35 + t80, t63 * qJDD(2) + 0.2e1 * t95, -0.2e1 * t125 * t109 + 0.2e1 * t69 * t110, t36, t64 * qJDD(2) - 0.2e1 * t95, -t37, 0, t139 * t71 + t78 * t69, -t139 * t69 + t78 * t71, -t124 * t140 - t32 * t97 - t79 - t85, t3 * t31 - t9 * t24 - g(1) * (t20 * pkin(4) - t19 * pkin(6) + t49) - g(2) * (t18 * pkin(4) - t17 * pkin(6) + t48) - g(3) * (t28 * pkin(4) + pkin(6) * t88 + t105) + t80 + t91 * t25 + (t9 * t65 - t91 * t67) * qJD(3) + t79 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t75, -t45 * qJD(2) - t106 + t22, 0, 0, 0, 0, 0, 0, -t67 * qJDD(2) - t65 * t75, t65 * qJDD(2) - t132, 0, t90 * qJD(2) + t5 * t67 + t6 * t65 - t106, 0, 0, 0, 0, 0, 0, (0.2e1 * t69 * t109 - t110) * t67 + (-t123 * t71 - t113) * t65, (0.2e1 * t101 + t111) * t67 + (t123 * t69 - t112) * t65, t124 * t132 - t65 * t97, (t91 * qJD(2) - t3) * t67 + (t79 - t119) * t65 - t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + t103, 0, 0, 0, 0, 0, 0, t37, t36, 0, -t91 * qJD(5) + t1 * t69 + t2 * t71 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t125 * t75, -t111, t107, -t110, qJDD(5), t76 * t69 + t83 * t71 + t120 + t56, t121 + (-qJDD(4) - t83) * t69 + t76 * t71, 0, 0;];
tau_reg = t13;
