% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:04
% EndTime: 2019-03-09 03:56:07
% DurationCPUTime: 1.11s
% Computational Cost: add. (1827->132), mult. (3745->232), div. (0->0), fcn. (3810->8), ass. (0->93)
t80 = sin(qJ(6));
t117 = qJD(6) * t80;
t128 = cos(qJ(5));
t78 = sin(pkin(10));
t79 = cos(pkin(10));
t82 = sin(qJ(3));
t84 = cos(qJ(3));
t59 = -t78 * t82 + t79 * t84;
t60 = -t78 * t84 - t79 * t82;
t81 = sin(qJ(5));
t33 = t128 * t59 + t81 * t60;
t115 = t84 * qJD(3);
t116 = t82 * qJD(3);
t54 = -t79 * t115 + t78 * t116;
t55 = -t78 * t115 - t79 * t116;
t20 = t33 * qJD(5) - t128 * t54 + t81 * t55;
t83 = cos(qJ(6));
t122 = t83 * t20;
t91 = t128 * t60 - t81 * t59;
t93 = -t117 * t91 - t122;
t18 = t91 * qJD(5) + t128 * t55 + t81 * t54;
t123 = t83 * t18;
t12 = -t33 * t117 + t123;
t125 = t80 * t20;
t74 = qJD(6) * t83;
t10 = -t74 * t91 + t125;
t127 = t33 * t18;
t136 = t91 * t20;
t95 = t80 * t18 + t33 * t74;
t30 = t33 ^ 2;
t133 = (t54 * t78 - t55 * t79) * pkin(3);
t109 = t79 * pkin(3) + pkin(4);
t130 = pkin(3) * t78;
t132 = -t128 * t109 + t81 * t130;
t85 = -pkin(1) - pkin(7);
t118 = qJ(4) - t85;
t62 = t118 * t82;
t63 = t118 * t84;
t35 = t78 * t62 - t79 * t63;
t28 = -t59 * pkin(8) + t35;
t36 = -t79 * t62 - t78 * t63;
t29 = t60 * pkin(8) + t36;
t15 = -t128 * t28 + t81 * t29;
t77 = t83 ^ 2;
t120 = t80 ^ 2 - t77;
t104 = t120 * qJD(6);
t131 = 2 * qJD(2);
t49 = -t84 * qJD(4) + t118 * t116;
t50 = -qJD(3) * t63 - t82 * qJD(4);
t27 = t78 * t49 + t79 * t50;
t103 = t54 * pkin(8) + t27;
t16 = t128 * t29 + t81 * t28;
t26 = t79 * t49 - t78 * t50;
t90 = t55 * pkin(8) - t26;
t5 = t16 * qJD(5) + t81 * t103 + t128 * t90;
t3 = t5 * t80;
t129 = t15 * t74 + t3;
t126 = t33 * t83;
t89 = t81 * t109 + t128 * t130;
t45 = t89 * qJD(5);
t51 = -pkin(5) + t132;
t121 = t45 * t80 + t51 * t74;
t119 = t82 * pkin(3) + qJ(2);
t66 = pkin(3) * t115 + qJD(2);
t114 = qJ(2) * qJD(3);
t112 = pkin(5) * t117;
t111 = pkin(5) * t74;
t110 = t80 * t74;
t107 = -0.4e1 * t80 * t126;
t106 = t51 * t117 - t45 * t83;
t43 = -t60 * pkin(4) + t119;
t37 = -t54 * pkin(4) + t66;
t17 = -pkin(5) * t91 - t33 * pkin(9) + t43;
t102 = t83 * t16 + t80 * t17;
t101 = t80 * t16 - t83 * t17;
t100 = t91 ^ 2 + t30;
t52 = pkin(9) + t89;
t99 = -t33 * t51 - t52 * t91;
t98 = t54 * t60 + t59 * t55;
t88 = -0.2e1 * t127 + 0.2e1 * t136;
t44 = t132 * qJD(5);
t87 = t18 * t51 - t20 * t52 + t33 * t45 - t44 * t91;
t86 = t26 * t59 - t27 * t60 + t35 * t55 - t36 * t54;
t65 = 0.2e1 * t110;
t61 = -0.2e1 * t104;
t13 = t15 * t117;
t8 = -t33 * t104 + t80 * t123;
t7 = t20 * pkin(5) - t18 * pkin(9) + t37;
t6 = qJD(6) * t107 - t120 * t18;
t4 = t15 * qJD(5) - t128 * t103 + t81 * t90;
t2 = -t102 * qJD(6) + t80 * t4 + t83 * t7;
t1 = t101 * qJD(6) + t83 * t4 - t80 * t7;
t9 = [0, 0, 0, 0, t131, qJ(2) * t131, -0.2e1 * t82 * t115, 0.2e1 * (t82 ^ 2 - t84 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t82 + 0.2e1 * t84 * t114, 0.2e1 * qJD(2) * t84 - 0.2e1 * t82 * t114, -0.2e1 * t86, 0.2e1 * t119 * t66 + 0.2e1 * t35 * t26 + 0.2e1 * t36 * t27, 0.2e1 * t127, 0.2e1 * t18 * t91 - 0.2e1 * t33 * t20, 0, 0, 0, 0.2e1 * t43 * t20 - 0.2e1 * t37 * t91, 0.2e1 * t43 * t18 + 0.2e1 * t37 * t33, -0.2e1 * t30 * t110 + 0.2e1 * t77 * t127, 0.2e1 * t30 * t104 + t18 * t107, -0.2e1 * t12 * t91 + 0.2e1 * t33 * t122, -0.2e1 * t33 * t125 + 0.2e1 * t91 * t95, -0.2e1 * t136, -0.2e1 * t101 * t20 + 0.2e1 * t95 * t15 - 0.2e1 * t2 * t91 + 0.2e1 * t33 * t3, -0.2e1 * t1 * t91 - 0.2e1 * t102 * t20 + 0.2e1 * t12 * t15 + 0.2e1 * t5 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t98, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t74 + t88 * t80, t100 * t117 + t88 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, 0, -t85 * t116, -t85 * t115, t133 (t26 * t79 + t27 * t78) * pkin(3), 0, 0, t18, -t20, 0, -t5, t4, t8, t6, t10, -t93, 0, t13 + (-t99 * qJD(6) - t5) * t83 + t87 * t80, t99 * t117 + t87 * t83 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t115, 0, -t133, 0, 0, 0, 0, 0, t18, -t20, 0, 0, 0, 0, 0, t12, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t45, 0.2e1 * t44, t65, t61, 0, 0, 0, 0.2e1 * t106, 0.2e1 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, t20, t18, 0, 0, 0, 0, 0, -t93, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t20, 0, -t5, t4, t8, t6, t10, -t93, 0, t13 + (-pkin(5) * t18 - pkin(9) * t20) * t80 + (-t5 + (-pkin(5) * t33 + pkin(9) * t91) * qJD(6)) * t83, -pkin(5) * t12 + t93 * pkin(9) + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t20, 0, 0, 0, 0, 0, t12, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44, t65, t61, 0, 0, 0, t106 - t112, -t111 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t61, 0, 0, 0, -0.2e1 * t112, -0.2e1 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t95, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t117, 0, t80 * t44 - t52 * t74, t52 * t117 + t83 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t117, 0, -pkin(9) * t74, pkin(9) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
