% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:57
% EndTime: 2022-01-23 09:20:58
% DurationCPUTime: 0.65s
% Computational Cost: add. (820->115), mult. (1822->164), div. (0->0), fcn. (1454->8), ass. (0->115)
t79 = sin(pkin(9));
t77 = t79 ^ 2;
t81 = cos(pkin(9));
t78 = t81 ^ 2;
t72 = t77 + t78;
t110 = qJD(1) + qJD(3);
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t59 = (t82 ^ 2 - t84 ^ 2) * t77;
t145 = t110 * t59;
t60 = t72 * t82;
t27 = t110 * t60;
t61 = t72 * t84;
t28 = t110 * t61;
t144 = t110 * t72;
t98 = t110 * t81;
t83 = sin(qJ(3));
t95 = cos(pkin(8)) * pkin(1) + pkin(2);
t88 = t83 * t95;
t138 = cos(qJ(3));
t140 = pkin(1) * sin(pkin(8));
t97 = t138 * t140;
t55 = t97 + t88;
t53 = qJ(4) + t55;
t143 = qJ(4) + t53;
t142 = t53 / 0.2e1;
t54 = -t138 * t95 + t83 * t140;
t141 = t54 / 0.2e1;
t139 = qJ(4) / 0.2e1;
t131 = t81 * t82;
t63 = -t81 * pkin(4) - t79 * pkin(7) - pkin(3);
t31 = t54 + t63;
t16 = t53 * t131 - t84 * t31;
t137 = t16 * t81;
t130 = t81 * t84;
t17 = -t53 * t130 - t82 * t31;
t136 = t17 * t81;
t124 = qJ(4) * t81;
t42 = t82 * t124 - t84 * t63;
t135 = t42 * t81;
t43 = -t84 * t124 - t82 * t63;
t134 = t43 * t81;
t133 = t77 * t82;
t132 = t77 * t84;
t129 = t82 * t55;
t128 = t84 * t55;
t13 = (-t54 * t130 + t129) * t81 - t54 * t132;
t58 = t61 * qJD(4);
t127 = t13 * qJD(3) + t58;
t20 = t72 * t54;
t67 = t72 * qJD(4);
t126 = -t20 * qJD(3) + t67;
t12 = (t54 * t131 + t128) * t81 + t54 * t133;
t57 = t60 * qJD(4);
t125 = -t12 * qJD(3) + t57;
t5 = (-pkin(3) + t54) * t55 - t53 * t20;
t123 = t5 * qJD(1);
t6 = -t53 * t133 - t137;
t122 = t6 * qJD(1);
t38 = t53 * t132;
t7 = -t38 + t136;
t121 = t7 * qJD(1);
t120 = t77 * qJ(4);
t119 = qJD(4) * t81;
t118 = qJD(5) * t82;
t117 = qJD(5) * t84;
t116 = t12 * qJD(1);
t115 = t13 * qJD(1);
t18 = t72 * t53;
t114 = t18 * qJD(1);
t113 = t20 * qJD(1);
t112 = t54 * qJD(1);
t111 = t55 * qJD(1);
t52 = t55 * qJD(3);
t109 = t82 * t132;
t108 = t81 * t118;
t107 = t81 * t117;
t106 = t82 * t119;
t105 = t84 * t119;
t104 = t82 * t141;
t103 = t84 * t141;
t101 = t128 / 0.2e1;
t73 = t84 * t120;
t100 = -t38 / 0.2e1 - t73 / 0.2e1;
t96 = qJD(5) * t109;
t94 = t82 * t98;
t93 = t84 * t98;
t1 = t101 + (t104 + t43 / 0.2e1 + t17 / 0.2e1) * t81 + t100;
t24 = -t73 + t134;
t92 = -t1 * qJD(1) - t24 * qJD(3);
t2 = (-t55 / 0.2e1 + (t139 + t142) * t77) * t82 + (t103 + t42 / 0.2e1 + t16 / 0.2e1) * t81;
t23 = -t82 * t120 - t135;
t91 = -t2 * qJD(1) + t23 * qJD(3);
t66 = t72 * qJ(4);
t85 = t97 / 0.2e1 + t88 / 0.2e1;
t8 = t85 + t143 * (-t78 / 0.2e1 - t77 / 0.2e1);
t90 = t8 * qJD(1) - t66 * qJD(3);
t89 = -qJD(5) + t98;
t87 = t89 * t82;
t86 = t89 * t84;
t65 = t79 * t107;
t64 = t79 * t108;
t56 = t59 * qJD(5);
t51 = t54 * qJD(3);
t48 = t110 * t109;
t45 = t79 * t93;
t44 = t79 * t94;
t30 = t79 * t86;
t29 = t79 * t87;
t22 = t107 - t28;
t21 = t108 - t27;
t9 = t78 * t139 + t120 / 0.2e1 + t85 + t72 * t142;
t4 = -t134 / 0.2e1 - t136 / 0.2e1 + t81 * t104 + t101 - t100;
t3 = -t135 / 0.2e1 - t137 / 0.2e1 + t81 * t103 - t129 / 0.2e1 - t143 * t133 / 0.2e1;
t10 = [0, 0, 0, 0, 0, -t52, t51, -t81 * t52, t126, t5 * qJD(3) + t18 * qJD(4), -t96, t56, t64, t65, 0, -t7 * qJD(5) + t125, t6 * qJD(5) + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t52 - t111, t51 + t112, -t55 * t98, -t113 + t126, t123 + (-t55 * pkin(3) - qJ(4) * t20) * qJD(3) + t9 * qJD(4), -t96, t56, t64, t65, 0, t4 * qJD(5) - t116 + t125, t3 * qJD(5) + t115 + t127; 0, 0, 0, 0, 0, 0, 0, 0, t144, t9 * qJD(3) + t114, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t145, t29, t30, 0, t4 * qJD(3) + t17 * qJD(5) - t121, t3 * qJD(3) + t16 * qJD(5) + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t117, t79 * t118; 0, 0, 0, 0, 0, t111, -t112, t81 * t111, t67 + t113, -t8 * qJD(4) - t123, -t96, t56, t64, t65, 0, -t1 * qJD(5) + t116 + t57, -t2 * qJD(5) - t115 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t67, t66 * qJD(4), -t96, t56, t64, t65, 0, -t24 * qJD(5) + t57, t23 * qJD(5) + t58; 0, 0, 0, 0, 0, 0, 0, 0, t144, -t90, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t145, t29, t30, 0, t43 * qJD(5) + t92, t42 * qJD(5) + t91; 0, 0, 0, 0, 0, 0, 0, 0, -t144, t8 * qJD(3) - t114, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t144, t90, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t145, -t44, -t45, 0, t1 * qJD(3) - t106 + t121, t2 * qJD(3) - t105 - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t145, -t44, -t45, 0, -t92 - t106, -t91 - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
