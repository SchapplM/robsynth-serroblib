% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR4
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
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:57
% EndTime: 2019-12-05 17:54:01
% DurationCPUTime: 0.77s
% Computational Cost: add. (1023->147), mult. (2552->217), div. (0->0), fcn. (1806->8), ass. (0->101)
t94 = cos(pkin(9));
t99 = cos(qJ(3));
t126 = t94 * t99;
t92 = sin(pkin(9));
t97 = sin(qJ(3));
t79 = -t92 * t97 + t126;
t72 = t79 * qJD(1);
t98 = cos(qJ(5));
t63 = t98 * t72;
t80 = t92 * t99 + t94 * t97;
t74 = t80 * qJD(1);
t96 = sin(qJ(5));
t29 = -t96 * t74 + t63;
t89 = qJD(3) + qJD(5);
t128 = t29 * t89;
t120 = qJD(5) * t96;
t73 = t80 * qJD(3);
t66 = qJD(1) * t73;
t117 = qJD(1) * qJD(3);
t113 = t97 * t117;
t67 = -t113 * t92 + t117 * t126;
t4 = qJD(5) * t63 - t120 * t74 - t96 * t66 + t98 * t67;
t141 = t4 - t128;
t108 = t96 * t72 + t98 * t74;
t140 = t108 * t29;
t85 = sin(pkin(8)) * pkin(1) + pkin(6);
t121 = qJ(4) + t85;
t129 = t108 * t89;
t5 = qJD(5) * t108 + t98 * t66 + t96 * t67;
t139 = -t5 + t129;
t138 = t108 ^ 2 - t29 ^ 2;
t87 = -cos(pkin(8)) * pkin(1) - pkin(2);
t105 = -t99 * pkin(3) + t87;
t103 = t105 * qJD(1);
t71 = qJD(4) + t103;
t38 = -t72 * pkin(4) + t71;
t132 = t72 * pkin(7);
t110 = t121 * qJD(1);
t54 = t97 * qJD(2) + t110 * t99;
t127 = t94 * t54;
t124 = qJD(3) * pkin(3);
t53 = t99 * qJD(2) - t110 * t97;
t48 = t53 + t124;
t19 = t92 * t48 + t127;
t9 = t19 + t132;
t137 = t9 * t120 - t38 * t29;
t135 = qJD(5) - t89;
t116 = qJD(1) * qJD(4);
t36 = t53 * qJD(3) + t99 * t116;
t37 = -t54 * qJD(3) - t97 * t116;
t10 = -t92 * t36 + t94 * t37;
t2 = -t67 * pkin(7) + t10;
t11 = t94 * t36 + t92 * t37;
t3 = -t66 * pkin(7) + t11;
t134 = -t38 * t108 + t98 * t2 - t96 * t3;
t133 = pkin(3) * t92;
t131 = t74 * pkin(7);
t107 = t98 * t79 - t96 * t80;
t76 = t79 * qJD(3);
t12 = qJD(5) * t107 - t96 * t73 + t98 * t76;
t130 = t12 * t89;
t44 = t92 * t54;
t21 = t94 * t53 - t44;
t111 = qJD(3) * t121;
t57 = t99 * qJD(4) - t111 * t97;
t58 = -t97 * qJD(4) - t111 * t99;
t23 = t94 * t57 + t92 * t58;
t77 = t121 * t97;
t78 = t121 * t99;
t35 = -t92 * t77 + t94 * t78;
t125 = t97 ^ 2 - t99 ^ 2;
t100 = qJD(3) ^ 2;
t123 = t100 * t97;
t122 = t100 * t99;
t82 = qJD(1) * t87;
t119 = t97 * qJD(1);
t115 = t97 * t124;
t18 = t94 * t48 - t44;
t20 = -t92 * t53 - t127;
t22 = -t92 * t57 + t94 * t58;
t34 = -t94 * t77 - t92 * t78;
t7 = qJD(3) * pkin(4) - t131 + t18;
t109 = -t96 * t7 - t98 * t9;
t40 = t96 * t79 + t98 * t80;
t104 = 0.2e1 * qJD(3) * t82;
t101 = qJD(1) ^ 2;
t86 = t94 * pkin(3) + pkin(4);
t84 = pkin(3) * t113;
t56 = t73 * pkin(4) + t115;
t55 = pkin(3) * t119 + t74 * pkin(4);
t52 = -t79 * pkin(4) + t105;
t43 = t66 * pkin(4) + t84;
t25 = t79 * pkin(7) + t35;
t24 = -t80 * pkin(7) + t34;
t17 = -t73 * pkin(7) + t23;
t16 = -t76 * pkin(7) + t22;
t15 = t21 - t131;
t14 = t20 - t132;
t13 = qJD(5) * t40 + t98 * t73 + t96 * t76;
t8 = t13 * t89;
t1 = [0, 0, 0, 0, 0.2e1 * t99 * t113, -0.2e1 * t125 * t117, t122, -t123, 0, t104 * t97 - t122 * t85, t104 * t99 + t123 * t85, -t10 * t80 + t11 * t79 - t18 * t76 - t19 * t73 - t22 * t74 + t23 * t72 - t34 * t67 - t35 * t66, t10 * t34 + t11 * t35 + t18 * t22 + t19 * t23 + (t71 + t103) * t115, t108 * t12 + t4 * t40, t107 * t4 - t108 * t13 + t12 * t29 - t40 * t5, t130, -t8, 0, -t56 * t29 + t52 * t5 - t43 * t107 + t38 * t13 + (t98 * t16 - t96 * t17 + (-t24 * t96 - t25 * t98) * qJD(5)) * t89, t56 * t108 + t52 * t4 + t43 * t40 + t38 * t12 - (t96 * t16 + t98 * t17 + (t24 * t98 - t25 * t96) * qJD(5)) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t122, -t80 * t66 - t79 * t67 + t76 * t72 + t73 * t74, t10 * t79 + t11 * t80 - t18 * t73 + t19 * t76, 0, 0, 0, 0, 0, -t8, -t130; 0, 0, 0, 0, -t97 * t101 * t99, t125 * t101, 0, 0, 0, -t82 * t119, -t82 * t99 * qJD(1), (t19 + t20) * t74 + (t18 - t21) * t72 + (-t66 * t92 - t67 * t94) * pkin(3), -t18 * t20 - t19 * t21 + (t10 * t94 + t11 * t92 - t119 * t71) * pkin(3), -t140, t138, t141, t139, 0, t55 * t29 - (t98 * t14 - t96 * t15) * t89 + ((-t98 * t133 - t86 * t96) * t89 + t109) * qJD(5) + t134, -t98 * t3 - t96 * t2 - t55 * t108 + (t96 * t14 + t98 * t15) * t89 + (-(-t96 * t133 + t86 * t98) * t89 - t98 * t7) * qJD(5) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 ^ 2 - t74 ^ 2, t18 * t74 - t19 * t72 + t84, 0, 0, 0, 0, 0, t5 + t129, t4 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t138, t141, t139, 0, t135 * t109 + t134, (-t9 * t89 - t2) * t96 + (-t135 * t7 - t3) * t98 + t137;];
tauc_reg = t1;
