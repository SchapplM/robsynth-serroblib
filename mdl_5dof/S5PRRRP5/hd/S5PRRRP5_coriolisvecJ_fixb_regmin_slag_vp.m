% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:22
% DurationCPUTime: 0.77s
% Computational Cost: add. (957->149), mult. (2441->221), div. (0->0), fcn. (1672->6), ass. (0->106)
t71 = sin(qJ(4));
t72 = sin(qJ(3));
t74 = cos(qJ(4));
t75 = cos(qJ(3));
t51 = t71 * t75 + t74 * t72;
t42 = t51 * qJD(2);
t114 = t42 * qJ(5);
t73 = sin(qJ(2));
t107 = t73 * qJD(1);
t58 = qJD(2) * pkin(6) + t107;
t94 = pkin(7) * qJD(2) + t58;
t35 = t94 * t75;
t28 = t71 * t35;
t34 = t94 * t72;
t31 = qJD(3) * pkin(3) - t34;
t96 = t74 * t31 - t28;
t139 = t114 - t96;
t68 = qJD(3) + qJD(4);
t138 = t68 * t51;
t19 = t138 * qJD(2);
t105 = qJD(2) * qJD(3);
t137 = -0.2e1 * t105;
t123 = t74 * t75;
t124 = t71 * t72;
t50 = -t123 + t124;
t37 = t50 * t73;
t76 = cos(qJ(2));
t106 = t76 * qJD(1);
t111 = qJD(3) * t72;
t24 = -t58 * t111 + (-pkin(7) * t111 + t75 * t106) * qJD(2);
t136 = (qJD(4) * t31 + t24) * t74;
t131 = pkin(6) + pkin(7);
t100 = qJD(3) * t131;
t52 = t72 * t100;
t53 = t75 * t100;
t54 = t131 * t72;
t55 = t131 * t75;
t88 = t71 * t54 - t74 * t55;
t135 = t88 * qJD(4) + t51 * t106 + t71 * t52 - t74 * t53;
t108 = qJD(4) * t74;
t109 = qJD(4) * t71;
t84 = t50 * t76;
t134 = -qJD(1) * t84 + t54 * t108 + t55 * t109 + t74 * t52 + t71 * t53;
t77 = qJD(3) ^ 2;
t78 = qJD(2) ^ 2;
t133 = (t77 + t78) * t73;
t132 = t42 ^ 2;
t5 = t68 * pkin(4) - t139;
t130 = t5 + t139;
t129 = -qJ(5) * t138 - t50 * qJD(5) - t134;
t110 = qJD(3) * t75;
t90 = t68 * t124;
t22 = -t75 * t108 - t74 * t110 + t90;
t128 = t22 * qJ(5) - t51 * qJD(5) + t135;
t102 = qJD(2) * t123;
t113 = qJD(2) * t72;
t103 = t71 * t113;
t40 = -t102 + t103;
t67 = -t75 * pkin(3) - pkin(2);
t47 = t67 * qJD(2) - t106;
t21 = t40 * pkin(4) + qJD(5) + t47;
t127 = t21 * t42;
t126 = t42 * t40;
t125 = t47 * t42;
t30 = t74 * t35;
t122 = t77 * t72;
t121 = t77 * t75;
t120 = -t74 * t34 - t28;
t99 = t75 * t105;
t119 = -qJD(4) * t102 - t74 * t99;
t104 = pkin(3) * t113;
t48 = qJD(2) * t107 + qJD(3) * t104;
t118 = t72 ^ 2 - t75 ^ 2;
t116 = qJD(2) * pkin(2);
t115 = t40 * qJ(5);
t112 = qJD(2) * t73;
t101 = -pkin(3) * t68 - t31;
t25 = -t58 * t110 + (-pkin(7) * t110 - t72 * t106) * qJD(2);
t98 = -t71 * t24 + t74 * t25;
t97 = -t35 * t109 + t71 * t25;
t95 = t71 * t34 - t30;
t91 = t76 * t137;
t14 = t19 * pkin(4) + t48;
t89 = -t71 * t31 - t30;
t87 = qJD(2) * t116;
t86 = t47 * t40 - t97;
t85 = pkin(3) * t111 - t107;
t82 = -0.2e1 * qJD(3) * t116;
t81 = t89 * qJD(4) + t98;
t66 = t74 * pkin(3) + pkin(4);
t39 = t40 ^ 2;
t36 = t51 * t73;
t18 = qJD(2) * t90 + t119;
t17 = -t50 * qJ(5) - t88;
t16 = -t51 * qJ(5) - t74 * t54 - t71 * t55;
t15 = -t39 + t132;
t13 = t42 * t68 - t19;
t12 = -t119 + (-t103 + t40) * t68;
t11 = t68 * t37 - t76 * t42;
t10 = -qJD(2) * t84 - t138 * t73;
t9 = -t114 + t120;
t8 = t95 + t115;
t7 = -t89 - t115;
t2 = t18 * qJ(5) - t42 * qJD(5) + t81;
t1 = -t19 * qJ(5) - t40 * qJD(5) + t136 + t97;
t3 = [0, 0, -t78 * t73, -t78 * t76, 0, 0, 0, 0, 0, -t75 * t133 + t72 * t91, t72 * t133 + t75 * t91, 0, 0, 0, 0, 0, t11 * t68 + t40 * t112 - t76 * t19, -t10 * t68 + t42 * t112 + t76 * t18, -t10 * t40 - t11 * t42 - t36 * t18 + t37 * t19, -t1 * t37 + t7 * t10 + t5 * t11 + t21 * t112 - t14 * t76 - t2 * t36; 0, 0, 0, 0, 0.2e1 * t72 * t99, t118 * t137, t121, -t122, 0, -pkin(6) * t121 + t72 * t82, pkin(6) * t122 + t75 * t82, -t18 * t51 - t42 * t22, -t138 * t42 + t18 * t50 - t51 * t19 + t22 * t40, -t22 * t68, -t138 * t68, 0, t135 * t68 + t138 * t47 + t67 * t19 + t85 * t40 + t48 * t50, t134 * t68 - t67 * t18 - t47 * t22 + t85 * t42 + t48 * t51, -t1 * t50 - t128 * t42 - t129 * t40 - t138 * t7 + t16 * t18 - t17 * t19 - t2 * t51 + t5 * t22, t1 * t17 + t2 * t16 + t14 * (t50 * pkin(4) + t67) + t129 * t7 + t128 * t5 + (pkin(4) * t138 + t85) * t21; 0, 0, 0, 0, -t72 * t78 * t75, t118 * t78, 0, 0, 0, t72 * t87, t75 * t87, t126, t15, t12, t13, 0, -t95 * t68 - t40 * t104 - t125 + (t101 * t71 - t30) * qJD(4) + t98, t120 * t68 - t42 * t104 + (t101 * qJD(4) - t24) * t74 + t86, t66 * t18 + (t7 + t8) * t42 + (-t5 + t9) * t40 + (-t19 * t71 + (-t40 * t74 + t42 * t71) * qJD(4)) * pkin(3), -pkin(4) * t127 + t2 * t66 - t5 * t8 - t7 * t9 + (-t21 * t113 + t1 * t71 + (-t5 * t71 + t7 * t74) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t15, t12, t13, 0, -t89 * t68 - t125 + t81, t96 * t68 - t136 + t86, pkin(4) * t18 - t130 * t40, t130 * t7 + (t2 - t127) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 - t132, t7 * t40 + t5 * t42 + t14;];
tauc_reg = t3;
