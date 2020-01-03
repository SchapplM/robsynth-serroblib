% Calculate inertial parameters regressor of coriolis matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:18
% EndTime: 2019-12-31 17:14:19
% DurationCPUTime: 0.82s
% Computational Cost: add. (553->142), mult. (1209->149), div. (0->0), fcn. (838->4), ass. (0->111)
t71 = cos(qJ(2));
t125 = t71 * pkin(1);
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t81 = -t70 * pkin(3) - t68 * qJ(4);
t40 = -pkin(2) + t81;
t35 = t40 - t125;
t27 = t35 * t68;
t37 = t40 * t68;
t117 = t27 / 0.2e1 + t37 / 0.2e1;
t110 = t70 * qJ(4);
t127 = t68 * pkin(3);
t42 = -t110 + t127;
t119 = t42 * t70;
t132 = t119 - t117;
t66 = t68 ^ 2;
t67 = t70 ^ 2;
t131 = t66 + t67;
t60 = -pkin(2) - t125;
t93 = pkin(2) / 0.2e1 - t60 / 0.2e1;
t130 = t93 * t68;
t100 = qJD(1) + qJD(2);
t56 = t67 - t66;
t129 = t100 * t56;
t128 = pkin(2) * t70;
t69 = sin(qJ(2));
t126 = t69 * pkin(1);
t124 = t35 * t42;
t123 = t35 * t70;
t122 = t40 * t70;
t121 = t42 * t40;
t120 = t42 * t68;
t118 = t60 * t70;
t115 = t131 * pkin(6) * t125;
t112 = pkin(1) * qJD(2);
t99 = t69 * t112;
t53 = t68 * t99;
t63 = t66 * qJD(4);
t114 = t63 - t53;
t113 = pkin(1) * qJD(1);
t59 = pkin(6) + t126;
t85 = t131 * t71;
t82 = t59 * t85;
t7 = (t35 * t69 + t82) * pkin(1);
t111 = t7 * qJD(1);
t109 = qJD(1) * t68;
t108 = qJD(2) * t68;
t10 = (t60 * t69 + t82) * pkin(1);
t107 = t10 * qJD(1);
t12 = t120 + t123;
t106 = t12 * qJD(1);
t13 = -t27 + t119;
t105 = t13 * qJD(1);
t36 = pkin(1) * t85;
t30 = t36 * qJD(1);
t104 = t56 * qJD(3);
t64 = t68 * qJD(3);
t65 = t70 * qJD(3);
t103 = t70 * qJD(4);
t102 = qJD(3) * qJ(4);
t98 = pkin(6) * t64;
t97 = t69 * t113;
t96 = pkin(6) * t65;
t95 = -t125 / 0.2e1;
t94 = t125 / 0.2e1;
t92 = qJD(1) * t124;
t91 = t35 * t109;
t90 = t60 * t109;
t89 = qJD(1) * t118;
t88 = t59 * t64;
t87 = t59 * t65;
t86 = t40 / 0.2e1 + t35 / 0.2e1;
t84 = pkin(1) * t100;
t83 = t70 * t99;
t14 = t120 + t122;
t49 = t70 * t94;
t4 = t86 * t70 + t120 + t49;
t80 = t4 * qJD(1) + t14 * qJD(2);
t15 = -t37 + t119;
t48 = t68 * t95;
t3 = t48 + t132;
t79 = t3 * qJD(1) + t15 * qJD(2);
t57 = t68 * t103;
t78 = t57 - t83;
t77 = qJD(3) * t42 - qJD(4) * t68;
t17 = t48 + t130;
t76 = pkin(2) * t108 + t17 * qJD(1);
t50 = t70 * t95;
t18 = t93 * t70 + t50;
t75 = t18 * qJD(1) + qJD(2) * t128;
t72 = (t110 / 0.2e1 - t127 / 0.2e1) * t125;
t1 = -t86 * t42 + t72;
t74 = t1 * qJD(1) - qJD(2) * t121;
t47 = t68 * t94;
t8 = t47 + t117;
t73 = t8 * qJD(1) + t40 * t108;
t28 = t81 * qJD(3) + t103;
t58 = t68 * t65;
t52 = t70 * t97;
t51 = t68 * t97;
t41 = t100 * t66;
t31 = t100 * t70 * t68;
t29 = t36 * qJD(2);
t20 = -t128 / 0.2e1 + t118 / 0.2e1 + t50;
t19 = t48 - t130;
t11 = t29 + t30;
t9 = t47 - t117;
t6 = t48 - t132;
t5 = -t123 / 0.2e1 - t120 - t122 / 0.2e1 + t49;
t2 = t121 / 0.2e1 + t124 / 0.2e1 + t72;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t71 * t112, 0, 0, t58, t104, 0, -t58, 0, 0, t60 * t64 - t83, t60 * t65 + t53, t29, t10 * qJD(2), t58, 0, -t104, 0, 0, -t58, -t13 * qJD(3) + t78, t29, -t12 * qJD(3) + t114, t7 * qJD(2) + t77 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t84, -t71 * t84, 0, 0, t58, t104, 0, -t58, 0, 0, t19 * qJD(3) - t52 - t83, t20 * qJD(3) + t51 + t53, t11, t107 + (-pkin(2) * t126 + t115) * qJD(2), t58, 0, -t104, 0, 0, -t58, t6 * qJD(3) - t52 + t78, t11, t5 * qJD(3) + t114 - t51, t111 + (t40 * t126 + t115) * qJD(2) + t2 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t129, t65, -t31, -t64, 0, t19 * qJD(2) - t87 + t90, t20 * qJD(2) + t88 + t89, 0, 0, t31, t65, -t129, 0, t64, -t31, t6 * qJD(2) - t105 - t87, t28, t5 * qJD(2) - t106 - t88, t2 * qJD(2) + t28 * t59 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t65, t41, t9 * qJD(2) + t87 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t71 * t113, 0, 0, t58, t104, 0, -t58, 0, 0, -t17 * qJD(3) + t52, -t18 * qJD(3) - t51, -t30, -t107, t58, 0, -t104, 0, 0, -t58, -t3 * qJD(3) + t52 + t57, -t30, -t4 * qJD(3) + t51 + t63, -t1 * qJD(3) - t8 * qJD(4) - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t104, 0, -t58, 0, 0, -pkin(2) * t64, -pkin(2) * t65, 0, 0, t58, 0, -t104, 0, 0, -t58, -t15 * qJD(3) + t57, 0, -t14 * qJD(3) + t63, t77 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t129, t65, -t31, -t64, 0, -t76 - t96, -t75 + t98, 0, 0, t31, t65, -t129, 0, t64, -t31, -t79 - t96, t28, -t80 - t98, t28 * pkin(6) - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t65, t41, -t73 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t129, 0, t31, 0, 0, t17 * qJD(2) - t90, t18 * qJD(2) - t89, 0, 0, -t31, 0, t129, 0, 0, t31, t3 * qJD(2) + t105, 0, t4 * qJD(2) + t106, t1 * qJD(2) - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t129, 0, t31, 0, 0, t76, t75, 0, 0, -t31, 0, t129, 0, 0, t31, t79, 0, t80, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t41, t8 * qJD(2) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t41, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
