% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:46
% EndTime: 2019-12-05 18:01:48
% DurationCPUTime: 0.54s
% Computational Cost: add. (948->100), mult. (1776->136), div. (0->0), fcn. (1395->6), ass. (0->92)
t120 = cos(qJ(3));
t123 = pkin(1) * sin(pkin(8));
t80 = sin(qJ(3));
t92 = cos(pkin(8)) * pkin(1) + pkin(2);
t51 = -t120 * t92 + t80 * t123;
t126 = t51 / 0.2e1;
t103 = qJD(1) + qJD(3);
t79 = sin(qJ(4));
t76 = t79 ^ 2;
t81 = cos(qJ(4));
t77 = t81 ^ 2;
t68 = t76 + t77;
t129 = t103 * t68;
t69 = t77 - t76;
t128 = t103 * t69;
t127 = -pkin(3) / 0.2e1;
t125 = t79 / 0.2e1;
t122 = pkin(4) * t79;
t121 = pkin(4) * t81;
t86 = t80 * t92;
t94 = t120 * t123;
t52 = t94 + t86;
t50 = pkin(7) + t52;
t110 = qJ(5) + t50;
t32 = t110 * t79;
t119 = t32 * t79;
t33 = t110 * t81;
t118 = t33 * t81;
t115 = -qJ(5) - pkin(7);
t62 = t115 * t79;
t117 = t62 * t79;
t63 = t115 * t81;
t116 = t63 * t81;
t19 = t68 * t51;
t65 = t68 * qJD(5);
t114 = -t19 * qJD(3) + t65;
t113 = pkin(3) * qJD(3);
t112 = qJD(4) * pkin(4);
t49 = -pkin(3) + t51;
t39 = t49 - t121;
t91 = -t118 - t119;
t3 = t39 * t52 + t91 * t51;
t111 = t3 * qJD(1);
t109 = qJD(1) * t91;
t108 = qJD(1) * t19;
t107 = qJD(1) * t49;
t106 = qJD(4) * t79;
t75 = qJD(4) * t81;
t105 = t51 * qJD(1);
t104 = t52 * qJD(1);
t48 = t52 * qJD(3);
t102 = pkin(4) * t75;
t101 = qJD(5) * t122;
t99 = t81 * t107;
t98 = t79 * t104;
t97 = t79 * t107;
t96 = t51 * t125;
t95 = t103 * t79;
t93 = t126 + pkin(3) / 0.2e1 - t49 / 0.2e1;
t90 = t116 + t117;
t26 = t118 / 0.2e1;
t10 = t26 - t118 / 0.2e1;
t4 = t39 * t122;
t89 = qJD(1) * t4 + qJD(2) * t10;
t82 = t86 / 0.2e1 + t94 / 0.2e1;
t5 = (t63 / 0.2e1 - t33 / 0.2e1) * t81 + (t62 / 0.2e1 - t32 / 0.2e1) * t79 + t82;
t88 = -qJD(1) * t5 - qJD(3) * t90;
t57 = -t116 / 0.2e1;
t31 = t57 + t116 / 0.2e1;
t87 = -qJD(1) * t10 - qJD(3) * t31;
t14 = t93 * t79;
t85 = t14 * qJD(1) + t79 * t113;
t15 = t93 * t81;
t84 = t15 * qJD(1) + t81 * t113;
t71 = -pkin(3) - t121;
t1 = (t126 - t71 / 0.2e1 - t39 / 0.2e1) * t122;
t13 = t71 * t122;
t83 = -qJD(1) * t1 + qJD(2) * t31 + qJD(3) * t13;
t74 = pkin(4) * t106;
t70 = t79 * t75;
t66 = t69 * qJD(4);
t59 = pkin(4) * t95;
t53 = t81 * t95;
t47 = t51 * qJD(3);
t42 = t79 * t48;
t29 = t31 * qJD(4);
t17 = (t127 + t49 / 0.2e1 + t126) * t81;
t16 = t49 * t125 + t79 * t127 + t96;
t9 = t10 * qJD(4);
t6 = t57 + t26 - t117 / 0.2e1 + t119 / 0.2e1 + t82;
t2 = pkin(4) * t96 + (t39 + t71) * t122 / 0.2e1;
t7 = [0, 0, 0, 0, 0, -t48, t47, t70, t66, 0, 0, 0, t49 * t106 - t81 * t48, t49 * t75 + t42, t114, qJD(3) * t3 + qJD(4) * t4 - qJD(5) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, -t48 - t104, t47 + t105, t70, t66, 0, 0, 0, -t103 * t81 * t52 + t16 * qJD(4), qJD(4) * t17 + t42 + t98, -t108 + t114, t111 + (t90 * t51 + t52 * t71) * qJD(3) + t2 * qJD(4) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t53, t128, t75, -t106, 0, t16 * qJD(3) - t50 * t75 + t97, t17 * qJD(3) + t50 * t106 + t99, -t102, qJD(3) * t2 - t33 * t112 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJD(3) * t6 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t75, 0, -t74 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t104, -t105, t70, t66, 0, 0, 0, -t14 * qJD(4) + t81 * t104, -qJD(4) * t15 - t98, t65 + t108, -qJD(4) * t1 - qJD(5) * t5 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, t70, t66, 0, 0, 0, -pkin(3) * t106, -pkin(3) * t75, t65, qJD(4) * t13 - qJD(5) * t90; 0, 0, 0, 0, 0, 0, 0, t53, t128, t75, -t106, 0, -pkin(7) * t75 - t85, pkin(7) * t106 - t84, -t102, t63 * t112 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t88; 0, 0, 0, 0, 0, 0, 0, -t53, -t128, 0, 0, 0, qJD(3) * t14 - t97, t15 * qJD(3) - t99, 0, qJD(3) * t1 - t101 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, -t53, -t128, 0, 0, 0, t85, t84, 0, -t83 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, qJD(3) * t5 + t109 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t74 - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
