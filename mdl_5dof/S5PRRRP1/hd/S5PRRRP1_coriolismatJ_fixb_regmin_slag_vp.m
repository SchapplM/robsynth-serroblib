% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP1
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
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:16
% EndTime: 2019-12-05 16:40:19
% DurationCPUTime: 0.49s
% Computational Cost: add. (584->94), mult. (1221->132), div. (0->0), fcn. (846->4), ass. (0->88)
t75 = cos(qJ(3));
t114 = t75 * pkin(2);
t63 = -pkin(3) - t114;
t119 = t63 / 0.2e1 - pkin(3) / 0.2e1;
t72 = sin(qJ(4));
t70 = t72 ^ 2;
t74 = cos(qJ(4));
t71 = t74 ^ 2;
t59 = t70 + t71;
t97 = qJD(2) + qJD(3);
t121 = t97 * t59;
t60 = t71 - t70;
t120 = t97 * t60;
t116 = pkin(4) * t72;
t73 = sin(qJ(3));
t115 = t73 * pkin(2);
t62 = pkin(7) + t115;
t102 = qJ(5) + t62;
t39 = t102 * t72;
t113 = t39 * t72;
t40 = t102 * t74;
t112 = t40 * t74;
t109 = -qJ(5) - pkin(7);
t51 = t109 * t72;
t111 = t51 * t72;
t52 = t109 * t74;
t110 = t52 * t74;
t38 = t59 * t114;
t53 = t59 * qJD(5);
t108 = t38 * qJD(3) + t53;
t107 = pkin(2) * qJD(2);
t106 = pkin(2) * qJD(3);
t105 = pkin(3) * qJD(3);
t104 = qJD(4) * pkin(4);
t13 = t112 + t113;
t64 = -t74 * pkin(4) - pkin(3);
t50 = t64 - t114;
t8 = (t13 * t75 + t50 * t73) * pkin(2);
t103 = t8 * qJD(2);
t101 = qJD(2) * t63;
t100 = qJD(4) * t72;
t69 = qJD(4) * t74;
t99 = t13 * qJD(2);
t98 = t38 * qJD(2);
t96 = qJD(5) * t116;
t95 = t73 * t106;
t94 = pkin(4) * t69;
t93 = t73 * t107;
t91 = t115 / 0.2e1;
t90 = -t114 / 0.2e1;
t89 = t72 * t101;
t88 = t74 * t101;
t87 = pkin(2) * t97;
t86 = t72 * t93;
t85 = t72 * t90;
t84 = t97 * t72;
t83 = t73 * t87;
t17 = -t110 - t111;
t33 = t112 / 0.2e1;
t11 = t33 - t112 / 0.2e1;
t3 = t50 * t116;
t82 = t11 * qJD(1) + t3 * qJD(2);
t4 = t91 + (t52 / 0.2e1 - t40 / 0.2e1) * t74 + (t51 / 0.2e1 - t39 / 0.2e1) * t72;
t81 = -t4 * qJD(2) + t17 * qJD(3);
t44 = -t110 / 0.2e1;
t15 = t44 + t110 / 0.2e1;
t80 = -t11 * qJD(2) - t15 * qJD(3);
t79 = t90 - t119;
t22 = t79 * t72;
t78 = t22 * qJD(2) + t72 * t105;
t23 = t79 * t74;
t77 = t23 * qJD(2) + t74 * t105;
t1 = (t90 - t64 / 0.2e1 - t50 / 0.2e1) * t116;
t9 = t64 * t116;
t76 = t15 * qJD(1) - t1 * qJD(2) + t9 * qJD(3);
t68 = pkin(4) * t100;
t61 = t72 * t69;
t58 = t72 * t95;
t54 = t60 * qJD(4);
t46 = pkin(4) * t84;
t37 = t74 * t84;
t25 = (t90 + t119) * t74;
t24 = t119 * t72 + t85;
t14 = t15 * qJD(4);
t10 = t11 * qJD(4);
t5 = t44 + t33 - t111 / 0.2e1 + t113 / 0.2e1 + t91;
t2 = pkin(4) * t85 + (t50 + t64) * t116 / 0.2e1;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t69, 0, -t68 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, -t95, -t75 * t106, t61, t54, 0, 0, 0, t63 * t100 - t74 * t95, t63 * t69 + t58, t108, t8 * qJD(3) + t3 * qJD(4) + t13 * qJD(5); 0, 0, 0, 0, 0, -t83, -t75 * t87, t61, t54, 0, 0, 0, t24 * qJD(4) - t74 * t83, t25 * qJD(4) + t58 + t86, t98 + t108, t103 + t2 * qJD(4) + t5 * qJD(5) + (t17 * t75 + t64 * t73) * t106; 0, 0, 0, 0, 0, 0, 0, t37, t120, t69, -t100, 0, t24 * qJD(3) - t62 * t69 + t89, t25 * qJD(3) + t62 * t100 + t88, -t94, t2 * qJD(3) - t40 * t104 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t5 * qJD(3) + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, t93, t75 * t107, t61, t54, 0, 0, 0, -t22 * qJD(4) + t74 * t93, -t23 * qJD(4) - t86, t53 - t98, -t1 * qJD(4) - t4 * qJD(5) - t103; 0, 0, 0, 0, 0, 0, 0, t61, t54, 0, 0, 0, -pkin(3) * t100, -pkin(3) * t69, t53, t9 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t37, t120, t69, -t100, 0, -pkin(7) * t69 - t78, pkin(7) * t100 - t77, -t94, t52 * t104 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80; 0, 0, 0, 0, 0, 0, 0, -t37, -t120, 0, 0, 0, t22 * qJD(3) - t89, t23 * qJD(3) - t88, 0, t1 * qJD(3) - t82 - t96; 0, 0, 0, 0, 0, 0, 0, -t37, -t120, 0, 0, 0, t78, t77, 0, -t76 - t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t4 * qJD(3) + t68 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t68 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
