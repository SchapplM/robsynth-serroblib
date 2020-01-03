% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:30
% EndTime: 2019-12-31 16:56:32
% DurationCPUTime: 0.58s
% Computational Cost: add. (286->112), mult. (753->175), div. (0->0), fcn. (561->4), ass. (0->100)
t51 = sin(qJ(3));
t47 = t51 ^ 2;
t118 = t47 / 0.2e1;
t50 = sin(qJ(4));
t117 = 0.2e1 * t50;
t46 = t50 ^ 2;
t52 = cos(qJ(4));
t48 = t52 ^ 2;
t32 = t48 - t46;
t53 = cos(qJ(3));
t109 = t52 * t53;
t69 = t109 * t117;
t55 = qJD(1) * t69 - t32 * qJD(3);
t115 = t51 * pkin(6);
t114 = t53 * pkin(3);
t49 = t53 ^ 2;
t43 = t49 * t52;
t25 = t114 + t115;
t113 = t50 * t25;
t54 = -pkin(1) - pkin(5);
t112 = t50 * t54;
t111 = t52 * t25;
t110 = t52 * t47;
t108 = t52 * t54;
t107 = t53 * t54;
t31 = t47 - t49;
t84 = t50 * t107;
t66 = t51 * pkin(3) - t53 * pkin(6);
t64 = qJ(2) + t66;
t9 = t51 * t112 - t52 * t64;
t1 = -t9 * t53 + (t84 + t111) * t51;
t106 = t1 * qJD(1);
t85 = t51 * t108;
t10 = t50 * t64 + t85;
t83 = t52 * t107;
t2 = t10 * t53 + (-t83 + t113) * t51;
t105 = t2 * qJD(1);
t5 = -t49 * t112 - t9 * t51;
t104 = t5 * qJD(1);
t6 = -t10 * t51 - t49 * t108;
t103 = t6 * qJD(1);
t102 = qJD(2) * t51;
t101 = qJD(3) * t51;
t100 = qJD(3) * t52;
t99 = qJD(3) * t53;
t98 = qJD(3) * t54;
t97 = qJD(4) * t50;
t96 = qJD(4) * t52;
t82 = 0.1e1 / 0.2e1 + t118;
t11 = (-t49 / 0.2e1 - t82) * t50;
t95 = t11 * qJD(1);
t12 = t43 / 0.2e1 + t82 * t52;
t94 = t12 * qJD(1);
t22 = t31 * t50;
t93 = t22 * qJD(1);
t23 = -t43 + t110;
t92 = t23 * qJD(1);
t91 = t31 * qJD(1);
t90 = t51 * qJD(1);
t89 = t53 * qJD(1);
t88 = t53 * qJD(4);
t87 = qJ(2) * qJD(3);
t86 = qJD(1) * qJ(2);
t81 = t51 * t97;
t80 = t50 * t88;
t79 = t51 * t96;
t78 = t52 * t88;
t77 = t50 * t99;
t76 = t50 * t96;
t75 = t50 * t100;
t74 = t52 * t99;
t73 = t51 * t99;
t72 = t51 * t89;
t71 = t51 * t86;
t70 = t53 * t86;
t67 = qJD(3) * t69;
t65 = (-qJD(4) - t90) * t53;
t63 = t115 / 0.2e1 + t114 / 0.2e1;
t58 = t25 / 0.2e1 + t63;
t8 = t58 * t52;
t62 = pkin(3) * t50 * qJD(3) + t8 * qJD(1);
t7 = t58 * t50;
t61 = pkin(3) * t100 - t7 * qJD(1);
t60 = t52 * t65;
t16 = (t46 / 0.2e1 - t48 / 0.2e1) * t53;
t59 = -t16 * qJD(1) + t75;
t57 = t50 * qJD(1) * t43 + t16 * qJD(3);
t21 = t32 * t49;
t56 = t21 * qJD(1) + t67;
t42 = t99 / 0.2e1;
t41 = t52 * t90;
t40 = t50 * t101;
t39 = t50 * t90;
t20 = (t90 + qJD(4) / 0.2e1) * t53;
t15 = t16 * qJD(4);
t14 = -t110 / 0.2e1 - t43 / 0.2e1 + t52 / 0.2e1;
t13 = (-0.1e1 / 0.2e1 + t118 + t49 / 0.2e1) * t50;
t4 = -t84 + t111 / 0.2e1 - t63 * t52;
t3 = -t83 - t113 / 0.2e1 + t63 * t50;
t17 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t73, t31 * qJD(3), 0, 0, 0, t53 * t87 + t102, qJD(2) * t53 - t51 * t87, -t48 * t73 - t49 * t76, -t21 * qJD(4) + t51 * t67, -t23 * qJD(3) - t51 * t80, t22 * qJD(3) - t51 * t78, t73, t1 * qJD(3) + t6 * qJD(4) + t52 * t102, -t2 * qJD(3) - t5 * qJD(4) - t50 * t102; 0, 0, 0, 0, qJD(1), t86, 0, 0, 0, 0, 0, t90, t89, 0, 0, 0, 0, 0, t14 * qJD(4) + t41, t13 * qJD(4) - t39; 0, 0, 0, 0, 0, 0, -t72, t91, -t101, -t99, 0, -t51 * t98 + t70, -t53 * t98 - t71, -t15 + (-t48 * t89 - t75) * t51, t55 * t51 - 0.2e1 * t53 * t76, t77 - t92, t74 + t93, t20, t106 + (t66 * t50 - t85) * qJD(3) + t4 * qJD(4), -t105 + (-pkin(6) * t109 + (pkin(3) * t52 + t112) * t51) * qJD(3) + t3 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56, t50 * t65, t60, t42, t14 * qJD(2) + t4 * qJD(3) - t10 * qJD(4) + t103, t13 * qJD(2) + t3 * qJD(3) + t9 * qJD(4) - t104; 0, 0, 0, 0, -qJD(1), -t86, 0, 0, 0, 0, 0, -t90, -t89, 0, 0, 0, 0, 0, -t12 * qJD(4) - t41, -t11 * qJD(4) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t99, 0, 0, 0, 0, 0, -t51 * t100 - t80, t40 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77 - t79 - t94, -t74 + t81 - t95; 0, 0, 0, 0, 0, 0, t72, -t91, 0, 0, 0, -t70, t71, t48 * t72 - t15, t60 * t117, t79 + t92, -t81 - t93, -t20, -t8 * qJD(4) - t106, t7 * qJD(4) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t32 * qJD(4), 0, 0, 0, -pkin(3) * t97, -pkin(3) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t55, t41 + t96, -t39 - t97, -t89 / 0.2e1, -pkin(6) * t96 - t62, pkin(6) * t97 - t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t56, (t50 * t89 - t100) * t51, t52 * t72 + t40, t42, t12 * qJD(2) + t8 * qJD(3) - t103, t11 * qJD(2) - t7 * qJD(3) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t55, -t41, t39, t89 / 0.2e1, t62, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t17;
