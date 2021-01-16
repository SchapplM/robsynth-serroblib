% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:54
% EndTime: 2021-01-15 15:04:57
% DurationCPUTime: 0.58s
% Computational Cost: add. (635->94), mult. (1459->149), div. (0->0), fcn. (1245->4), ass. (0->104)
t68 = sin(qJ(4));
t119 = qJ(3) * t68;
t67 = cos(pkin(8));
t66 = sin(pkin(8));
t118 = qJ(5) * t66;
t48 = -t67 * pkin(3) - t66 * pkin(6) - pkin(2);
t69 = cos(qJ(4));
t39 = t69 * t48;
t72 = -t69 * t118 + t39;
t20 = (-pkin(4) - t119) * t67 + t72;
t93 = t67 * t119;
t23 = -t72 + t93;
t126 = -t20 - t23;
t64 = t68 ^ 2;
t65 = t69 ^ 2;
t125 = t64 + t65;
t123 = t67 * pkin(4);
t122 = t20 * t68;
t121 = t66 * t68;
t120 = t66 * t69;
t62 = t66 ^ 2;
t54 = t67 ^ 2 + t62;
t80 = t20 / 0.2e1 + t23 / 0.2e1;
t2 = (t123 / 0.2e1 + t80) * t120;
t117 = t2 * qJD(2);
t3 = t126 * t121;
t116 = t3 * qJD(2);
t92 = -t123 / 0.2e1;
t4 = (t92 + t80) * t68;
t115 = t4 * qJD(2);
t114 = t62 * qJ(3);
t29 = -t69 * t67 * qJ(3) - t68 * t48;
t24 = -t68 * t118 - t29;
t8 = (t20 * t69 + t24 * t68) * t66;
t113 = t8 * qJD(2);
t94 = t62 * t68 * t69;
t46 = (pkin(4) * t68 + qJ(3)) * t66;
t95 = t46 * t120;
t9 = -pkin(4) * t94 - t24 * t67 - t95;
t112 = t9 * qJD(2);
t111 = qJD(2) * t68;
t110 = qJD(2) * t69;
t109 = qJD(3) * t67;
t108 = qJD(4) * t68;
t107 = qJD(4) * t69;
t106 = qJD(5) * t66;
t10 = t65 * t62 * pkin(4) - t46 * t121 - t23 * t67;
t105 = t10 * qJD(2);
t28 = -t39 + t93;
t21 = -t68 * t114 - t28 * t67;
t104 = t21 * qJD(2);
t22 = -t69 * t114 + t29 * t67;
t103 = t22 * qJD(2);
t102 = t24 * qJD(4);
t79 = t64 / 0.2e1 + t65 / 0.2e1;
t25 = (-0.1e1 / 0.2e1 + t79) * t67 * t66;
t101 = t25 * qJD(2);
t33 = (0.1e1 / 0.2e1 + t79) * t66;
t100 = t33 * qJD(2);
t40 = t125 * t62;
t99 = t40 * qJD(2);
t41 = (t64 - t65) * t62;
t98 = t41 * qJD(2);
t42 = t54 * t68;
t36 = t42 * qJD(2);
t43 = t54 * t69;
t38 = t43 * qJD(2);
t51 = t54 * qJ(3);
t97 = t51 * qJD(2);
t96 = t54 * qJD(2);
t58 = t66 * t108;
t91 = t67 * t108;
t90 = t66 * t107;
t89 = t67 * t107;
t88 = t66 * t111;
t87 = t68 * t106;
t86 = t67 * t111;
t85 = t68 * t109;
t84 = t66 * t110;
t83 = t69 * t106;
t82 = t67 * t110;
t81 = t69 * t109;
t78 = pkin(4) * t90;
t77 = pkin(4) * t84;
t76 = t66 * t86;
t75 = qJD(2) * t94;
t74 = t66 * t82;
t73 = qJD(2) * t67 - qJD(4);
t6 = pkin(4) * t95 + t126 * t24;
t71 = -t2 * qJD(1) + t6 * qJD(2);
t7 = t46 * t66 + (t24 * t69 - t122) * t67;
t70 = -t25 * qJD(1) - t7 * qJD(2);
t45 = t73 * t69;
t44 = t73 * t68;
t37 = t43 * qJD(3);
t35 = t42 * qJD(3);
t34 = (0.1e1 / 0.2e1 - t125 / 0.2e1) * t66;
t31 = t66 * t45;
t30 = -t58 + t76;
t27 = t89 - t38;
t26 = t91 - t36;
t5 = -t122 / 0.2e1 + (-t23 / 0.2e1 + t92) * t68;
t1 = t25 * qJD(3) - t2 * qJD(4);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t58, -t90, t58, 0, -t78 - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, t54 * qJD(3), t51 * qJD(3), -qJD(4) * t94, t41 * qJD(4), t66 * t91, t66 * t89, 0, -t22 * qJD(4) + t35, t21 * qJD(4) + t37, -t9 * qJD(4) + t67 * t83 + t35, t10 * qJD(4) - t67 * t87 + t37, -t3 * qJD(4) + t40 * qJD(5), t7 * qJD(3) + t6 * qJD(4) - t8 * qJD(5); 0, 0, 0, 0, 0, t96, t97, 0, 0, 0, 0, 0, t36, t38, t36, t38, 0, t5 * qJD(4) + t34 * qJD(5) - t70; 0, 0, 0, 0, 0, 0, 0, -t75, t98, t30, t31, 0, t29 * qJD(4) - t103, t28 * qJD(4) + t104, -t102 - t112, t23 * qJD(4) + t105, pkin(4) * t58 - t116, -pkin(4) * t102 + t5 * qJD(3) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t76, t99, t34 * qJD(3) - t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101; 0, 0, 0, 0, 0, -t96, -t97, 0, 0, 0, 0, 0, t26, t27, t26, t27, 0, -t4 * qJD(4) - t33 * qJD(5) + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, t44, t45, 0, -pkin(4) * t108 - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117; 0, 0, 0, 0, 0, 0, 0, t75, -t98, -t76, -t74, 0, -t85 + t103, -t81 - t104, -t83 - t85 + t112, -t81 + t87 - t105, t116, -pkin(4) * t83 + t4 * qJD(3) - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t82, -t86, -t82, 0, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t88, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t30, -t99, t33 * qJD(3) + t113 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t88, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
