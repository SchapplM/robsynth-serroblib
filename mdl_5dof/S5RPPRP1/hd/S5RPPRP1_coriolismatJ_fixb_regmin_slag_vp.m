% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:46
% EndTime: 2022-01-23 09:12:48
% DurationCPUTime: 0.60s
% Computational Cost: add. (861->97), mult. (1685->148), div. (0->0), fcn. (1471->6), ass. (0->106)
t60 = sin(pkin(7)) * pkin(1) + qJ(3);
t68 = sin(qJ(4));
t123 = t60 * t68;
t67 = cos(pkin(8));
t66 = sin(pkin(8));
t118 = qJ(5) * t66;
t69 = cos(qJ(4));
t70 = -cos(pkin(7)) * pkin(1) - pkin(2) - t66 * pkin(6) - t67 * pkin(3);
t34 = t69 * t70;
t73 = -t69 * t118 + t34;
t22 = (-pkin(4) - t123) * t67 + t73;
t94 = t67 * t123;
t23 = -t73 + t94;
t128 = -t22 - t23;
t64 = t68 ^ 2;
t65 = t69 ^ 2;
t127 = t64 + t65;
t126 = -t23 / 0.2e1;
t124 = t22 * t68;
t62 = t66 ^ 2;
t122 = t62 * t68;
t121 = t66 * t68;
t120 = t66 * t69;
t119 = t69 * t60;
t54 = t67 ^ 2 + t62;
t81 = -t22 / 0.2e1 + t126;
t93 = -t67 * pkin(4) / 0.2e1;
t2 = (t93 + t81) * t120;
t117 = t2 * qJD(1);
t4 = t128 * t121;
t116 = t4 * qJD(1);
t5 = (t93 - t81) * t68;
t115 = t5 * qJD(1);
t26 = -t67 * t119 - t68 * t70;
t24 = -t68 * t118 - t26;
t8 = (t22 * t69 + t24 * t68) * t66;
t114 = t8 * qJD(1);
t95 = t69 * t122;
t37 = (pkin(4) * t68 + t60) * t66;
t96 = t37 * t120;
t9 = -pkin(4) * t95 - t24 * t67 - t96;
t113 = t9 * qJD(1);
t112 = qJD(1) * t68;
t111 = qJD(1) * t69;
t110 = qJD(3) * t67;
t109 = qJD(4) * t68;
t108 = qJD(4) * t69;
t107 = qJD(5) * t66;
t10 = t65 * t62 * pkin(4) - t37 * t121 - t23 * t67;
t106 = t10 * qJD(1);
t25 = -t34 + t94;
t13 = -t60 * t122 - t25 * t67;
t105 = t13 * qJD(1);
t14 = -t62 * t119 + t26 * t67;
t104 = t14 * qJD(1);
t103 = t24 * qJD(4);
t80 = t64 / 0.2e1 + t65 / 0.2e1;
t27 = (-0.1e1 / 0.2e1 + t80) * t67 * t66;
t102 = t27 * qJD(1);
t31 = t54 * t60;
t101 = t31 * qJD(1);
t35 = (0.1e1 / 0.2e1 + t80) * t66;
t100 = t35 * qJD(1);
t42 = (t64 - t65) * t62;
t99 = t42 * qJD(1);
t43 = t54 * t68;
t39 = t43 * qJD(1);
t44 = t127 * t62;
t98 = t44 * qJD(1);
t45 = t54 * t69;
t41 = t45 * qJD(1);
t97 = t54 * qJD(1);
t58 = t66 * t109;
t92 = t67 * t109;
t91 = t66 * t108;
t90 = t67 * t108;
t89 = t66 * t112;
t88 = t68 * t107;
t87 = t67 * t112;
t86 = t68 * t110;
t85 = t66 * t111;
t84 = t69 * t107;
t83 = t67 * t111;
t82 = t69 * t110;
t79 = pkin(4) * t91;
t78 = pkin(4) * t85;
t77 = t66 * t87;
t76 = qJD(1) * t95;
t75 = t66 * t83;
t74 = qJD(1) * t67 - qJD(4);
t3 = pkin(4) * t96 + t128 * t24;
t72 = t3 * qJD(1) + t2 * qJD(2);
t7 = t37 * t66 + (t24 * t69 - t124) * t67;
t71 = -t7 * qJD(1) - t27 * qJD(2);
t47 = t74 * t69;
t46 = t74 * t68;
t40 = t45 * qJD(3);
t38 = t43 * qJD(3);
t36 = (0.1e1 / 0.2e1 - t127 / 0.2e1) * t66;
t33 = t66 * t47;
t32 = -t58 + t77;
t30 = t90 - t41;
t29 = t92 - t39;
t6 = -t124 / 0.2e1 + (t126 + t93) * t68;
t1 = t27 * qJD(3) + t2 * qJD(4);
t11 = [0, 0, 0, 0, 0, t54 * qJD(3), t31 * qJD(3), -qJD(4) * t95, t42 * qJD(4), t66 * t92, t66 * t90, 0, -t14 * qJD(4) + t38, t13 * qJD(4) + t40, -t9 * qJD(4) + t67 * t84 + t38, t10 * qJD(4) - t67 * t88 + t40, -t4 * qJD(4) + t44 * qJD(5), t7 * qJD(3) + t3 * qJD(4) - t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, t97, t101, 0, 0, 0, 0, 0, t39, t41, t39, t41, 0, t6 * qJD(4) + t36 * qJD(5) - t71; 0, 0, 0, 0, 0, 0, 0, -t76, t99, t32, t33, 0, t26 * qJD(4) - t104, t25 * qJD(4) + t105, -t103 - t113, t23 * qJD(4) + t106, pkin(4) * t58 - t116, -pkin(4) * t103 + t6 * qJD(3) + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t77, t98, t36 * qJD(3) - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t58, -t91, t58, 0, -t79 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t97, -t101, 0, 0, 0, 0, 0, t29, t30, t29, t30, 0, -t5 * qJD(4) - t35 * qJD(5) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t47, t46, t47, 0, -pkin(4) * t109 - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100; 0, 0, 0, 0, 0, 0, 0, t76, -t99, -t77, -t75, 0, -t86 + t104, -t82 - t105, -t84 - t86 + t113, -t82 + t88 - t106, t116, -pkin(4) * t84 + t5 * qJD(3) - t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t83, -t87, -t83, 0, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t89, 0, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, -t98, t35 * qJD(3) + t114 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t89, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
