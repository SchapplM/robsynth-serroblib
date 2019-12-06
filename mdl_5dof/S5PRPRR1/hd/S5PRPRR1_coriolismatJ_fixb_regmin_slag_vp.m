% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:08
% EndTime: 2019-12-05 15:43:11
% DurationCPUTime: 0.65s
% Computational Cost: add. (820->81), mult. (1705->108), div. (0->0), fcn. (1980->6), ass. (0->73)
t123 = qJD(4) + qJD(5);
t108 = cos(qJ(4));
t59 = sin(pkin(9));
t60 = cos(pkin(9));
t62 = sin(qJ(4));
t48 = -t108 * t60 + t59 * t62;
t61 = sin(qJ(5));
t103 = t61 * t48;
t107 = cos(qJ(5));
t50 = t108 * t59 + t60 * t62;
t44 = t107 * t50;
t111 = t44 - t103;
t74 = t107 * t48 + t50 * t61;
t116 = t74 * qJD(2);
t122 = t111 * t116;
t102 = pkin(6) + qJ(3);
t52 = t102 * t59;
t53 = t102 * t60;
t63 = -t108 * t53 + t52 * t62;
t26 = -pkin(7) * t48 - t63;
t73 = t108 * t52 + t53 * t62;
t64 = -pkin(7) * t50 - t73;
t121 = t123 * (-t107 * t26 - t61 * t64);
t114 = -t111 ^ 2 + t74 ^ 2;
t120 = t114 * qJD(2);
t112 = t107 * t64;
t69 = -t112 / 0.2e1;
t117 = qJD(3) * t74;
t104 = t61 * t26;
t115 = -t112 + t104;
t90 = t74 * qJD(5);
t66 = qJD(4) * t74 + t90;
t68 = t44 / 0.2e1;
t110 = pkin(4) * t50;
t109 = pkin(4) * t61;
t54 = t59 ^ 2 + t60 ^ 2;
t2 = t69 + t112 / 0.2e1;
t101 = t2 * qJD(2);
t56 = -pkin(3) * t60 - pkin(2);
t38 = pkin(4) * t48 + t56;
t6 = -t110 * t74 - t111 * t38;
t99 = t6 * qJD(2);
t7 = -t110 * t111 + t38 * t74;
t98 = t7 * qJD(2);
t96 = qJD(2) * t38;
t15 = 0.2e1 * t68 - t103;
t94 = t15 * qJD(2);
t27 = t48 ^ 2 - t50 ^ 2;
t93 = t27 * qJD(2);
t32 = t68 - t44 / 0.2e1;
t92 = t32 * qJD(2);
t30 = t32 * qJD(5);
t91 = t111 * qJD(2);
t87 = t111 * qJD(5);
t86 = t48 * qJD(2);
t45 = t48 * qJD(4);
t85 = t50 * qJD(2);
t46 = t50 * qJD(4);
t51 = t54 * qJ(3);
t84 = t51 * qJD(2);
t83 = t54 * qJD(2);
t80 = t74 * t96;
t79 = t111 * t96;
t78 = t48 * t85;
t72 = t107 * qJD(4);
t71 = t107 * qJD(5);
t70 = qJD(2) * t56 + qJD(3);
t67 = t32 * qJD(1);
t65 = qJD(4) * t111 + qJD(5) * t15;
t29 = t32 * qJD(4);
t9 = -qJD(4) * t15 - t87;
t3 = t104 + 0.2e1 * t69;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, -t65, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t54 * qJD(3), t51 * qJD(3), -t48 * t46, t27 * qJD(4), 0, 0, 0, t56 * t46, -t56 * t45, -t66 * t111, t123 * t114, 0, 0, 0, -qJD(4) * t6 + t38 * t87, -qJD(4) * t7 - t38 * t90; 0, 0, 0, 0, 0, 0, t83, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t78, t93, -t45, -t46, 0, qJD(4) * t63 + t56 * t85, qJD(4) * t73 - t56 * t86, -t122, t120, -t66, -t65, 0, -t99 + t121, qJD(4) * t115 + qJD(5) * t3 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t120, -t66, t9, 0, qJD(3) * t32 + t121 + t79, qJD(4) * t3 + qJD(5) * t115 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t83, -t84, 0, 0, 0, 0, 0, t46, -t45, 0, 0, 0, 0, 0, t65, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t86, 0, 0, 0, 0, 0, t91, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0; 0, 0, 0, 0, 0, 0, 0, 0, t78, -t93, 0, 0, 0, -t70 * t50, t70 * t48, t122, -t120, 0, -t30, 0, -qJD(3) * t111 + t99, qJD(5) * t2 + t117 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t86, 0, 0, 0, 0, 0, -t91, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t109, -pkin(4) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, 0, -t109 * t123 - t67, t101 + (-t72 - t71) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t120, 0, t29, 0, -qJD(3) * t15 - t79, -qJD(4) * t2 + t117 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, qJD(4) * t109 + t67, pkin(4) * t72 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
