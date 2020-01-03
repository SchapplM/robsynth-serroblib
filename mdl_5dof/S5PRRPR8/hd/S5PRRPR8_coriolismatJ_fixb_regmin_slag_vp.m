% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:44
% EndTime: 2019-12-31 17:42:45
% DurationCPUTime: 0.48s
% Computational Cost: add. (710->86), mult. (1582->136), div. (0->0), fcn. (1790->8), ass. (0->84)
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t60 = -t69 ^ 2 + t72 ^ 2;
t90 = qJD(2) + qJD(3);
t119 = t90 * t60;
t67 = sin(pkin(9));
t70 = sin(qJ(3));
t103 = t67 * t70;
t73 = cos(qJ(3));
t65 = t73 * pkin(2) + pkin(3);
t68 = cos(pkin(9));
t47 = -pkin(2) * t103 + t68 * t65;
t45 = -pkin(4) - t47;
t63 = -t68 * pkin(3) - pkin(4);
t118 = t45 + t63;
t53 = (t68 * t73 - t103) * pkin(2);
t113 = -t53 / 0.2e1;
t82 = t113 - t63 / 0.2e1 - t45 / 0.2e1;
t22 = t82 * t72;
t71 = sin(qJ(2));
t74 = cos(qJ(2));
t55 = t70 * t74 + t73 * t71;
t104 = t67 * t55;
t54 = -t70 * t71 + t73 * t74;
t49 = t68 * t54;
t86 = -t49 + t104;
t79 = t86 / 0.2e1 - t104 / 0.2e1 + t49 / 0.2e1;
t15 = t79 * t72;
t91 = t15 * qJD(1);
t94 = qJD(3) * t72;
t117 = t22 * qJD(2) - t63 * t94 + t91;
t21 = t82 * t69;
t14 = t79 * t69;
t92 = t14 * qJD(1);
t95 = qJD(3) * t69;
t116 = t21 * qJD(2) - t63 * t95 + t92;
t96 = qJD(2) * t72;
t115 = -t45 * t96 + t91;
t97 = qJD(2) * t69;
t114 = -t45 * t97 + t92;
t112 = -t69 / 0.2e1;
t111 = t69 / 0.2e1;
t110 = -t72 / 0.2e1;
t109 = t72 / 0.2e1;
t81 = t67 * t54 + t68 * t55;
t108 = t81 * t47;
t107 = t81 * t68;
t102 = t68 * t70;
t48 = pkin(2) * t102 + t67 * t65;
t106 = t86 * t48;
t105 = t86 * t67;
t101 = pkin(2) * qJD(3);
t100 = pkin(3) * qJD(3);
t99 = qJD(2) * pkin(2);
t93 = qJD(5) * t69;
t66 = qJD(5) * t72;
t52 = (t67 * t73 + t102) * pkin(2);
t87 = t52 * t97;
t85 = pkin(2) * t90;
t84 = t90 * t72;
t83 = t90 * t69;
t75 = t106 / 0.2e1 + t81 * t113 + t108 / 0.2e1 - t86 * t52 / 0.2e1;
t76 = (-t105 / 0.2e1 - t107 / 0.2e1) * pkin(3);
t1 = t76 + t75;
t18 = -t47 * t52 + t48 * t53;
t80 = t1 * qJD(1) - t18 * qJD(2);
t62 = t67 * pkin(3) + pkin(7);
t61 = t69 * t66;
t59 = t60 * qJD(5);
t51 = t72 * t83;
t46 = pkin(7) + t48;
t42 = t52 * t95;
t37 = t90 * t54;
t36 = t90 * t55;
t24 = t118 * t109 + t53 * t110;
t23 = t118 * t111 + t53 * t112;
t17 = (t109 - t110) * t86;
t16 = (t111 - t112) * t86;
t12 = t15 * qJD(5);
t10 = t14 * qJD(5);
t4 = t16 * qJD(5) - t81 * t84;
t3 = t17 * qJD(5) + t81 * t83;
t2 = t76 - t75;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t71, -qJD(2) * t74, 0, -t36, -t37, (-t106 - t108) * qJD(2) + t2 * qJD(3), 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, -t36, -t37, t2 * qJD(2) + (-t105 - t107) * t100, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 * t16 - t66 * t81, t90 * t17 + t81 * t93; 0, 0, 0, 0, 0, 0, 0, -t1 * qJD(3), 0, 0, 0, 0, 0, -t10, -t12; 0, 0, 0, 0, 0, -t70 * t101, -t73 * t101, t18 * qJD(3), t61, t59, 0, 0, 0, t45 * t93 - t52 * t94, t45 * t66 + t42; 0, 0, 0, 0, 0, -t70 * t85, -t73 * t85, (-t52 * t68 + t53 * t67) * t100 - t80, t61, t59, 0, 0, 0, t23 * qJD(5) - t52 * t84, t24 * qJD(5) + t42 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, t119, t66, -t93, 0, t23 * qJD(3) - t46 * t66 - t114, t24 * qJD(3) + t46 * t93 - t115; 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2), 0, 0, 0, 0, 0, -t10, -t12; 0, 0, 0, 0, 0, t70 * t99, t73 * t99, t80, t61, t59, 0, 0, 0, -t21 * qJD(5) + t52 * t96, -t22 * qJD(5) - t87; 0, 0, 0, 0, 0, 0, 0, 0, t61, t59, 0, 0, 0, t63 * t93, t63 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, t119, t66, -t93, 0, -t62 * t66 - t116, t62 * t93 - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 * t14, t90 * t15; 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t119, 0, 0, 0, t21 * qJD(3) + t114, t22 * qJD(3) + t115; 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t119, 0, 0, 0, t116, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
