% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:38
% EndTime: 2022-01-23 09:14:39
% DurationCPUTime: 0.65s
% Computational Cost: add. (960->83), mult. (1845->110), div. (0->0), fcn. (2120->8), ass. (0->74)
t125 = qJD(4) + qJD(5);
t109 = cos(qJ(4));
t60 = sin(pkin(9));
t61 = cos(pkin(9));
t64 = sin(qJ(4));
t51 = -t109 * t61 + t60 * t64;
t63 = sin(qJ(5));
t104 = t63 * t51;
t108 = cos(qJ(5));
t53 = t109 * t60 + t61 * t64;
t46 = t108 * t53;
t113 = t46 - t104;
t75 = t108 * t51 + t53 * t63;
t118 = t75 * qJD(1);
t124 = t113 * t118;
t56 = sin(pkin(8)) * pkin(1) + qJ(3);
t110 = pkin(6) + t56;
t47 = t110 * t60;
t48 = t110 * t61;
t65 = -t109 * t48 + t47 * t64;
t24 = -pkin(7) * t51 - t65;
t76 = t109 * t47 + t48 * t64;
t66 = -pkin(7) * t53 - t76;
t123 = t125 * (-t108 * t24 - t63 * t66);
t116 = -t113 ^ 2 + t75 ^ 2;
t122 = t116 * qJD(1);
t114 = t108 * t66;
t71 = -t114 / 0.2e1;
t119 = qJD(3) * t75;
t105 = t63 * t24;
t117 = -t114 + t105;
t92 = t75 * qJD(5);
t68 = qJD(4) * t75 + t92;
t70 = t46 / 0.2e1;
t112 = pkin(4) * t53;
t111 = pkin(4) * t63;
t55 = t60 ^ 2 + t61 ^ 2;
t2 = t71 + t114 / 0.2e1;
t103 = t2 * qJD(1);
t54 = -cos(pkin(8)) * pkin(1) - t61 * pkin(3) - pkin(2);
t38 = pkin(4) * t51 + t54;
t6 = -t112 * t75 - t113 * t38;
t101 = t6 * qJD(1);
t7 = -t112 * t113 + t38 * t75;
t100 = t7 * qJD(1);
t98 = qJD(1) * t38;
t26 = 0.2e1 * t70 - t104;
t96 = t26 * qJD(1);
t27 = t51 ^ 2 - t53 ^ 2;
t95 = t27 * qJD(1);
t32 = t70 - t46 / 0.2e1;
t94 = t32 * qJD(1);
t30 = t32 * qJD(5);
t93 = t113 * qJD(1);
t89 = t113 * qJD(5);
t39 = t55 * t56;
t88 = t39 * qJD(1);
t87 = t51 * qJD(1);
t49 = t51 * qJD(4);
t86 = t53 * qJD(1);
t50 = t53 * qJD(4);
t85 = t55 * qJD(1);
t82 = t75 * t98;
t81 = t113 * t98;
t80 = t51 * t86;
t74 = t108 * qJD(4);
t73 = t108 * qJD(5);
t72 = qJD(1) * t54 + qJD(3);
t69 = t32 * qJD(2);
t67 = qJD(4) * t113 + qJD(5) * t26;
t29 = t32 * qJD(4);
t9 = -qJD(4) * t26 - t89;
t3 = t105 + 0.2e1 * t71;
t1 = [0, 0, 0, 0, 0, t55 * qJD(3), t39 * qJD(3), -t51 * t50, t27 * qJD(4), 0, 0, 0, t54 * t50, -t54 * t49, -t68 * t113, t125 * t116, 0, 0, 0, -qJD(4) * t6 + t38 * t89, -qJD(4) * t7 - t38 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t85, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0; 0, 0, 0, 0, 0, 0, 0, -t80, t95, -t49, -t50, 0, qJD(4) * t65 + t54 * t86, qJD(4) * t76 - t54 * t87, -t124, t122, -t68, -t67, 0, -t101 + t123, qJD(4) * t117 + qJD(5) * t3 - t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t122, -t68, t9, 0, qJD(3) * t32 + t123 + t81, qJD(4) * t3 + qJD(5) * t117 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, 0, 0, 0, 0, -t67, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t68; 0, 0, 0, 0, 0, -t85, -t88, 0, 0, 0, 0, 0, t50, -t49, 0, 0, 0, 0, 0, t67, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t87, 0, 0, 0, 0, 0, t93, -t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t118; 0, 0, 0, 0, 0, 0, 0, t80, -t95, 0, 0, 0, -t72 * t53, t72 * t51, t124, -t122, 0, -t30, 0, -qJD(3) * t113 + t101, qJD(5) * t2 + t100 + t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t87, 0, 0, 0, 0, 0, -t93, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t111, -pkin(4) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, 0, -t111 * t125 - t69, t103 + (-t74 - t73) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t122, 0, t29, 0, -qJD(3) * t26 - t81, -qJD(4) * t2 + t119 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, qJD(4) * t111 + t69, pkin(4) * t74 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
