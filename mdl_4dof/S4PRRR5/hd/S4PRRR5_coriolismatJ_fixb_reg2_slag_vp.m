% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:49
% EndTime: 2019-12-31 16:33:49
% DurationCPUTime: 0.52s
% Computational Cost: add. (502->83), mult. (1297->125), div. (0->0), fcn. (1234->6), ass. (0->72)
t70 = sin(qJ(4));
t68 = t70 ^ 2;
t73 = cos(qJ(4));
t69 = t73 ^ 2;
t60 = t69 - t68;
t94 = qJD(2) + qJD(3);
t113 = t94 * t60;
t102 = t68 + t69;
t74 = cos(qJ(3));
t106 = t74 * pkin(2);
t64 = -pkin(3) - t106;
t109 = t64 / 0.2e1;
t90 = -t106 / 0.2e1;
t112 = t109 - pkin(3) / 0.2e1 + t90;
t110 = pkin(3) / 0.2e1;
t105 = cos(qJ(2));
t71 = sin(qJ(3));
t72 = sin(qJ(2));
t53 = t105 * t71 + t74 * t72;
t108 = t53 * pkin(3);
t107 = t71 * pkin(2);
t52 = -t105 * t74 + t71 * t72;
t104 = t52 * t71;
t103 = t53 * t64;
t101 = pkin(2) * qJD(3);
t100 = pkin(3) * qJD(3);
t99 = qJD(2) * pkin(2);
t5 = (0.1e1 - t102) * t53 * t52;
t98 = t5 * qJD(1);
t97 = qJD(2) * t64;
t85 = t102 * t74;
t51 = pkin(2) * t85;
t96 = t51 * qJD(2);
t95 = t70 * qJD(4);
t67 = t73 * qJD(4);
t93 = t71 * t101;
t92 = t71 * t99;
t89 = t70 * t97;
t88 = t73 * t97;
t87 = t69 / 0.2e1 + t68 / 0.2e1;
t86 = t102 * t52;
t84 = pkin(2) * t94;
t83 = t70 * t92;
t27 = t94 * t53;
t26 = t94 * t52;
t82 = t94 * t70;
t63 = pkin(6) + t107;
t81 = t87 * t63;
t80 = t87 * t74;
t79 = t71 * t84;
t19 = (t63 * t85 + t64 * t71) * pkin(2);
t2 = (pkin(2) * t80 + t109 + t110) * t53 + (t107 / 0.2e1 - t81 + t87 * pkin(6)) * t52;
t78 = -t2 * qJD(1) - t19 * qJD(2);
t77 = t90 + t110 - t64 / 0.2e1;
t30 = t77 * t70;
t76 = qJD(2) * t30 + t100 * t70;
t31 = t77 * t73;
t75 = qJD(2) * t31 + t100 * t73;
t61 = t70 * t67;
t59 = t70 * t93;
t56 = t60 * qJD(4);
t50 = t73 * t82;
t48 = t51 * qJD(3);
t33 = t112 * t73;
t32 = t112 * t70;
t24 = t52 * t73;
t23 = t52 * t70;
t8 = t102 * t26;
t7 = t23 * qJD(4) - t27 * t73;
t6 = t24 * qJD(4) + t53 * t82;
t1 = t103 / 0.2e1 - t108 / 0.2e1 + (t104 / 0.2e1 + t53 * t80) * pkin(2) + (-t81 - t102 * pkin(6) / 0.2e1) * t52;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * qJD(2), -t105 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, (-t53 * t74 - t104) * t99, 0, 0, 0, 0, 0, 0, t7, t6, -t8, t98 + (-t63 * t86 + t103) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6, -t8, t98 + t1 * qJD(2) + (-pkin(6) * t86 - t108) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t94 - t53 * t67, t24 * t94 + t53 * t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t2 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t74 * t101, 0, 0, t61, t56, 0, -t61, 0, 0, t64 * t95 - t73 * t93, t64 * t67 + t59, t48, t19 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t74 * t84, 0, 0, t61, t56, 0, -t61, 0, 0, t32 * qJD(4) - t73 * t79, qJD(4) * t33 + t59 + t83, t48 + t96, (-pkin(3) * t71 + pkin(6) * t85) * t101 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t113, t67, -t50, -t95, 0, qJD(3) * t32 - t63 * t67 + t89, qJD(3) * t33 + t63 * t95 + t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t2 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t74 * t99, 0, 0, t61, t56, 0, -t61, 0, 0, -qJD(4) * t30 + t73 * t92, -qJD(4) * t31 - t83, -t96, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t56, 0, -t61, 0, 0, -pkin(3) * t95, -pkin(3) * t67, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t113, t67, -t50, -t95, 0, -pkin(6) * t67 - t76, pkin(6) * t95 - t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t113, 0, t50, 0, 0, qJD(3) * t30 - t89, qJD(3) * t31 - t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t113, 0, t50, 0, 0, t76, t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
