% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:55
% EndTime: 2019-12-05 16:06:59
% DurationCPUTime: 0.63s
% Computational Cost: add. (1059->86), mult. (2192->128), div. (0->0), fcn. (2300->4), ass. (0->69)
t113 = cos(pkin(8));
t77 = sin(pkin(8));
t78 = sin(qJ(3));
t79 = cos(qJ(3));
t64 = t113 * t78 + t77 * t79;
t61 = t64 ^ 2;
t62 = -t113 * t79 + t77 * t78;
t34 = t62 ^ 2 + t61;
t136 = t34 * qJD(2);
t135 = t34 * qJD(4);
t134 = qJ(4) + pkin(6);
t68 = t134 * t78;
t69 = t134 * t79;
t43 = t113 * t68 + t77 * t69;
t120 = t77 * t68;
t67 = t113 * t69;
t91 = t67 - t120;
t87 = t43 * t64 - t62 * t91;
t131 = qJD(4) * t87;
t127 = t87 * qJD(2);
t89 = t67 / 0.2e1;
t125 = t64 * pkin(4);
t124 = t78 * pkin(3);
t123 = t62 * t77;
t71 = t77 * pkin(3) + qJ(5);
t122 = t71 * t62;
t73 = -t113 * pkin(3) - pkin(4);
t121 = t73 * t64;
t119 = qJD(3) * pkin(3);
t115 = t62 * qJ(5);
t92 = -t79 * pkin(3) - pkin(2);
t32 = t62 * pkin(4) - t64 * qJ(5) + t92;
t33 = t115 + t124 + t125;
t9 = t32 * t64 + t33 * t62;
t114 = t9 * qJD(2);
t112 = qJD(2) * t79;
t10 = t32 * t62 - t33 * t64;
t111 = t10 * qJD(2);
t74 = t124 / 0.2e1;
t17 = t74 + (pkin(4) / 0.2e1 - t73 / 0.2e1) * t64 + (qJ(5) / 0.2e1 + t71 / 0.2e1) * t62;
t107 = t17 * qJD(2);
t90 = t113 * t64;
t82 = -t123 / 0.2e1 - t90 / 0.2e1;
t24 = (-t78 / 0.2e1 + t82) * pkin(3);
t106 = t24 * qJD(2);
t103 = t61 * qJD(2);
t102 = t62 * qJD(2);
t55 = t62 * qJD(3);
t101 = t64 * qJD(2);
t56 = t64 * qJD(3);
t100 = t64 * qJD(5);
t70 = -t78 ^ 2 + t79 ^ 2;
t99 = t70 * qJD(2);
t98 = t78 * qJD(3);
t97 = t79 * qJD(3);
t96 = pkin(2) * t78 * qJD(2);
t95 = pkin(2) * t112;
t94 = t62 * t101;
t93 = t78 * t112;
t3 = t32 * t33;
t86 = t3 * qJD(2);
t8 = t92 * t124;
t85 = t8 * qJD(2);
t41 = t89 - t67 / 0.2e1;
t83 = t41 * qJD(2) + t71 * qJD(3);
t31 = 0.2e1 * t89 - t120;
t23 = t82 * pkin(3) + t74;
t18 = -t122 / 0.2e1 + t121 / 0.2e1 + t74 + t115 / 0.2e1 + t125 / 0.2e1;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, 0, (-t90 - t123) * t119, -t56, 0, -t55, (t121 - t122) * qJD(3) + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t78 * t97, t70 * qJD(3), 0, 0, 0, -pkin(2) * t98, -pkin(2) * t97, t135, t8 * qJD(3) + t131, t9 * qJD(3) - t62 * t100, t135, t10 * qJD(3) + t61 * qJD(5), t3 * qJD(3) - t32 * t100 + t131; 0, 0, 0, 0, t93, t99, t97, -t98, 0, -pkin(6) * t97 - t96, pkin(6) * t98 - t95, (t113 * t62 - t64 * t77) * t119, (-t113 * t91 - t43 * t77) * t119 + t23 * qJD(4) + t85, -qJD(3) * t91 + t114, (-t73 * t62 - t71 * t64) * qJD(3) - qJD(5) * t62, -qJD(3) * t43 + t111, (-t43 * t71 + t73 * t91) * qJD(3) + t18 * qJD(4) + t31 * qJD(5) + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t23 * qJD(3) + t127, 0, t136, 0, t18 * qJD(3) + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t55, t103, t31 * qJD(3) - t32 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t93, -t99, 0, 0, 0, t96, t95, 0, t24 * qJD(4) - t85, -t64 * qJD(4) - t114, 0, -t62 * qJD(4) - t111, -t17 * qJD(4) + t41 * qJD(5) - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t71 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t101, 0, -t102, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t24 * qJD(3) - t127, t56, -t136, t55, t17 * qJD(3) - t100 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, t101, 0, t102, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t103, -t41 * qJD(3) + (qJD(2) * t32 + qJD(4)) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
