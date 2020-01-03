% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:15
% EndTime: 2019-12-31 18:09:17
% DurationCPUTime: 0.60s
% Computational Cost: add. (1324->88), mult. (2457->130), div. (0->0), fcn. (2565->6), ass. (0->71)
t116 = cos(pkin(8));
t79 = sin(pkin(8));
t81 = sin(qJ(3));
t82 = cos(qJ(3));
t68 = t116 * t81 + t79 * t82;
t65 = t68 ^ 2;
t66 = -t116 * t82 + t79 * t81;
t39 = t66 ^ 2 + t65;
t139 = t39 * qJD(1);
t138 = t39 * qJD(4);
t74 = sin(pkin(7)) * pkin(1) + pkin(6);
t137 = qJ(4) + t74;
t63 = t137 * t81;
t64 = t137 * t82;
t35 = t116 * t63 + t79 * t64;
t123 = t79 * t63;
t53 = t116 * t64;
t95 = t53 - t123;
t91 = t35 * t68 - t66 * t95;
t134 = qJD(4) * t91;
t130 = t91 * qJD(1);
t93 = t53 / 0.2e1;
t128 = t68 * pkin(4);
t127 = t81 * pkin(3);
t126 = t66 * t79;
t72 = t79 * pkin(3) + qJ(5);
t125 = t72 * t66;
t75 = -t116 * pkin(3) - pkin(4);
t124 = t75 * t68;
t122 = qJD(3) * pkin(3);
t118 = t66 * qJ(5);
t76 = -cos(pkin(7)) * pkin(1) - pkin(2);
t87 = -t82 * pkin(3) + t76;
t32 = t66 * pkin(4) - t68 * qJ(5) + t87;
t38 = t118 + t127 + t128;
t9 = t32 * t68 + t38 * t66;
t117 = t9 * qJD(1);
t115 = qJD(1) * t82;
t10 = t32 * t66 - t38 * t68;
t114 = t10 * qJD(1);
t77 = t127 / 0.2e1;
t17 = t77 + (pkin(4) / 0.2e1 - t75 / 0.2e1) * t68 + (qJ(5) / 0.2e1 + t72 / 0.2e1) * t66;
t110 = t17 * qJD(1);
t94 = t116 * t68;
t85 = -t126 / 0.2e1 - t94 / 0.2e1;
t31 = (-t81 / 0.2e1 + t85) * pkin(3);
t109 = t31 * qJD(1);
t106 = t65 * qJD(1);
t105 = t66 * qJD(1);
t59 = t66 * qJD(3);
t104 = t68 * qJD(1);
t60 = t68 * qJD(3);
t103 = t68 * qJD(5);
t71 = -t81 ^ 2 + t82 ^ 2;
t102 = t71 * qJD(1);
t101 = t81 * qJD(3);
t100 = t82 * qJD(3);
t99 = t66 * t104;
t98 = t76 * t81 * qJD(1);
t97 = t76 * t115;
t96 = t81 * t115;
t3 = t32 * t38;
t90 = t3 * qJD(1);
t8 = t87 * t127;
t89 = t8 * qJD(1);
t33 = t93 - t53 / 0.2e1;
t86 = t33 * qJD(1) + t72 * qJD(3);
t30 = t85 * pkin(3) + t77;
t19 = 0.2e1 * t93 - t123;
t18 = -t125 / 0.2e1 + t124 / 0.2e1 + t77 + t118 / 0.2e1 + t128 / 0.2e1;
t1 = [0, 0, 0, 0, t81 * t100, t71 * qJD(3), 0, 0, 0, t76 * t101, t76 * t100, t138, t8 * qJD(3) + t134, t9 * qJD(3) - t66 * t103, t138, t10 * qJD(3) + t65 * qJD(5), t3 * qJD(3) - t32 * t103 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t96, t102, t100, -t101, 0, -t74 * t100 + t98, t74 * t101 + t97, (t116 * t66 - t68 * t79) * t122, (-t116 * t95 - t35 * t79) * t122 + t30 * qJD(4) + t89, -qJD(3) * t95 + t117, (-t75 * t66 - t72 * t68) * qJD(3) - qJD(5) * t66, -qJD(3) * t35 + t114, (-t35 * t72 + t75 * t95) * qJD(3) + t18 * qJD(4) + t19 * qJD(5) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t30 * qJD(3) + t130, 0, t139, 0, t18 * qJD(3) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t59, t106, t19 * qJD(3) - t32 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t100, 0, (-t94 - t126) * t122, -t60, 0, -t59, (t124 - t125) * qJD(3) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, -t96, -t102, 0, 0, 0, -t98, -t97, 0, t31 * qJD(4) - t89, -t68 * qJD(4) - t117, 0, -t66 * qJD(4) - t114, -t17 * qJD(4) + t33 * qJD(5) - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t72 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t104, 0, -t105, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t31 * qJD(3) - t130, t60, -t139, t59, t17 * qJD(3) - t103 - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t104, 0, t105, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, -t106, -t33 * qJD(3) + (qJD(1) * t32 + qJD(4)) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
