% Calculate inertial parameters regressor of coriolis matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:48
% EndTime: 2019-12-31 17:53:50
% DurationCPUTime: 1.00s
% Computational Cost: add. (877->102), mult. (1805->113), div. (0->0), fcn. (1916->4), ass. (0->76)
t124 = sin(qJ(4));
t125 = cos(qJ(4));
t83 = sin(pkin(7));
t84 = cos(pkin(7));
t62 = t83 * t124 + t84 * t125;
t127 = t62 ^ 2;
t64 = -t84 * t124 + t83 * t125;
t60 = t64 ^ 2;
t20 = -t60 + t127;
t139 = t20 * qJD(1);
t138 = t20 * qJD(4);
t123 = -pkin(6) + qJ(2);
t71 = t123 * t84;
t92 = t123 * t83;
t38 = t124 * t71 - t125 * t92;
t39 = t124 * t92 + t125 * t71;
t87 = -t38 * t64 + t39 * t62;
t130 = t87 * qJD(1);
t126 = -t83 / 0.2e1;
t128 = -t125 * t64 / 0.2e1 - t62 * t124 / 0.2e1;
t85 = t126 + t128;
t137 = qJD(3) * t85 - t130;
t86 = t126 - t128;
t136 = qJD(3) * t86 + t130;
t135 = qJD(2) * t85;
t134 = qJD(2) * t86;
t131 = t85 * qJD(1);
t129 = t87 * qJD(2);
t81 = t83 ^ 2;
t75 = t84 ^ 2 + t81;
t69 = -t84 * pkin(2) - t83 * qJ(3) - pkin(1);
t51 = t84 * pkin(3) - t69;
t88 = t62 * pkin(4) - t64 * qJ(5);
t18 = t88 + t51;
t89 = t64 * pkin(4) + t62 * qJ(5);
t1 = t18 * t89;
t122 = t1 * qJD(1);
t5 = t18 * t64 + t62 * t89;
t121 = t5 * qJD(1);
t6 = t18 * t62 - t64 * t89;
t120 = t6 * qJD(1);
t16 = -t60 - t127;
t117 = t16 * qJD(1);
t116 = t89 * qJD(1);
t111 = t38 * qJD(4);
t35 = t39 * qJD(4);
t110 = t60 * qJD(1);
t109 = t62 * qJD(1);
t53 = t62 * qJD(4);
t108 = t64 * qJD(1);
t57 = t64 * qJD(4);
t107 = t64 * qJD(5);
t70 = t75 * qJ(2);
t106 = t70 * qJD(1);
t105 = t75 * qJD(1);
t104 = t81 * qJD(1);
t103 = t83 * qJD(1);
t102 = t83 * qJD(3);
t101 = qJD(4) * qJ(5);
t100 = t18 * t103;
t99 = t51 * t108;
t98 = t51 * t103;
t37 = t62 * t108;
t36 = t62 * t57;
t97 = t64 * t102;
t96 = t64 * t103;
t95 = t84 * t103;
t91 = t125 * qJD(4);
t90 = t124 * qJD(4);
t72 = t75 * qJD(2);
t66 = t70 * qJD(2);
t55 = t64 * qJD(2);
t41 = t62 * t102;
t40 = t62 * t103;
t13 = t16 * qJD(2);
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t66, 0, 0, 0, 0, 0, 0, t84 * t102, t72, t81 * qJD(3), -t69 * t102 + t66, -t36, t138, 0, t36, 0, 0, t51 * t57 + t41, -t51 * t53 + t97, t13, t51 * t102 + t129, -t36, 0, -t138, 0, 0, t36, t5 * qJD(4) - t62 * t107 + t41, t13, t6 * qJD(4) + t60 * qJD(5) - t97, t129 + t1 * qJD(4) + (t102 - t107) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t106, 0, 0, 0, 0, 0, 0, 0, t105, 0, t106, 0, 0, 0, 0, 0, 0, 0, 0, t117, t136, 0, 0, 0, 0, 0, 0, 0, t117, 0, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, t104, -t69 * t103, 0, 0, 0, 0, 0, 0, t40, t96, 0, t98 + t134, 0, 0, 0, 0, 0, 0, t40, 0, -t96, t100 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t139, -t53, t37, -t57, 0, -t35 + t99, -t51 * t109 + t111, 0, 0, -t37, -t53, -t139, 0, t57, t37, -t35 + t121, qJD(4) * t88 - t62 * qJD(5), -t111 + t120, t122 + (-t39 * pkin(4) - t38 * qJ(5)) * qJD(4) + t39 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t53, t110, -t18 * t108 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t106, 0, 0, 0, 0, 0, 0, 0, -t105, 0, -t106 - t102, 0, 0, 0, 0, 0, 0, -t57, t53, -t117, t137, 0, 0, 0, 0, 0, 0, -t57, -t117, -t53, -qJD(4) * t89 + t107 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t109, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t109, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, 0, -t104, (qJD(1) * t69 + qJD(2)) * t83, 0, 0, 0, 0, 0, 0, -t40, -t96, 0, -t98 - t135, 0, 0, 0, 0, 0, 0, -t40, 0, t96, -t100 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t91, 0, 0, 0, 0, 0, 0, 0, 0, -t90, 0, t91, (-t124 * pkin(4) + t125 * qJ(5)) * qJD(4) + t124 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t139, 0, -t37, 0, 0, t55 - t99, (qJD(1) * t51 - qJD(2)) * t62, 0, 0, t37, 0, t139, 0, 0, -t37, t55 - t121, 0, t62 * qJD(2) - t120, qJD(2) * t89 - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t109, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, t109, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t110, (qJD(1) * t18 - qJD(2)) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
