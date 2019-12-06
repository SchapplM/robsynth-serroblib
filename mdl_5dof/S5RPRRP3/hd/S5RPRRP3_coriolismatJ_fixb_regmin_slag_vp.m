% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:57
% EndTime: 2019-12-05 18:04:00
% DurationCPUTime: 0.90s
% Computational Cost: add. (1768->102), mult. (3202->148), div. (0->0), fcn. (3310->6), ass. (0->92)
t89 = sin(pkin(8)) * pkin(1) + pkin(6);
t153 = pkin(7) + t89;
t96 = sin(qJ(3));
t80 = t153 * t96;
t95 = sin(qJ(4));
t150 = t95 * t80;
t152 = cos(qJ(4));
t97 = cos(qJ(3));
t81 = t153 * t97;
t72 = t152 * t81;
t163 = -t72 + t150;
t83 = -t152 * t97 + t95 * t96;
t40 = t83 * qJ(5) + t163;
t155 = t40 * pkin(4);
t118 = t152 * pkin(3);
t92 = t118 + pkin(4);
t102 = t118 / 0.2e1 - t92 / 0.2e1;
t122 = qJD(3) + qJD(4);
t29 = t152 * t80 + t95 * t81;
t85 = t152 * t96 + t95 * t97;
t165 = -t85 * qJ(5) - t29;
t164 = t122 * t83;
t162 = t122 * t29;
t82 = t83 ^ 2;
t161 = t85 ^ 2;
t105 = -t72 / 0.2e1;
t158 = pkin(3) * t95;
t157 = pkin(4) * t83;
t156 = pkin(4) * t85;
t154 = t96 * pkin(3);
t151 = t92 * t85;
t149 = t95 * t83;
t142 = qJD(1) * t85;
t90 = -cos(pkin(8)) * pkin(1) - pkin(2);
t87 = -t97 * pkin(3) + t90;
t141 = qJD(1) * t87;
t140 = qJD(1) * t97;
t11 = -t165 * t85 + t40 * t83;
t139 = t11 * qJD(1);
t111 = -t149 / 0.2e1;
t112 = -t151 / 0.2e1;
t75 = -t156 / 0.2e1;
t25 = t112 + t75 + (t111 - t96 / 0.2e1) * pkin(3);
t136 = t25 * qJD(1);
t36 = (-pkin(4) / 0.2e1 - t102) * t83;
t133 = t36 * qJD(1);
t47 = t83 * t154 + t87 * t85;
t132 = t47 * qJD(1);
t48 = t85 * t154 - t87 * t83;
t131 = t48 * qJD(1);
t49 = t82 - t161;
t130 = t49 * qJD(1);
t51 = t105 + t72 / 0.2e1;
t129 = t51 * qJD(1);
t55 = t82 + t161;
t128 = t55 * qJD(1);
t127 = t83 * qJD(4);
t126 = t85 * qJD(4);
t88 = -t96 ^ 2 + t97 ^ 2;
t125 = t88 * qJD(1);
t124 = t96 * qJD(3);
t123 = t97 * qJD(3);
t121 = pkin(4) * t142;
t120 = pkin(4) * t126;
t119 = qJD(4) * t158;
t117 = t83 * t141;
t116 = t85 * t141;
t115 = t90 * t96 * qJD(1);
t114 = t90 * t140;
t113 = t96 * t140;
t108 = t152 * qJD(3);
t107 = t152 * qJD(4);
t54 = t122 * t85;
t61 = t87 + t157;
t9 = t61 * t156;
t104 = t9 * qJD(1);
t8 = t61 * (t154 + t156);
t103 = t8 * qJD(1);
t100 = t102 * t85;
t98 = t102 * t40;
t3 = -t155 / 0.2e1 - t98;
t76 = t156 / 0.2e1;
t38 = t76 + t100;
t73 = (t118 - t92) * t158;
t99 = -t3 * qJD(1) - t38 * qJD(2) - t73 * qJD(3);
t60 = t83 * t142;
t37 = t75 + t100;
t35 = t157 / 0.2e1 - t102 * t83;
t30 = 0.2e1 * t105 + t150;
t24 = pkin(3) * t111 + t112 + t154 / 0.2e1 + t76;
t2 = t155 / 0.2e1 - t98;
t1 = [0, 0, 0, 0, t96 * t123, t88 * qJD(3), 0, 0, 0, t90 * t124, t90 * t123, -t83 * t54, t122 * t49, 0, 0, 0, t47 * qJD(3) + t87 * t126, t48 * qJD(3) - t87 * t127, t55 * qJD(5), t8 * qJD(3) + t9 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t113, t125, t123, -t124, 0, -t89 * t123 + t115, t89 * t124 + t114, -t60, t130, -t164, -t54, 0, qJD(3) * t163 + t30 * qJD(4) + t132, t131 + t162, (-t85 * t158 + t92 * t83) * qJD(3) + t35 * qJD(4), (t158 * t165 + t40 * t92) * qJD(3) + t2 * qJD(4) + t24 * qJD(5) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t130, -t164, -t54, 0, t30 * qJD(3) + qJD(4) * t163 + t116, -t117 + t162, pkin(4) * t127 + t35 * qJD(3), t2 * qJD(3) + qJD(4) * t155 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t24 * qJD(3) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t123, 0, 0, 0, 0, 0, -t54, t164, 0, (-pkin(3) * t149 - t151) * qJD(3) + t37 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t164, 0, t37 * qJD(3) - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t113, -t125, 0, 0, 0, -t115, -t114, t60, -t130, 0, 0, 0, t51 * qJD(4) - t132, -t131, t36 * qJD(4), t3 * qJD(4) + t25 * qJD(5) - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -pkin(3) * t107, 0, t73 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t158 + t129, (-t108 - t107) * pkin(3), t133, -pkin(4) * t119 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t130, 0, 0, 0, -t51 * qJD(3) - t116, t117, -t36 * qJD(3), -t3 * qJD(3) - qJD(5) * t156 - t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t158 - t129, pkin(3) * t108, -t133, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t25 * qJD(3) + t120 - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
