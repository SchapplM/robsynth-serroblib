% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:00
% EndTime: 2019-12-05 18:00:04
% DurationCPUTime: 0.94s
% Computational Cost: add. (1801->126), mult. (2960->163), div. (0->0), fcn. (3119->4), ass. (0->100)
t99 = -pkin(1) - pkin(6);
t170 = -pkin(7) + t99;
t98 = cos(qJ(3));
t88 = t170 * t98;
t95 = sin(qJ(4));
t153 = t95 * t88;
t96 = sin(qJ(3));
t87 = t170 * t96;
t97 = cos(qJ(4));
t80 = t97 * t87;
t167 = -t80 - t153;
t83 = t95 * t98 + t97 * t96;
t44 = t83 * qJ(5) + t167;
t156 = t44 * pkin(4);
t122 = qJD(3) + qJD(4);
t132 = t83 * qJD(4);
t155 = t97 * pkin(3);
t92 = pkin(4) + t155;
t161 = -t92 / 0.2e1;
t107 = t155 / 0.2e1 + t161;
t158 = pkin(4) * t83;
t30 = t158 / 0.2e1 - t107 * t83;
t174 = pkin(4) * t132 + t30 * qJD(3);
t159 = pkin(3) * t95;
t85 = -t95 * t96 + t97 * t98;
t173 = t30 * qJD(4) - (t85 * t159 - t92 * t83) * qJD(3);
t168 = (pkin(4) / 0.2e1 + t107) * t83;
t172 = qJD(3) * t168;
t171 = qJD(4) * t168;
t40 = t95 * t87 - t97 * t88;
t169 = -t85 * qJ(5) - t40;
t166 = t122 * t40;
t82 = t83 ^ 2;
t164 = t85 ^ 2;
t162 = -t44 / 0.2e1;
t114 = -t80 / 0.2e1;
t6 = (t162 + t44 / 0.2e1) * t85;
t160 = t6 * qJD(3);
t157 = pkin(4) * t85;
t154 = t98 * pkin(3);
t151 = pkin(3) * qJD(4);
t150 = qJD(3) * pkin(3);
t148 = t6 * qJD(1);
t12 = -t169 * t85 + t44 * t83;
t145 = t12 * qJD(1);
t115 = -t95 * t83 / 0.2e1;
t25 = (t161 - pkin(4) / 0.2e1) * t85 + (t115 - t98 / 0.2e1) * pkin(3);
t144 = t25 * qJD(1);
t109 = -t82 / 0.2e1 - t164 / 0.2e1;
t27 = -0.1e1 / 0.2e1 + t109;
t143 = t27 * qJD(1);
t142 = t168 * qJD(1);
t49 = t82 - t164;
t139 = t49 * qJD(1);
t90 = t96 * pkin(3) + qJ(2);
t50 = t83 * t154 + t90 * t85;
t138 = t50 * qJD(1);
t51 = t85 * t154 - t90 * t83;
t137 = t51 * qJD(1);
t53 = t114 + t80 / 0.2e1;
t136 = t53 * qJD(1);
t56 = t82 + t164;
t135 = t56 * qJD(1);
t67 = t90 + t158;
t134 = t67 * qJD(1);
t133 = t83 * qJD(1);
t131 = t85 * qJD(1);
t130 = t85 * qJD(4);
t89 = t96 ^ 2 - t98 ^ 2;
t129 = t89 * qJD(1);
t128 = t96 * qJD(1);
t127 = t96 * qJD(3);
t126 = t98 * qJD(1);
t125 = t98 * qJD(3);
t124 = qJ(2) * qJD(3);
t123 = qJD(1) * qJ(2);
t120 = pkin(4) * t131;
t119 = t95 * t151;
t118 = t90 * t133;
t117 = t90 * t131;
t116 = t96 * t126;
t112 = t96 * t123;
t111 = t98 * t123;
t110 = pkin(3) * t122;
t55 = t122 * t85;
t10 = t67 * (t154 + t157);
t105 = t10 * qJD(1) + t6 * qJD(2);
t11 = t67 * t157;
t104 = t11 * qJD(1);
t100 = t162 * t155 - t161 * t44;
t2 = -t156 / 0.2e1 + t100;
t71 = (-t92 + t155) * t159;
t101 = -t2 * qJD(1) - qJD(2) * t168 - t71 * qJD(3);
t63 = t83 * t131;
t54 = t122 * t83;
t41 = 0.2e1 * t114 - t153;
t26 = 0.1e1 / 0.2e1 + t109;
t24 = pkin(3) * t115 + t85 * t161 + t154 / 0.2e1 + t157 / 0.2e1;
t1 = t156 / 0.2e1 + t100;
t3 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t96 * t125, t89 * qJD(3), 0, 0, 0, qJD(2) * t96 + t98 * t124, qJD(2) * t98 - t96 * t124, -t83 * t55, t122 * t49, 0, 0, 0, t83 * qJD(2) + t50 * qJD(3) + t90 * t130, t85 * qJD(2) + t51 * qJD(3) - t90 * t132, t56 * qJD(5), t67 * qJD(2) + t10 * qJD(3) + t11 * qJD(4) + t12 * qJD(5); 0, 0, 0, 0, qJD(1), t123, 0, 0, 0, 0, 0, t128, t126, 0, 0, 0, 0, 0, t133, t131, 0, t26 * qJD(5) + t134 + t160; 0, 0, 0, 0, 0, 0, -t116, t129, -t127, -t125, 0, -t99 * t127 + t111, -t99 * t125 - t112, -t63, t139, -t54, -t55, 0, t167 * qJD(3) + t41 * qJD(4) + t138, t137 + t166, t173, (t159 * t169 + t44 * t92) * qJD(3) + t1 * qJD(4) + t24 * qJD(5) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t139, -t54, -t55, 0, t41 * qJD(3) + t167 * qJD(4) + t117, -t118 + t166, t174, t1 * qJD(3) + qJD(4) * t156 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t26 * qJD(2) + t24 * qJD(3) + t145; 0, 0, 0, 0, -qJD(1), -t123, 0, 0, 0, 0, 0, -t128, -t126, 0, 0, 0, 0, 0, -t133, -t131, 0, t27 * qJD(5) - t134 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t125, 0, 0, 0, 0, 0, -t54, -t55, 0, t148 - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, -t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143; 0, 0, 0, 0, 0, 0, t116, -t129, 0, 0, 0, -t111, t112, t63, -t139, 0, 0, 0, t53 * qJD(4) - t138, -t137, -t171, t2 * qJD(4) + t25 * qJD(5) - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t97 * t151, 0, t71 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95 * t110 + t136, -t97 * t110, -t142, -pkin(4) * t119 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t139, 0, 0, 0, -t53 * qJD(3) - t117, t118, t172, -t2 * qJD(3) - qJD(5) * t157 - t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t150 - t136, t97 * t150, t142, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, pkin(4) * t130 - t27 * qJD(2) - t25 * qJD(3) - t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
