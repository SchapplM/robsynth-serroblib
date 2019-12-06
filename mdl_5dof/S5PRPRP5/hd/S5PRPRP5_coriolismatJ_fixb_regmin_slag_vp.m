% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:44
% EndTime: 2019-12-05 15:38:47
% DurationCPUTime: 0.81s
% Computational Cost: add. (882->119), mult. (2160->179), div. (0->0), fcn. (2321->6), ass. (0->105)
t100 = cos(qJ(2));
t152 = cos(qJ(4));
t96 = sin(pkin(8));
t97 = cos(pkin(8));
t98 = sin(qJ(4));
t101 = t152 * t97 - t98 * t96;
t154 = t101 * t100;
t118 = -t154 / 0.2e1;
t75 = t152 * t96 + t98 * t97;
t52 = t100 * t75;
t116 = t52 / 0.2e1;
t94 = t96 ^ 2;
t95 = t97 ^ 2;
t87 = t94 + t95;
t71 = t75 ^ 2;
t99 = sin(qJ(2));
t153 = t99 / 0.2e1;
t151 = pkin(6) + qJ(3);
t113 = t100 * t152;
t110 = -t113 / 0.2e1;
t143 = t100 * t98;
t115 = -t143 / 0.2e1;
t145 = t96 * t110 + t97 * t115;
t27 = t116 + t145;
t137 = t27 * qJD(1);
t150 = -t75 * qJD(3) + t137;
t112 = t151 * t96;
t86 = t151 * t97;
t47 = -t98 * t112 + t152 * t86;
t44 = t47 * qJD(4);
t149 = -t137 - t44;
t114 = t143 / 0.2e1;
t148 = t97 * t110 + t96 * t114;
t109 = t113 / 0.2e1;
t147 = t97 * t109 + t96 * t115;
t146 = t96 * t109 + t97 * t114;
t108 = -pkin(4) * t101 - t75 * qJ(5);
t91 = -t97 * pkin(3) - pkin(2);
t34 = t108 + t91;
t43 = pkin(4) * t75 - qJ(5) * t101;
t7 = -t101 * t43 + t34 * t75;
t142 = t7 * qJD(2);
t51 = t75 * t99;
t53 = t101 * t99;
t92 = t100 * t99;
t10 = t154 * t53 + t51 * t52 - t92;
t141 = t10 * qJD(1);
t140 = t43 * qJD(2);
t70 = t101 ^ 2;
t21 = t70 - t71;
t139 = t21 * qJD(2);
t117 = -t52 / 0.2e1;
t26 = t117 + t146;
t138 = t26 * qJD(1);
t119 = t154 / 0.2e1;
t29 = t119 + t148;
t136 = t29 * qJD(1);
t45 = t70 + t71;
t135 = t45 * qJD(2);
t46 = t152 * t112 + t98 * t86;
t134 = t46 * qJD(4);
t111 = t87 * t100;
t48 = t111 * t99 - t92;
t133 = t48 * qJD(1);
t132 = t51 * qJD(4);
t131 = t53 * qJD(4);
t130 = t71 * qJD(2);
t129 = t101 * qJD(2);
t64 = t101 * qJD(4);
t128 = t75 * qJD(2);
t66 = t75 * qJD(4);
t127 = t75 * qJD(5);
t126 = t87 * qJD(2);
t125 = t99 * qJD(2);
t124 = qJD(4) * qJ(5);
t123 = t75 * t125;
t122 = t101 * t128;
t121 = t91 * t128;
t120 = -t100 * t43 / 0.2e1;
t102 = pkin(4) * t116 + qJ(5) * t118;
t2 = t120 + t102;
t3 = t34 * t43;
t107 = t2 * qJD(1) + t3 * qJD(2);
t103 = -t51 * t75 / 0.2e1 - t53 * t101 / 0.2e1;
t12 = t153 + t103;
t9 = t101 * t47 + t46 * t75;
t106 = -t12 * qJD(1) + t9 * qJD(2);
t28 = t118 + t147;
t8 = -t101 * t34 - t43 * t75;
t105 = t28 * qJD(1) - t8 * qJD(2);
t55 = (0.1e1 / 0.2e1 - t94 / 0.2e1 - t95 / 0.2e1) * t99;
t77 = t87 * qJ(3);
t104 = t55 * qJD(1) - t77 * qJD(2);
t56 = (0.1e1 + t87) * t153;
t33 = t116 + t146;
t32 = t117 + t145;
t31 = t118 + t148;
t30 = t119 + t147;
t24 = t27 * qJD(4);
t22 = t27 * qJD(2);
t15 = t32 * qJD(4) - t101 * t125;
t14 = t32 * qJD(2) - t131;
t13 = t153 - t103;
t1 = t120 - t102;
t4 = [0, 0, 0, 0, 0, 0, 0, t48 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(2); 0, 0, -t125, -t100 * qJD(2), -t97 * t125, t96 * t125, qJD(2) * t111, t133 + (-t99 * pkin(2) + qJ(3) * t111) * qJD(2) + t56 * qJD(3), 0, 0, 0, 0, 0, t15, t31 * qJD(4) + t123, t15, (t101 * t154 + t52 * t75) * qJD(2), t30 * qJD(4) - t123, t141 + (t154 * t47 + t99 * t34 + t52 * t46) * qJD(2) + t13 * qJD(3) + t1 * qJD(4) + t33 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t56 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t31 * qJD(2) + t132, t14, 0, t30 * qJD(2) - t132, t1 * qJD(2) + (-t53 * pkin(4) - t51 * qJ(5)) * qJD(4) + t53 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * qJD(2) + t131; 0, 0, 0, 0, 0, 0, 0, -t55 * qJD(3) - t133, 0, 0, 0, 0, 0, -t24, -t29 * qJD(4), -t24, 0, -t28 * qJD(4), -t12 * qJD(3) + t2 * qJD(4) - t26 * qJD(5) - t141; 0, 0, 0, 0, 0, 0, t87 * qJD(3), t77 * qJD(3), t101 * t66, t21 * qJD(4), 0, 0, 0, t91 * t66, t91 * t64, t7 * qJD(4) + t101 * t127, t45 * qJD(3), t8 * qJD(4) + t71 * qJD(5), t9 * qJD(3) + t3 * qJD(4) - t127 * t34; 0, 0, 0, 0, 0, 0, t126, -t104, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, t106; 0, 0, 0, 0, 0, 0, 0, 0, t122, t139, t64, -t66, 0, t121 + t149, t129 * t91 + t134 - t136, t142 + t149, qJD(4) * t108 + qJD(5) * t101, -t105 - t134, (-t47 * pkin(4) - t46 * qJ(5)) * qJD(4) + t47 * qJD(5) + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t64, t130, -t128 * t34 - t138 + t44; 0, 0, 0, 0, 0, 0, 0, t55 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * qJD(2); 0, 0, 0, 0, 0, 0, -t126, t104, 0, 0, 0, 0, 0, t66, t64, t66, -t135, -t64, qJD(4) * t43 - t106 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t129, t128, 0, -t129, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t29 * qJD(2), t22, 0, t28 * qJD(2), -t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t139, 0, 0, 0, -t121 + t150, t136 - (qJD(2) * t91 + qJD(3)) * t101, -t142 + t150, 0, qJD(3) * t101 + t105, -qJD(3) * t43 - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t129, -t128, 0, t129, -t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, 0, -t130, t138 + (qJD(2) * t34 + qJD(3)) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
