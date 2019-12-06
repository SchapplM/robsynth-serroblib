% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:45
% EndTime: 2019-12-05 15:54:50
% DurationCPUTime: 1.17s
% Computational Cost: add. (1091->106), mult. (2532->182), div. (0->0), fcn. (2889->8), ass. (0->107)
t194 = qJD(4) + qJD(5);
t100 = sin(qJ(5));
t173 = -t100 / 0.2e1;
t103 = cos(qJ(5));
t198 = -t103 / 0.2e1;
t105 = cos(qJ(2));
t104 = cos(qJ(4));
t98 = sin(pkin(9));
t158 = t104 * t98;
t101 = sin(qJ(4));
t99 = cos(pkin(9));
t162 = t101 * t99;
t81 = t158 + t162;
t65 = t105 * t81;
t157 = t104 * t99;
t163 = t101 * t98;
t175 = t157 - t163;
t67 = t105 * t175;
t111 = -t65 * t173 + t67 * t198;
t50 = t100 * t81 - t103 * t175;
t189 = t105 * t50;
t191 = -t189 / 0.2e1 + t111;
t204 = qJD(1) * t191;
t110 = t67 * t173 + t65 * t198;
t164 = t100 * t175;
t76 = t103 * t81;
t180 = t76 + t164;
t182 = t105 * t180;
t192 = t182 / 0.2e1 + t110;
t203 = qJD(1) * t192;
t202 = t191 * qJD(2);
t201 = t192 * qJD(2);
t190 = t189 / 0.2e1 + t111;
t102 = sin(qJ(2));
t64 = t81 * t102;
t66 = t175 * t102;
t200 = qJD(2) * t190 + t194 * (t100 * t66 + t103 * t64);
t193 = -t182 / 0.2e1 + t110;
t199 = qJD(2) * t193 + t194 * (t100 * t64 - t103 * t66);
t136 = t50 * qJD(5);
t9 = -t50 * qJD(4) - t136;
t171 = pkin(6) + qJ(3);
t87 = t171 * t98;
t88 = t171 * t99;
t120 = t101 * t88 + t104 * t87;
t112 = -t81 * pkin(7) - t120;
t117 = t101 * t87 - t104 * t88;
t36 = pkin(7) * t175 - t117;
t197 = t194 * (-t100 * t112 - t103 * t36);
t196 = t194 * (t100 * t36 - t103 * t112);
t181 = -t180 ^ 2 + t50 ^ 2;
t195 = t181 * qJD(2);
t184 = qJD(2) * t50;
t183 = qJD(3) * t50;
t137 = t180 * qJD(2);
t96 = t98 ^ 2;
t97 = t99 ^ 2;
t90 = t97 + t96;
t122 = t76 / 0.2e1;
t174 = pkin(4) * t81;
t170 = pkin(4) * qJD(5);
t169 = qJD(4) * pkin(4);
t94 = -t99 * pkin(3) - pkin(2);
t63 = -pkin(4) * t175 + t94;
t149 = qJD(2) * t63;
t23 = 0.2e1 * t122 + t164;
t143 = t23 * qJD(2);
t37 = t175 ^ 2 - t81 ^ 2;
t142 = t37 * qJD(2);
t109 = -t162 / 0.2e1 - t158 / 0.2e1;
t43 = (t81 / 0.2e1 + t109) * t105;
t141 = t43 * qJD(1);
t108 = -t157 / 0.2e1 + t163 / 0.2e1;
t44 = (t175 / 0.2e1 + t108) * t105;
t140 = t44 * qJD(1);
t48 = t122 - t76 / 0.2e1;
t139 = t48 * qJD(2);
t138 = t48 * qJD(5);
t133 = t180 * qJD(5);
t54 = (-0.1e1 + t90) * t105 * t102;
t132 = t54 * qJD(1);
t131 = t175 * qJD(2);
t77 = t175 * qJD(4);
t130 = t81 * qJD(2);
t129 = t81 * qJD(4);
t128 = t90 * qJD(2);
t127 = t102 * qJD(2);
t126 = t50 * t137;
t125 = t180 * t184;
t124 = t175 * t130;
t121 = t90 * t105;
t119 = pkin(4) * t194;
t118 = qJD(2) * t94 + qJD(3);
t6 = -t50 * t174 - t180 * t63;
t116 = -t6 * qJD(2) - t203;
t7 = -t174 * t180 + t50 * t63;
t115 = -t7 * qJD(2) - t204;
t68 = (0.1e1 / 0.2e1 - t97 / 0.2e1 - t96 / 0.2e1) * t102;
t82 = t90 * qJ(3);
t114 = t68 * qJD(1) - t82 * qJD(2);
t113 = qJD(4) * t180 + t23 * qJD(5);
t107 = t149 * t180 - t203;
t106 = -t50 * t149 - t204;
t69 = (0.1e1 + t90) * t102 / 0.2e1;
t46 = (-t81 / 0.2e1 + t109) * t105;
t45 = (-t175 / 0.2e1 + t108) * t105;
t1 = [0, 0, 0, 0, 0, 0, 0, t54 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t127, -t105 * qJD(2), -t99 * t127, t98 * t127, qJD(2) * t121, t132 + (-t102 * pkin(2) + qJ(3) * t121) * qJD(2) + t69 * qJD(3), 0, 0, 0, 0, 0, t46 * qJD(4) - t127 * t175, t45 * qJD(4) + t81 * t127, 0, 0, 0, 0, 0, t50 * t127 + t194 * t193, t127 * t180 + t194 * t190; 0, 0, 0, 0, 0, 0, 0, t69 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * qJD(2) - t66 * qJD(4), t45 * qJD(2) + t64 * qJD(4), 0, 0, 0, 0, 0, t199, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t200; 0, 0, 0, 0, 0, 0, 0, -t68 * qJD(3) - t132, 0, 0, 0, 0, 0, -t43 * qJD(4), -t44 * qJD(4), 0, 0, 0, 0, 0, -t194 * t192, -t194 * t191; 0, 0, 0, 0, 0, 0, t90 * qJD(3), t82 * qJD(3), t175 * t129, t37 * qJD(4), 0, 0, 0, t94 * t129, t94 * t77, t9 * t180, t194 * t181, 0, 0, 0, -t6 * qJD(4) + t63 * t133, -t7 * qJD(4) - t63 * t136; 0, 0, 0, 0, 0, 0, t128, -t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0; 0, 0, 0, 0, 0, 0, 0, 0, t124, t142, t77, -t129, 0, t117 * qJD(4) + t94 * t130 - t141, t120 * qJD(4) + t94 * t131 - t140, -t125, t195, t9, -t113, 0, t116 + t197, t115 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t195, t9, -t23 * qJD(4) - t133, 0, t48 * qJD(3) + t107 + t197, t106 + t196; 0, 0, 0, 0, 0, 0, 0, t68 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t128, t114, 0, 0, 0, 0, 0, t129, t77, 0, 0, 0, 0, 0, t113, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t131, 0, 0, 0, 0, 0, t137, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * qJD(2), t44 * qJD(2), 0, 0, 0, 0, 0, t201, t202; 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t142, 0, 0, 0, -t118 * t81 + t141, -t118 * t175 + t140, t125, -t195, 0, -t138, 0, -qJD(3) * t180 - t116, -t115 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t131, 0, 0, 0, 0, 0, -t137, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t170, -t103 * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, -t100 * t119, -t103 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t195, 0, t48 * qJD(4), 0, -t23 * qJD(3) - t107, -t106 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t100 * t169, t103 * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
