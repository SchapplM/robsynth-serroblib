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
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 15:22:43
% EndTime: 2021-01-15 15:22:45
% DurationCPUTime: 0.65s
% Computational Cost: add. (1103->96), mult. (2290->136), div. (0->0), fcn. (2412->4), ass. (0->77)
t127 = cos(pkin(8));
t86 = sin(pkin(8));
t87 = sin(qJ(3));
t88 = cos(qJ(3));
t72 = t127 * t87 + t86 * t88;
t69 = t72 ^ 2;
t70 = -t127 * t88 + t86 * t87;
t148 = t70 ^ 2 + t69;
t151 = qJD(4) * t148;
t150 = t148 * qJD(2);
t149 = -qJ(4) - pkin(6);
t76 = t149 * t87;
t134 = t86 * t76;
t77 = t149 * t88;
t75 = t127 * t77;
t100 = -t75 + t134;
t46 = -t127 * t76 - t86 * t77;
t96 = -t100 * t70 + t46 * t72;
t145 = qJD(4) * t96;
t141 = t96 * qJD(2);
t98 = -t75 / 0.2e1;
t139 = t72 * pkin(4);
t138 = t87 * pkin(3);
t137 = t70 * t86;
t79 = t86 * pkin(3) + qJ(5);
t136 = t79 * t70;
t81 = -t127 * pkin(3) - pkin(4);
t135 = t81 * t72;
t133 = qJD(3) * pkin(3);
t129 = t70 * qJ(5);
t82 = -t88 * pkin(3) - pkin(2);
t34 = t70 * pkin(4) - t72 * qJ(5) + t82;
t35 = t129 + t138 + t139;
t9 = t34 * t72 + t35 * t70;
t128 = t9 * qJD(2);
t126 = qJD(2) * t88;
t10 = t34 * t70 - t35 * t72;
t125 = t10 * qJD(2);
t83 = t138 / 0.2e1;
t17 = t83 + (pkin(4) / 0.2e1 - t81 / 0.2e1) * t72 + (qJ(5) / 0.2e1 + t79 / 0.2e1) * t70;
t121 = t17 * qJD(2);
t99 = t127 * t72;
t91 = -t137 / 0.2e1 - t99 / 0.2e1;
t24 = (-t87 / 0.2e1 + t91) * pkin(3);
t120 = t24 * qJD(2);
t32 = t70 * t138 + t82 * t72;
t119 = t32 * qJD(2);
t33 = t72 * t138 - t82 * t70;
t118 = t33 * qJD(2);
t115 = t46 * qJD(3);
t114 = t69 * qJD(2);
t113 = t70 * qJD(2);
t112 = t70 * qJD(3);
t111 = t70 * qJD(4);
t110 = t72 * qJD(2);
t109 = t72 * qJD(3);
t108 = t72 * qJD(5);
t78 = -t87 ^ 2 + t88 ^ 2;
t107 = t78 * qJD(2);
t106 = t87 * qJD(3);
t105 = t88 * qJD(3);
t104 = pkin(2) * t87 * qJD(2);
t103 = pkin(2) * t126;
t102 = t70 * t110;
t101 = t87 * t126;
t3 = t34 * t35;
t95 = t3 * qJD(2);
t8 = t82 * t138;
t94 = t8 * qJD(2);
t44 = t98 + t75 / 0.2e1;
t92 = t44 * qJD(2) + t79 * qJD(3);
t64 = t72 * qJD(4);
t36 = t100 * qJD(3);
t31 = 0.2e1 * t98 + t134;
t23 = t91 * pkin(3) + t83;
t18 = t135 / 0.2e1 - t136 / 0.2e1 + t83 + t129 / 0.2e1 + t139 / 0.2e1;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, -t109, t112, 0, (-t99 - t137) * t133, -t109, 0, -t112, (t135 - t136) * qJD(3) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t87 * t105, t78 * qJD(3), 0, 0, 0, -pkin(2) * t106, -pkin(2) * t105, t32 * qJD(3), t33 * qJD(3), t151, t8 * qJD(3) + t145, t9 * qJD(3) - t70 * t108, t151, t10 * qJD(3) + t69 * qJD(5), t3 * qJD(3) - t34 * t108 + t145; 0, 0, 0, 0, t101, t107, t105, -t106, 0, -pkin(6) * t105 - t104, pkin(6) * t106 - t103, -t36 + t119, t115 + t118, (t127 * t70 - t72 * t86) * t133, (-t100 * t127 - t46 * t86) * t133 + t23 * qJD(4) + t94, -t36 + t128, (-t81 * t70 - t79 * t72) * qJD(3) - qJD(5) * t70, -t115 + t125, (t100 * t81 - t46 * t79) * qJD(3) + t18 * qJD(4) + t31 * qJD(5) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t23 * qJD(3) + t141, 0, t150, 0, t18 * qJD(3) + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t112, t114, t31 * qJD(3) - t34 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t101, -t107, 0, 0, 0, t104, t103, -t64 - t119, t111 - t118, 0, t24 * qJD(4) - t94, -t64 - t128, 0, -t111 - t125, -t17 * qJD(4) + t44 * qJD(5) - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t79 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t113, 0, t120, -t110, 0, -t113, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t112, -t150, -t24 * qJD(3) - t141, t109, -t150, t112, t17 * qJD(3) - t108 - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t113, 0, -t120, t110, 0, t113, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, -t114, -t44 * qJD(3) + (qJD(2) * t34 + qJD(4)) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
