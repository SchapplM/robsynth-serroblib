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
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 11:14:07
% EndTime: 2021-01-15 11:14:09
% DurationCPUTime: 0.71s
% Computational Cost: add. (1378->98), mult. (2565->138), div. (0->0), fcn. (2687->6), ass. (0->79)
t130 = cos(pkin(8));
t88 = sin(pkin(8));
t90 = sin(qJ(3));
t91 = cos(qJ(3));
t76 = t130 * t90 + t88 * t91;
t73 = t76 ^ 2;
t74 = -t130 * t91 + t88 * t90;
t42 = t74 ^ 2 + t73;
t153 = t42 * qJD(1);
t152 = t42 * qJD(4);
t83 = sin(pkin(7)) * pkin(1) + pkin(6);
t151 = qJ(4) + t83;
t71 = t151 * t90;
t137 = t88 * t71;
t72 = t151 * t91;
t56 = t130 * t72;
t103 = t56 - t137;
t38 = t130 * t71 + t88 * t72;
t99 = -t103 * t74 + t38 * t76;
t148 = qJD(4) * t99;
t144 = t99 * qJD(1);
t101 = t56 / 0.2e1;
t142 = t76 * pkin(4);
t141 = t90 * pkin(3);
t140 = t74 * t88;
t81 = t88 * pkin(3) + qJ(5);
t139 = t81 * t74;
t84 = -t130 * pkin(3) - pkin(4);
t138 = t84 * t76;
t136 = qJD(3) * pkin(3);
t132 = t74 * qJ(5);
t85 = -cos(pkin(7)) * pkin(1) - pkin(2);
t78 = -t91 * pkin(3) + t85;
t34 = t74 * pkin(4) - t76 * qJ(5) + t78;
t41 = t132 + t141 + t142;
t9 = t34 * t76 + t41 * t74;
t131 = t9 * qJD(1);
t129 = qJD(1) * t91;
t10 = t34 * t74 - t41 * t76;
t128 = t10 * qJD(1);
t86 = t141 / 0.2e1;
t17 = t86 + (pkin(4) / 0.2e1 - t84 / 0.2e1) * t76 + (qJ(5) / 0.2e1 + t81 / 0.2e1) * t74;
t124 = t17 * qJD(1);
t102 = t130 * t76;
t94 = -t140 / 0.2e1 - t102 / 0.2e1;
t31 = (-t90 / 0.2e1 + t94) * pkin(3);
t123 = t31 * qJD(1);
t32 = t74 * t141 + t78 * t76;
t122 = t32 * qJD(1);
t33 = t76 * t141 - t78 * t74;
t121 = t33 * qJD(1);
t120 = t38 * qJD(3);
t117 = t73 * qJD(1);
t116 = t74 * qJD(1);
t115 = t74 * qJD(3);
t114 = t74 * qJD(4);
t113 = t76 * qJD(1);
t112 = t76 * qJD(3);
t111 = t76 * qJD(5);
t80 = -t90 ^ 2 + t91 ^ 2;
t110 = t80 * qJD(1);
t109 = t90 * qJD(3);
t108 = t91 * qJD(3);
t107 = t74 * t113;
t106 = t85 * t90 * qJD(1);
t105 = t85 * t129;
t104 = t90 * t129;
t3 = t34 * t41;
t98 = t3 * qJD(1);
t8 = t78 * t141;
t97 = t8 * qJD(1);
t36 = t101 - t56 / 0.2e1;
t95 = t36 * qJD(1) + t81 * qJD(3);
t68 = t76 * qJD(4);
t35 = t103 * qJD(3);
t30 = t94 * pkin(3) + t86;
t19 = 0.2e1 * t101 - t137;
t18 = -t139 / 0.2e1 + t138 / 0.2e1 + t86 + t132 / 0.2e1 + t142 / 0.2e1;
t1 = [0, 0, 0, 0, t90 * t108, t80 * qJD(3), 0, 0, 0, t85 * t109, t85 * t108, t32 * qJD(3), t33 * qJD(3), t152, t8 * qJD(3) + t148, t9 * qJD(3) - t111 * t74, t152, t10 * qJD(3) + t73 * qJD(5), t3 * qJD(3) - t111 * t34 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t104, t110, t108, -t109, 0, -t83 * t108 + t106, t109 * t83 + t105, -t35 + t122, t120 + t121, (t130 * t74 - t76 * t88) * t136, (-t103 * t130 - t38 * t88) * t136 + t30 * qJD(4) + t97, -t35 + t131, (-t84 * t74 - t81 * t76) * qJD(3) - qJD(5) * t74, -t120 + t128, (t103 * t84 - t38 * t81) * qJD(3) + t18 * qJD(4) + t19 * qJD(5) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t30 * qJD(3) + t144, 0, t153, 0, t18 * qJD(3) + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t115, t117, t19 * qJD(3) - t113 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t108, -t112, t115, 0, (-t102 - t140) * t136, -t112, 0, -t115, (t138 - t139) * qJD(3) + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112; 0, 0, 0, 0, -t104, -t110, 0, 0, 0, -t106, -t105, -t68 - t122, t114 - t121, 0, t31 * qJD(4) - t97, -t68 - t131, 0, -t114 - t128, -t17 * qJD(4) + t36 * qJD(5) - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t81 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t116, 0, t123, -t113, 0, -t116, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t115, -t153, -t31 * qJD(3) - t144, t112, -t153, t115, t17 * qJD(3) - t111 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, -t116, 0, -t123, t113, 0, t116, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, -t117, -t36 * qJD(3) + (qJD(1) * t34 + qJD(4)) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
