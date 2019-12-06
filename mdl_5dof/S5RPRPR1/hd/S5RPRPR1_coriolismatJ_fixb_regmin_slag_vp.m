% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:47
% EndTime: 2019-12-05 17:47:51
% DurationCPUTime: 1.10s
% Computational Cost: add. (1476->100), mult. (2577->137), div. (0->0), fcn. (2899->6), ass. (0->90)
t185 = qJD(3) + qJD(5);
t101 = cos(qJ(5));
t148 = sin(pkin(8));
t100 = sin(qJ(3));
t103 = -pkin(1) - pkin(6);
t170 = -qJ(4) + t103;
t88 = t170 * t100;
t102 = cos(qJ(3));
t89 = t170 * t102;
t98 = cos(pkin(8));
t157 = -t148 * t88 + t98 * t89;
t85 = -t148 * t100 + t98 * t102;
t169 = -t85 * pkin(7) + t157;
t53 = -t148 * t89 - t98 * t88;
t86 = -t98 * t100 - t148 * t102;
t38 = t86 * pkin(7) - t53;
t99 = sin(qJ(5));
t191 = t185 * (-t101 * t38 - t99 * t169);
t190 = t185 * (-t101 * t169 + t99 * t38);
t23 = t101 * t86 - t99 * t85;
t136 = t23 * qJD(5);
t15 = t23 * qJD(3) + t136;
t158 = t99 * t86;
t79 = t101 * t85;
t172 = t79 + t158;
t133 = t172 * qJD(5);
t175 = qJD(3) * t172;
t182 = -t133 - t175;
t173 = -t172 ^ 2 + t23 ^ 2;
t181 = t173 * qJD(1);
t178 = qJD(4) * t23;
t176 = t23 * qJD(1);
t174 = t172 * qJD(1);
t156 = qJD(3) * pkin(3);
t171 = (t148 * t85 + t86 * t98) * t156;
t168 = t85 ^ 2;
t167 = t86 ^ 2;
t116 = t79 / 0.2e1;
t165 = t102 * pkin(3);
t95 = t100 * pkin(3) + qJ(2);
t63 = -t86 * pkin(4) + t95;
t64 = t85 * pkin(4) + t165;
t12 = t172 * t63 - t23 * t64;
t147 = t12 * qJD(1);
t13 = t172 * t64 + t23 * t63;
t146 = t13 * qJD(1);
t16 = -t157 * t85 - t53 * t86;
t144 = t16 * qJD(1);
t110 = -t168 / 0.2e1 - t167 / 0.2e1;
t18 = -0.1e1 / 0.2e1 + t110;
t143 = t18 * qJD(1);
t22 = 0.2e1 * t116 + t158;
t141 = t22 * qJD(1);
t105 = t148 * t86 / 0.2e1 - t98 * t85 / 0.2e1;
t41 = (-t102 / 0.2e1 + t105) * pkin(3);
t140 = t41 * qJD(1);
t45 = t116 - t79 / 0.2e1;
t139 = t45 * qJD(1);
t138 = t45 * qJD(5);
t56 = t167 + t168;
t130 = t56 * qJD(1);
t90 = t100 ^ 2 - t102 ^ 2;
t129 = t90 * qJD(1);
t128 = t95 * qJD(1);
t127 = qJD(1) * qJ(2);
t126 = t100 * qJD(1);
t125 = t100 * qJD(3);
t124 = t102 * qJD(1);
t123 = t102 * qJD(3);
t122 = t23 * t174;
t121 = t172 * t176;
t120 = t63 * t176;
t119 = t63 * t174;
t114 = pkin(3) * t148;
t113 = qJ(2) * t126;
t112 = qJ(2) * t124;
t111 = t100 * t124;
t9 = t95 * t165;
t109 = t9 * qJD(1);
t94 = t98 * pkin(3) + pkin(4);
t76 = -t101 * t94 + t99 * t114;
t108 = t76 * qJD(3);
t107 = t22 * qJD(5) + t175;
t77 = t101 * t114 + t99 * t94;
t104 = t77 * qJD(3);
t66 = t77 * qJD(5);
t65 = t76 * qJD(5);
t40 = t165 / 0.2e1 + t105 * pkin(3);
t17 = 0.1e1 / 0.2e1 + t110;
t1 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t100 * t123, t90 * qJD(3), 0, 0, 0, qJ(2) * t123 + qJD(2) * t100, -qJ(2) * t125 + qJD(2) * t102, t56 * qJD(4), t95 * qJD(2) + t9 * qJD(3) + t16 * qJD(4), t15 * t172, t185 * t173, 0, 0, 0, -qJD(2) * t23 + t12 * qJD(3) + t63 * t133, qJD(2) * t172 + t13 * qJD(3) + t136 * t63; 0, 0, 0, 0, qJD(1), t127, 0, 0, 0, 0, 0, t126, t124, 0, t17 * qJD(4) + t128, 0, 0, 0, 0, 0, -t176, t174; 0, 0, 0, 0, 0, 0, -t111, t129, -t125, -t123, 0, -t103 * t125 + t112, -t103 * t123 - t113, -t171, (t148 * t157 + t53 * t98) * t156 + t40 * qJD(4) + t109, t121, t181, t15, -t107, 0, t147 + t191, t146 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t17 * qJD(2) + t40 * qJD(3) + t144, 0, 0, 0, 0, 0, t138, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t181, t15, -t22 * qJD(3) - t133, 0, t45 * qJD(4) + t119 + t191, t120 + t190; 0, 0, 0, 0, -qJD(1), -t127, 0, 0, 0, 0, 0, -t126, -t124, 0, t18 * qJD(4) - t128, 0, 0, 0, 0, 0, t176, -t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t123, 0, t171, 0, 0, 0, 0, 0, t15, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185 * t23, t182; 0, 0, 0, 0, 0, 0, t111, -t129, 0, 0, 0, -t112, t113, 0, t41 * qJD(4) - t109, -t121, -t181, 0, -t138, 0, -qJD(4) * t172 - t147, -t146 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, 0, 0, 0, 0, 0, -t174, -t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, -t104 - t66, t108 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t18 * qJD(2) - t41 * qJD(3) - t144, 0, 0, 0, 0, 0, t107, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, 0, 0, 0, 0, 0, t174, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t181, 0, t45 * qJD(3), 0, -t22 * qJD(4) - t119, -t120 - t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t104, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
