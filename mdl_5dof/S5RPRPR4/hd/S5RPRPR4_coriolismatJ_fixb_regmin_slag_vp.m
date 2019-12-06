% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:54:00
% EndTime: 2019-12-05 17:54:04
% DurationCPUTime: 0.85s
% Computational Cost: add. (1446->100), mult. (2743->142), div. (0->0), fcn. (3029->8), ass. (0->92)
t153 = cos(qJ(5));
t88 = sin(pkin(8)) * pkin(1) + pkin(6);
t140 = qJ(4) + t88;
t96 = sin(qJ(3));
t81 = t140 * t96;
t97 = cos(qJ(3));
t82 = t140 * t97;
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t158 = -t93 * t81 - t92 * t82;
t84 = t92 * t97 + t93 * t96;
t162 = -t84 * pkin(7) + t158;
t165 = t153 * t162;
t106 = -t165 / 0.2e1;
t41 = t92 * t81 - t93 * t82;
t83 = -t92 * t96 + t93 * t97;
t36 = t83 * pkin(7) - t41;
t167 = t153 * t36;
t107 = -t167 / 0.2e1;
t95 = sin(qJ(5));
t168 = t95 * t36;
t170 = -t165 + t168;
t166 = t95 * t162;
t169 = -t167 - t166;
t145 = t95 * t83;
t78 = t153 * t84;
t157 = t78 + t145;
t108 = -t153 * t83 + t95 * t84;
t160 = t108 * qJD(1);
t164 = t157 * t160;
t159 = t108 ^ 2 - t157 ^ 2;
t163 = t159 * qJD(1);
t161 = qJD(4) * t108;
t128 = t108 * qJD(5);
t101 = qJD(3) * t108 + t128;
t105 = t78 / 0.2e1;
t156 = pkin(3) * t92;
t154 = t96 * pkin(3);
t152 = t83 * t92;
t151 = t84 * t93;
t144 = qJD(3) * pkin(3);
t90 = -cos(pkin(8)) * pkin(1) - pkin(2);
t102 = -t97 * pkin(3) + t90;
t60 = -t83 * pkin(4) + t102;
t139 = qJD(1) * t60;
t138 = qJD(1) * t97;
t61 = t84 * pkin(4) + t154;
t10 = t108 * t61 + t157 * t60;
t137 = t10 * qJD(1);
t11 = -t108 * t60 + t157 * t61;
t136 = t11 * qJD(1);
t13 = -t158 * t84 - t41 * t83;
t134 = t13 * qJD(1);
t38 = 0.2e1 * t105 + t145;
t132 = t38 * qJD(1);
t99 = t152 / 0.2e1 - t151 / 0.2e1;
t40 = (-t96 / 0.2e1 + t99) * pkin(3);
t131 = t40 * qJD(1);
t49 = t105 - t78 / 0.2e1;
t130 = t49 * qJD(1);
t47 = t49 * qJD(5);
t129 = t157 * qJD(1);
t125 = t157 * qJD(5);
t55 = t83 ^ 2 + t84 ^ 2;
t124 = t55 * qJD(1);
t86 = -t96 ^ 2 + t97 ^ 2;
t123 = t86 * qJD(1);
t122 = t96 * qJD(3);
t121 = t97 * qJD(3);
t118 = t108 * t139;
t117 = t157 * t139;
t116 = t90 * t96 * qJD(1);
t115 = t90 * t138;
t114 = t96 * t138;
t8 = t102 * t154;
t104 = t8 * qJD(1);
t3 = t106 + t165 / 0.2e1;
t89 = t93 * pkin(3) + pkin(4);
t75 = -t153 * t89 + t95 * t156;
t103 = -t3 * qJD(1) - t75 * qJD(3);
t100 = qJD(3) * t157 + t38 * qJD(5);
t2 = t107 + t167 / 0.2e1;
t76 = t153 * t156 + t95 * t89;
t98 = t2 * qJD(1) + t49 * qJD(2) + t76 * qJD(3);
t63 = t76 * qJD(5);
t62 = t75 * qJD(5);
t46 = t49 * qJD(3);
t39 = t154 / 0.2e1 + t99 * pkin(3);
t14 = -t38 * qJD(3) - t125;
t5 = 0.2e1 * t107 - t166;
t4 = t168 + 0.2e1 * t106;
t1 = [0, 0, 0, 0, t96 * t121, t86 * qJD(3), 0, 0, 0, t90 * t122, t90 * t121, t55 * qJD(4), t8 * qJD(3) + t13 * qJD(4), -t101 * t157, (qJD(3) + qJD(5)) * t159, 0, 0, 0, t10 * qJD(3) + t60 * t125, t11 * qJD(3) - t60 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t114, t123, t121, -t122, 0, -t88 * t121 + t116, t88 * t122 + t115, (-t83 * t93 - t84 * t92) * t144, t39 * qJD(4) + (t158 * t92 + t41 * t93) * t144 + t104, -t164, t163, -t101, -t100, 0, t169 * qJD(3) + t5 * qJD(5) + t137, qJD(3) * t170 + t4 * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t39 * qJD(3) + t134, 0, 0, 0, 0, 0, t47, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t163, -t101, t14, 0, t5 * qJD(3) + t49 * qJD(4) + t169 * qJD(5) + t117, t4 * qJD(3) + qJD(5) * t170 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121, 0, (-t151 + t152) * t144, 0, 0, 0, 0, 0, -t100, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t101; 0, 0, 0, 0, -t114, -t123, 0, 0, 0, -t116, -t115, 0, t40 * qJD(4) - t104, t164, -t163, 0, -t47, 0, -qJD(4) * t157 - t2 * qJD(5) - t137, t3 * qJD(5) - t136 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, 0, 0, 0, -t129, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, 0, -t63 - t98, -t103 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t40 * qJD(3) - t134, 0, 0, 0, 0, 0, t100, -t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, 0, 0, 0, 0, t129, -t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t163, 0, t46, 0, t2 * qJD(3) - t38 * qJD(4) - t117, -t3 * qJD(3) + t118 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, t98, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
