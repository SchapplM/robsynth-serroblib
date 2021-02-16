% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:05
% EndTime: 2021-01-15 11:24:08
% DurationCPUTime: 0.85s
% Computational Cost: add. (1434->124), mult. (2413->150), div. (0->0), fcn. (2592->4), ass. (0->94)
t142 = sin(pkin(7));
t95 = -pkin(1) - pkin(6);
t161 = -qJ(4) + t95;
t94 = cos(qJ(3));
t84 = t161 * t94;
t107 = t142 * t84;
t93 = sin(qJ(3));
t83 = t161 * t93;
t92 = cos(pkin(7));
t75 = t92 * t83;
t160 = t107 + t75;
t46 = t142 * t83 - t92 * t84;
t79 = -t142 * t93 + t92 * t94;
t80 = t142 * t94 + t92 * t93;
t104 = -t160 * t80 + t46 * t79;
t154 = t104 * qJD(1);
t77 = t80 ^ 2;
t78 = t79 ^ 2;
t106 = -t78 / 0.2e1 - t77 / 0.2e1;
t28 = -0.1e1 / 0.2e1 + t106;
t167 = -t28 * qJD(2) - t154;
t27 = 0.1e1 / 0.2e1 + t106;
t166 = t27 * qJD(2) + t154;
t163 = t77 + t78;
t165 = qJD(4) * t163;
t164 = t163 * qJD(1);
t148 = qJD(3) * pkin(3);
t162 = (t142 * t79 - t80 * t92) * t148;
t153 = t104 * qJD(4);
t121 = t80 * qJD(5);
t86 = t142 * pkin(3) + qJ(5);
t89 = -t92 * pkin(3) - pkin(4);
t152 = -(t79 * t86 + t80 * t89) * qJD(3) - t121;
t110 = t75 / 0.2e1;
t151 = -t80 / 0.2e1;
t150 = t79 * pkin(4);
t149 = t94 * pkin(3);
t143 = t80 * qJ(5);
t88 = t93 * pkin(3) + qJ(2);
t42 = t80 * pkin(4) - t79 * qJ(5) + t88;
t43 = t143 + t149 + t150;
t15 = t42 * t79 + t43 * t80;
t141 = t15 * qJD(1);
t16 = t42 * t80 - t43 * t79;
t140 = t16 * qJD(1);
t91 = t149 / 0.2e1;
t19 = t91 + (qJ(5) / 0.2e1 + t86 / 0.2e1) * t80 + (pkin(4) / 0.2e1 - t89 / 0.2e1) * t79;
t137 = t19 * qJD(1);
t136 = t28 * qJD(1);
t98 = t142 * t151 - t92 * t79 / 0.2e1;
t39 = (-t94 / 0.2e1 + t98) * pkin(3);
t135 = t39 * qJD(1);
t40 = t79 * t149 - t88 * t80;
t134 = t40 * qJD(1);
t41 = t80 * t149 + t88 * t79;
t133 = t41 * qJD(1);
t132 = t42 * qJD(1);
t131 = t46 * qJD(3);
t128 = t78 * qJD(1);
t127 = t79 * qJD(1);
t126 = t79 * qJD(2);
t66 = t79 * qJD(3);
t125 = t79 * qJD(5);
t124 = t80 * qJD(1);
t123 = t80 * qJD(3);
t122 = t80 * qJD(4);
t85 = t93 ^ 2 - t94 ^ 2;
t120 = t85 * qJD(1);
t119 = t88 * qJD(1);
t118 = t93 * qJD(1);
t117 = t93 * qJD(3);
t116 = t94 * qJD(1);
t115 = t94 * qJD(3);
t114 = qJ(2) * qJD(3);
t113 = qJD(1) * qJ(2);
t112 = t79 * t124;
t111 = t93 * t116;
t109 = t93 * t113;
t108 = t94 * t113;
t1 = t42 * t43;
t102 = t1 * qJD(1);
t8 = t88 * t149;
t101 = t8 * qJD(1);
t45 = t110 - t75 / 0.2e1;
t100 = t45 * qJD(1) + t86 * qJD(3);
t73 = t79 * qJD(4);
t69 = t80 * qJD(2);
t44 = t160 * qJD(3);
t38 = pkin(3) * t98 + t91;
t31 = 0.2e1 * t110 + t107;
t23 = t28 * qJD(4);
t21 = t27 * qJD(4);
t20 = t89 * t79 / 0.2e1 + t86 * t151 + t91 + t143 / 0.2e1 + t150 / 0.2e1;
t2 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t93 * t115, t85 * qJD(3), 0, 0, 0, qJD(2) * t93 + t94 * t114, qJD(2) * t94 - t93 * t114, t41 * qJD(3) + t69, t40 * qJD(3) + t126, t165, t88 * qJD(2) + t8 * qJD(3) + t153, t15 * qJD(3) - t79 * t121 + t69, t165, t16 * qJD(3) + t78 * qJD(5) - t126, t1 * qJD(3) + t153 + (qJD(2) - t125) * t42; 0, 0, 0, 0, qJD(1), t113, 0, 0, 0, 0, 0, t118, t116, t124, t127, 0, t21 + t119, t124, 0, -t127, t21 + t132; 0, 0, 0, 0, 0, 0, -t111, t120, -t117, -t115, 0, -t95 * t117 + t108, -t95 * t115 - t109, -t44 + t133, t131 + t134, -t162, (-t142 * t46 - t160 * t92) * t148 + t38 * qJD(4) + t101, -t44 + t141, t152, -t131 + t140, (t160 * t89 - t46 * t86) * qJD(3) + t20 * qJD(4) + t31 * qJD(5) + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t38 * qJD(3) + t166, 0, t164, 0, t20 * qJD(3) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t123, t128, t31 * qJD(3) - t42 * t127; 0, 0, 0, 0, -qJD(1), -t113, 0, 0, 0, 0, 0, -t118, -t116, -t124, -t127, 0, t23 - t119, -t124, 0, t127, t23 - t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t115, -t123, -t66, 0, t162, -t123, 0, t66, -t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123; 0, 0, 0, 0, 0, 0, t111, -t120, 0, 0, 0, -t108, t109, -t73 - t133, t122 - t134, 0, t39 * qJD(4) - t101, -t73 - t141, 0, -t122 - t140, -t19 * qJD(4) + t45 * qJD(5) - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t86 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t124, 0, t135, -t127, 0, -t124, -t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t123, -t164, -t39 * qJD(3) + t167, t66, -t164, t123, t19 * qJD(3) - t125 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, 0, 0, 0, -t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t124, 0, -t135, t127, 0, t124, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, -t128, -t45 * qJD(3) + (qJD(4) + t132) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
