% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:58
% EndTime: 2019-12-05 17:12:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->31), mult. (287->76), div. (0->0), fcn. (197->8), ass. (0->34)
t125 = qJD(2) ^ 2;
t136 = t125 / 0.2e1;
t135 = cos(qJ(4));
t120 = sin(qJ(3));
t121 = sin(qJ(2));
t113 = qJD(2) * pkin(6) + t121 * qJD(1);
t127 = pkin(7) * qJD(2) + t113;
t107 = qJD(3) * pkin(3) - t127 * t120;
t123 = cos(qJ(3));
t109 = t127 * t123;
t119 = sin(qJ(4));
t134 = t119 * t107 + t135 * t109;
t133 = qJD(2) * t120;
t132 = qJD(2) * t123;
t131 = qJD(3) * t113;
t124 = cos(qJ(2));
t130 = t124 * qJD(1);
t129 = qJD(1) * qJD(2);
t128 = qJD(2) * qJD(3);
t117 = qJD(3) + qJD(4);
t126 = t135 * t107 - t119 * t109;
t112 = -t130 + (-pkin(3) * t123 - pkin(2)) * qJD(2);
t122 = cos(qJ(5));
t118 = sin(qJ(5));
t116 = qJD(5) + t117;
t114 = -qJD(2) * pkin(2) - t130;
t111 = (t119 * t123 + t135 * t120) * qJD(2);
t110 = t119 * t133 - t135 * t132;
t103 = t110 * pkin(4) + t112;
t102 = -t118 * t110 + t122 * t111;
t101 = t122 * t110 + t118 * t111;
t100 = -t110 * pkin(8) + t134;
t99 = t117 * pkin(4) - t111 * pkin(8) + t126;
t1 = [qJD(1) ^ 2 / 0.2e1, t136, t124 * t129, -t121 * t129, t120 ^ 2 * t136, t120 * t125 * t123, t120 * t128, t123 * t128, qJD(3) ^ 2 / 0.2e1, -t114 * t132 - t120 * t131, t114 * t133 - t123 * t131, t111 ^ 2 / 0.2e1, -t111 * t110, t111 * t117, -t110 * t117, t117 ^ 2 / 0.2e1, t112 * t110 + t126 * t117, t112 * t111 - t134 * t117, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t116, -t101 * t116, t116 ^ 2 / 0.2e1, (-t118 * t100 + t122 * t99) * t116 + t103 * t101, -(t122 * t100 + t118 * t99) * t116 + t103 * t102;];
T_reg = t1;
