% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:34
% EndTime: 2019-12-05 19:00:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (215->32), mult. (309->76), div. (0->0), fcn. (197->8), ass. (0->34)
t113 = qJD(1) + qJD(2);
t111 = t113 ^ 2;
t132 = t111 / 0.2e1;
t131 = cos(qJ(4));
t130 = pkin(1) * qJD(1);
t116 = sin(qJ(3));
t124 = sin(qJ(2)) * t130;
t107 = t113 * pkin(7) + t124;
t122 = pkin(8) * t113 + t107;
t101 = qJD(3) * pkin(3) - t122 * t116;
t119 = cos(qJ(3));
t102 = t122 * t119;
t115 = sin(qJ(4));
t129 = t115 * t101 + t131 * t102;
t128 = t113 * t116;
t127 = t113 * t119;
t126 = qJD(3) * t116;
t125 = qJD(3) * t119;
t112 = qJD(3) + qJD(4);
t123 = cos(qJ(2)) * t130;
t121 = t131 * t101 - t115 * t102;
t106 = -t123 + (-pkin(3) * t119 - pkin(2)) * t113;
t118 = cos(qJ(5));
t114 = sin(qJ(5));
t110 = qJD(5) + t112;
t108 = -t113 * pkin(2) - t123;
t105 = (t115 * t119 + t131 * t116) * t113;
t104 = t115 * t128 - t131 * t127;
t97 = t104 * pkin(4) + t106;
t96 = -t114 * t104 + t118 * t105;
t95 = t118 * t104 + t114 * t105;
t94 = -t104 * pkin(9) + t129;
t93 = t112 * pkin(4) - t105 * pkin(9) + t121;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t132, t113 * t123, -t113 * t124, t116 ^ 2 * t132, t116 * t111 * t119, t113 * t126, t113 * t125, qJD(3) ^ 2 / 0.2e1, -t107 * t126 - t108 * t127, -t107 * t125 + t108 * t128, t105 ^ 2 / 0.2e1, -t105 * t104, t105 * t112, -t104 * t112, t112 ^ 2 / 0.2e1, t106 * t104 + t121 * t112, t106 * t105 - t129 * t112, t96 ^ 2 / 0.2e1, -t96 * t95, t96 * t110, -t95 * t110, t110 ^ 2 / 0.2e1, t97 * t95 + (-t114 * t94 + t118 * t93) * t110, t97 * t96 - (t114 * t93 + t118 * t94) * t110;];
T_reg = t1;
