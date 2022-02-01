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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:08:21
% EndTime: 2022-01-20 12:08:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (215->32), mult. (309->76), div. (0->0), fcn. (197->8), ass. (0->34)
t114 = qJD(1) + qJD(2);
t112 = t114 ^ 2;
t133 = t112 / 0.2e1;
t132 = cos(qJ(4));
t131 = pkin(1) * qJD(1);
t117 = sin(qJ(3));
t125 = sin(qJ(2)) * t131;
t108 = t114 * pkin(7) + t125;
t123 = pkin(8) * t114 + t108;
t102 = qJD(3) * pkin(3) - t123 * t117;
t120 = cos(qJ(3));
t103 = t123 * t120;
t116 = sin(qJ(4));
t130 = t116 * t102 + t132 * t103;
t129 = t114 * t117;
t128 = t114 * t120;
t127 = qJD(3) * t117;
t126 = qJD(3) * t120;
t113 = qJD(3) + qJD(4);
t124 = cos(qJ(2)) * t131;
t122 = t132 * t102 - t116 * t103;
t107 = -t124 + (-pkin(3) * t120 - pkin(2)) * t114;
t119 = cos(qJ(5));
t115 = sin(qJ(5));
t111 = qJD(5) + t113;
t109 = -t114 * pkin(2) - t124;
t106 = (t116 * t120 + t132 * t117) * t114;
t105 = t116 * t129 - t132 * t128;
t98 = t105 * pkin(4) + t107;
t97 = -t115 * t105 + t119 * t106;
t96 = t119 * t105 + t115 * t106;
t95 = -t105 * pkin(9) + t130;
t94 = t113 * pkin(4) - t106 * pkin(9) + t122;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t133, t114 * t124, -t114 * t125, t117 ^ 2 * t133, t117 * t112 * t120, t114 * t127, t114 * t126, qJD(3) ^ 2 / 0.2e1, -t108 * t127 - t109 * t128, -t108 * t126 + t109 * t129, t106 ^ 2 / 0.2e1, -t106 * t105, t106 * t113, -t105 * t113, t113 ^ 2 / 0.2e1, t107 * t105 + t122 * t113, t107 * t106 - t130 * t113, t97 ^ 2 / 0.2e1, -t97 * t96, t97 * t111, -t96 * t111, t111 ^ 2 / 0.2e1, t98 * t96 + (-t115 * t95 + t119 * t94) * t111, t98 * t97 - (t115 * t94 + t119 * t95) * t111;];
T_reg = t1;
