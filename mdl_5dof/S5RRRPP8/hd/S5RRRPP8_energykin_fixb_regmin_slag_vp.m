% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:23
% EndTime: 2019-12-31 21:09:23
% DurationCPUTime: 0.12s
% Computational Cost: add. (183->36), mult. (394->76), div. (0->0), fcn. (219->4), ass. (0->30)
t125 = qJD(1) ^ 2;
t137 = t125 / 0.2e1;
t136 = pkin(3) + qJ(5);
t124 = cos(qJ(2));
t135 = t124 * t125;
t122 = sin(qJ(2));
t110 = (-pkin(2) * t124 - pkin(7) * t122 - pkin(1)) * qJD(1);
t132 = t124 * qJD(1);
t116 = pkin(6) * t132 + qJD(2) * pkin(7);
t121 = sin(qJ(3));
t123 = cos(qJ(3));
t134 = t121 * t110 + t123 * t116;
t133 = qJD(1) * t122;
t131 = qJD(1) * qJD(2);
t130 = t122 * t131;
t129 = t124 * t131;
t115 = -qJD(2) * pkin(2) + pkin(6) * t133;
t128 = t123 * t110 - t121 * t116;
t117 = -qJD(3) + t132;
t107 = t117 * qJ(4) - t134;
t127 = qJD(4) - t128;
t112 = t121 * qJD(2) + t123 * t133;
t126 = -t112 * qJ(4) + t115;
t111 = -t123 * qJD(2) + t121 * t133;
t108 = t111 * pkin(3) + t126;
t106 = t117 * pkin(3) + t127;
t105 = t136 * t111 + t126;
t104 = -t111 * pkin(4) + qJD(5) - t107;
t103 = t112 * pkin(4) + t136 * t117 + t127;
t1 = [t137, 0, 0, t122 ^ 2 * t137, t122 * t135, t130, t129, qJD(2) ^ 2 / 0.2e1, pkin(1) * t135 - pkin(6) * t130, -t125 * pkin(1) * t122 - pkin(6) * t129, t112 ^ 2 / 0.2e1, -t112 * t111, -t112 * t117, t111 * t117, t117 ^ 2 / 0.2e1, t115 * t111 - t128 * t117, t115 * t112 + t134 * t117, t106 * t112 + t107 * t111, -t106 * t117 - t108 * t111, t107 * t117 - t108 * t112, t108 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1, t103 * t112 - t104 * t111, -t104 * t117 - t105 * t112, t103 * t117 + t105 * t111, t105 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1;];
T_reg = t1;
