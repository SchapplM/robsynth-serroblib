% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:11
% EndTime: 2019-03-09 04:52:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (236->42), mult. (429->84), div. (0->0), fcn. (219->4), ass. (0->30)
t145 = -pkin(4) - pkin(5);
t134 = qJD(1) ^ 2;
t144 = t134 / 0.2e1;
t143 = cos(qJ(4));
t142 = t134 * qJ(2);
t132 = sin(qJ(3));
t133 = cos(qJ(3));
t121 = (pkin(3) * t132 - pkin(8) * t133 + qJ(2)) * qJD(1);
t127 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t122 = qJD(3) * pkin(8) + t132 * t127;
t131 = sin(qJ(4));
t141 = t131 * t121 + t143 * t122;
t140 = qJD(1) * t133;
t139 = qJD(3) * t127;
t138 = qJD(1) * qJD(3);
t128 = t132 * qJD(1) + qJD(4);
t116 = t128 * qJ(5) + t141;
t137 = t143 * t121 - t131 * t122;
t123 = -qJD(3) * pkin(3) - t133 * t127;
t136 = qJD(5) - t137;
t125 = t131 * qJD(3) + t143 * t140;
t135 = t125 * qJ(5) - t123;
t129 = -qJD(1) * pkin(1) + qJD(2);
t124 = -t143 * qJD(3) + t131 * t140;
t117 = t124 * pkin(4) - t135;
t115 = -t128 * pkin(4) + t136;
t114 = t145 * t124 + qJD(6) + t135;
t113 = t124 * qJ(6) + t116;
t112 = -t125 * qJ(6) + t145 * t128 + t136;
t1 = [t144, 0, 0, t129 * qJD(1), t142, qJ(2) ^ 2 * t144 + t129 ^ 2 / 0.2e1, t133 ^ 2 * t144, -t133 * t134 * t132, t133 * t138, -t132 * t138, qJD(3) ^ 2 / 0.2e1, t132 * t142 + t133 * t139, -t132 * t139 + t133 * t142, t125 ^ 2 / 0.2e1, -t125 * t124, t125 * t128, -t124 * t128, t128 ^ 2 / 0.2e1, t123 * t124 + t137 * t128, t123 * t125 - t141 * t128, -t115 * t128 + t117 * t124, t115 * t125 - t116 * t124, t116 * t128 - t117 * t125, t116 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1, -t112 * t128 - t114 * t124, t113 * t128 + t114 * t125, -t112 * t125 + t113 * t124, t113 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1;];
T_reg  = t1;
