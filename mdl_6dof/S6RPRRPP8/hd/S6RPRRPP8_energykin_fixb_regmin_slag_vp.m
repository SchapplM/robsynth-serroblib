% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP8
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:52
% EndTime: 2019-03-09 04:55:52
% DurationCPUTime: 0.07s
% Computational Cost: add. (236->42), mult. (429->84), div. (0->0), fcn. (219->4), ass. (0->30)
t132 = qJD(1) ^ 2;
t142 = t132 / 0.2e1;
t141 = pkin(4) + qJ(6);
t140 = t132 * qJ(2);
t129 = sin(qJ(3));
t131 = cos(qJ(3));
t119 = (pkin(3) * t129 - pkin(8) * t131 + qJ(2)) * qJD(1);
t124 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t120 = qJD(3) * pkin(8) + t129 * t124;
t128 = sin(qJ(4));
t130 = cos(qJ(4));
t139 = t128 * t119 + t130 * t120;
t138 = qJD(1) * t131;
t137 = qJD(3) * t124;
t136 = qJD(1) * qJD(3);
t135 = t130 * t119 - t128 * t120;
t125 = t129 * qJD(1) + qJD(4);
t114 = -t125 * qJ(5) - t139;
t134 = qJD(5) - t135;
t121 = -qJD(3) * pkin(3) - t131 * t124;
t123 = t128 * qJD(3) + t130 * t138;
t133 = -t123 * qJ(5) + t121;
t126 = -qJD(1) * pkin(1) + qJD(2);
t122 = -t130 * qJD(3) + t128 * t138;
t115 = t122 * pkin(4) + t133;
t113 = -t125 * pkin(4) + t134;
t112 = t141 * t122 + t133;
t111 = -t122 * pkin(5) + qJD(6) - t114;
t110 = t123 * pkin(5) - t141 * t125 + t134;
t1 = [t142, 0, 0, t126 * qJD(1), t140, qJ(2) ^ 2 * t142 + t126 ^ 2 / 0.2e1, t131 ^ 2 * t142, -t131 * t132 * t129, t131 * t136, -t129 * t136, qJD(3) ^ 2 / 0.2e1, t129 * t140 + t131 * t137, -t129 * t137 + t131 * t140, t123 ^ 2 / 0.2e1, -t123 * t122, t123 * t125, -t122 * t125, t125 ^ 2 / 0.2e1, t121 * t122 + t135 * t125, t121 * t123 - t139 * t125, t113 * t123 + t114 * t122, t113 * t125 - t115 * t122, -t114 * t125 - t115 * t123, t115 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1, t110 * t123 - t111 * t122, t111 * t125 - t112 * t123, -t110 * t125 + t112 * t122, t112 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1;];
T_reg  = t1;
