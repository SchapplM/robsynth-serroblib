% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:31
% EndTime: 2019-03-09 03:32:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (195->44), mult. (352->88), div. (0->0), fcn. (155->4), ass. (0->29)
t137 = qJD(1) ^ 2;
t146 = t137 / 0.2e1;
t145 = t137 * qJ(2);
t126 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t136 = cos(qJ(3));
t116 = qJD(4) + (pkin(4) * qJD(1) - t126) * t136 + (-pkin(3) - pkin(8)) * qJD(3);
t134 = sin(qJ(3));
t142 = qJD(1) * t134;
t143 = pkin(3) * t142 + qJD(1) * qJ(2);
t118 = (pkin(8) * t134 - qJ(4) * t136) * qJD(1) + t143;
t133 = sin(qJ(5));
t135 = cos(qJ(5));
t144 = t133 * t116 + t135 * t118;
t122 = -qJD(3) * qJ(4) - t134 * t126;
t141 = qJD(3) * t126;
t140 = t136 * qJD(1);
t139 = qJD(1) * qJD(3);
t138 = t135 * t116 - t133 * t118;
t119 = -pkin(4) * t142 - t122;
t130 = -qJD(1) * pkin(1) + qJD(2);
t128 = qJD(5) + t140;
t124 = t135 * qJD(3) + t133 * t142;
t123 = t133 * qJD(3) - t135 * t142;
t121 = -qJ(4) * t140 + t143;
t120 = -qJD(3) * pkin(3) - t136 * t126 + qJD(4);
t114 = t123 * pkin(5) - t124 * qJ(6) + t119;
t113 = t128 * qJ(6) + t144;
t112 = -t128 * pkin(5) + qJD(6) - t138;
t1 = [t146, 0, 0, t130 * qJD(1), t145, qJ(2) ^ 2 * t146 + t130 ^ 2 / 0.2e1, t136 ^ 2 * t146, -t136 * t137 * t134, t136 * t139, -t134 * t139, qJD(3) ^ 2 / 0.2e1, t134 * t145 + t136 * t141, -t134 * t141 + t136 * t145 (t120 * t136 + t122 * t134) * qJD(1), t120 * qJD(3) - t121 * t142, -t122 * qJD(3) - t121 * t140, t121 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t128, -t123 * t128, t128 ^ 2 / 0.2e1, t119 * t123 + t138 * t128, t119 * t124 - t144 * t128, -t112 * t128 + t114 * t123, t112 * t124 - t113 * t123, t113 * t128 - t114 * t124, t113 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1;];
T_reg  = t1;
