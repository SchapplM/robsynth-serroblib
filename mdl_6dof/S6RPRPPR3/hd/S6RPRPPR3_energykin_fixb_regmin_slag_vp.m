% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (149->41), mult. (323->92), div. (0->0), fcn. (146->6), ass. (0->31)
t149 = -pkin(3) - pkin(4);
t139 = qJD(1) ^ 2;
t148 = t139 / 0.2e1;
t132 = sin(pkin(9));
t123 = (pkin(1) * t132 + pkin(7)) * qJD(1);
t136 = sin(qJ(3));
t138 = cos(qJ(3));
t147 = t136 * qJD(2) + t138 * t123;
t133 = cos(pkin(9));
t124 = (-pkin(1) * t133 - pkin(2)) * qJD(1);
t146 = qJD(1) * t138;
t145 = t136 * qJD(1);
t144 = qJD(1) * qJD(3);
t117 = qJD(3) * qJ(4) + t147;
t143 = t138 * qJD(2) - t136 * t123;
t118 = -pkin(3) * t146 - qJ(4) * t145 + t124;
t142 = qJD(4) - t143;
t114 = pkin(4) * t146 + qJD(5) - t118;
t115 = qJ(5) * t146 - t117;
t141 = -qJ(5) * t145 + t142;
t137 = cos(qJ(6));
t135 = sin(qJ(6));
t125 = qJD(6) + t145;
t122 = -t135 * qJD(3) - t137 * t146;
t121 = -t137 * qJD(3) + t135 * t146;
t116 = -qJD(3) * pkin(3) + t142;
t113 = qJD(3) * pkin(5) - t115;
t112 = t149 * qJD(3) + t141;
t111 = (-pkin(8) + t149) * qJD(3) + t141;
t110 = (pkin(5) * t136 + pkin(8) * t138) * qJD(1) + t114;
t1 = [t148, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t132 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t139, t136 ^ 2 * t148, t136 * t139 * t138, t136 * t144, t138 * t144, qJD(3) ^ 2 / 0.2e1, t143 * qJD(3) - t124 * t146, -t147 * qJD(3) + t124 * t145, -t116 * qJD(3) - t118 * t146 (t116 * t136 + t117 * t138) * qJD(1), t117 * qJD(3) - t118 * t145, t117 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t116 ^ 2 / 0.2e1, -t115 * qJD(3) + t114 * t145, t112 * qJD(3) - t114 * t146 (-t112 * t136 + t115 * t138) * qJD(1), t112 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1 + t114 ^ 2 / 0.2e1, t122 ^ 2 / 0.2e1, t122 * t121, t122 * t125, t121 * t125, t125 ^ 2 / 0.2e1 (t137 * t110 - t135 * t111) * t125 - t113 * t121 -(t135 * t110 + t137 * t111) * t125 + t113 * t122;];
T_reg  = t1;
