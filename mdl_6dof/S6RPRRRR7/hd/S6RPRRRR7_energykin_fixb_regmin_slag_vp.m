% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x34]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:54
% EndTime: 2019-03-09 07:17:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (320->46), mult. (655->97), div. (0->0), fcn. (451->8), ass. (0->39)
t157 = qJD(1) ^ 2;
t166 = t157 / 0.2e1;
t165 = t157 * qJ(2);
t151 = sin(qJ(4));
t152 = sin(qJ(3));
t155 = cos(qJ(4));
t156 = cos(qJ(3));
t140 = (-t151 * t152 + t155 * t156) * qJD(1);
t147 = qJD(3) + qJD(4);
t143 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t160 = -pkin(8) * qJD(1) + t143;
t136 = qJD(3) * pkin(3) + t160 * t156;
t138 = t160 * t152;
t159 = t155 * t136 - t151 * t138;
t124 = t147 * pkin(4) - t140 * pkin(9) + t159;
t139 = (t151 * t156 + t152 * t155) * qJD(1);
t163 = t151 * t136 + t155 * t138;
t126 = -t139 * pkin(9) + t163;
t150 = sin(qJ(5));
t154 = cos(qJ(5));
t164 = t150 * t124 + t154 * t126;
t141 = (pkin(3) * t152 + qJ(2)) * qJD(1);
t162 = qJD(3) * t143;
t161 = qJD(1) * qJD(3);
t130 = t154 * t139 + t150 * t140;
t132 = t139 * pkin(4) + t141;
t158 = t154 * t124 - t150 * t126;
t153 = cos(qJ(6));
t149 = sin(qJ(6));
t146 = qJD(5) + t147;
t145 = -qJD(1) * pkin(1) + qJD(2);
t131 = -t150 * t139 + t154 * t140;
t129 = qJD(6) + t130;
t128 = t153 * t131 + t149 * t146;
t127 = t149 * t131 - t153 * t146;
t122 = t130 * pkin(5) - t131 * pkin(10) + t132;
t121 = t146 * pkin(10) + t164;
t120 = -t146 * pkin(5) - t158;
t1 = [t166, 0, 0, t145 * qJD(1), t165, qJ(2) ^ 2 * t166 + t145 ^ 2 / 0.2e1, t156 ^ 2 * t166, -t156 * t157 * t152, t156 * t161, -t152 * t161, qJD(3) ^ 2 / 0.2e1, t152 * t165 + t156 * t162, -t152 * t162 + t156 * t165, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t147, -t139 * t147, t147 ^ 2 / 0.2e1, t141 * t139 + t159 * t147, t141 * t140 - t163 * t147, t131 ^ 2 / 0.2e1, -t131 * t130, t131 * t146, -t130 * t146, t146 ^ 2 / 0.2e1, t132 * t130 + t158 * t146, t132 * t131 - t164 * t146, t128 ^ 2 / 0.2e1, -t128 * t127, t128 * t129, -t127 * t129, t129 ^ 2 / 0.2e1 (-t149 * t121 + t153 * t122) * t129 + t120 * t127 -(t153 * t121 + t149 * t122) * t129 + t120 * t128;];
T_reg  = t1;
