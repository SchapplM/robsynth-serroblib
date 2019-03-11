% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:46
% EndTime: 2019-03-09 02:26:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (208->42), mult. (391->91), div. (0->0), fcn. (214->8), ass. (0->37)
t154 = qJD(1) ^ 2;
t163 = t154 / 0.2e1;
t139 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t145 = sin(pkin(10));
t146 = cos(pkin(10));
t160 = qJ(2) * qJD(1);
t133 = t145 * t139 + t146 * t160;
t131 = -qJD(1) * pkin(7) + t133;
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t161 = t150 * qJD(3) + t153 * t131;
t124 = qJD(4) * pkin(8) + t161;
t132 = t146 * t139 - t145 * t160;
t130 = qJD(1) * pkin(3) - t132;
t125 = (pkin(4) * t153 + pkin(8) * t150) * qJD(1) + t130;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t162 = t152 * t124 + t149 * t125;
t159 = qJD(1) * t150;
t158 = t153 * qJD(1);
t157 = qJD(1) * qJD(4);
t156 = -t149 * t124 + t152 * t125;
t155 = t153 * qJD(3) - t150 * t131;
t140 = qJD(5) + t158;
t123 = -qJD(4) * pkin(4) - t155;
t151 = cos(qJ(6));
t148 = sin(qJ(6));
t143 = -qJD(1) * pkin(1) + qJD(2);
t138 = qJD(6) + t140;
t136 = t149 * qJD(4) - t152 * t159;
t135 = t152 * qJD(4) + t149 * t159;
t127 = t148 * t135 + t151 * t136;
t126 = -t151 * t135 + t148 * t136;
t119 = -t135 * pkin(5) + t123;
t118 = t135 * pkin(9) + t162;
t117 = t140 * pkin(5) - t136 * pkin(9) + t156;
t1 = [t163, 0, 0, -t143 * qJD(1), t154 * qJ(2), qJ(2) ^ 2 * t163 + t143 ^ 2 / 0.2e1, -t132 * qJD(1), t133 * qJD(1), t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t150 ^ 2 * t163, t150 * t154 * t153, -t150 * t157, -t153 * t157, qJD(4) ^ 2 / 0.2e1, t155 * qJD(4) + t130 * t158, -t161 * qJD(4) - t130 * t159, t136 ^ 2 / 0.2e1, t136 * t135, t136 * t140, t135 * t140, t140 ^ 2 / 0.2e1, -t123 * t135 + t156 * t140, t123 * t136 - t162 * t140, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t138, -t126 * t138, t138 ^ 2 / 0.2e1 (t151 * t117 - t148 * t118) * t138 + t119 * t126 -(t148 * t117 + t151 * t118) * t138 + t119 * t127;];
T_reg  = t1;
