% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:00
% EndTime: 2019-03-09 03:13:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (205->42), mult. (423->92), div. (0->0), fcn. (212->6), ass. (0->34)
t163 = -pkin(3) - pkin(8);
t150 = qJD(1) ^ 2;
t162 = t150 / 0.2e1;
t144 = sin(pkin(9));
t138 = (pkin(1) * t144 + pkin(7)) * qJD(1);
t147 = sin(qJ(3));
t149 = cos(qJ(3));
t155 = t149 * qJD(2) - t147 * t138;
t154 = qJD(4) - t155;
t158 = t147 * qJD(1);
t127 = pkin(4) * t158 + t163 * qJD(3) + t154;
t145 = cos(pkin(9));
t156 = -pkin(1) * t145 - pkin(2);
t152 = -qJ(4) * t147 + t156;
t130 = (t163 * t149 + t152) * qJD(1);
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t161 = t146 * t127 + t148 * t130;
t160 = t147 * qJD(2) + t149 * t138;
t159 = qJD(1) * t149;
t157 = qJD(1) * qJD(3);
t132 = -qJD(3) * qJ(4) - t160;
t129 = pkin(4) * t159 - t132;
t153 = t148 * t127 - t146 * t130;
t140 = qJD(5) + t158;
t139 = t156 * qJD(1);
t137 = t148 * qJD(3) - t146 * t159;
t136 = t146 * qJD(3) + t148 * t159;
t133 = (-pkin(3) * t149 + t152) * qJD(1);
t131 = -qJD(3) * pkin(3) + t154;
t125 = t136 * pkin(5) - t137 * qJ(6) + t129;
t124 = t140 * qJ(6) + t161;
t123 = -t140 * pkin(5) + qJD(6) - t153;
t1 = [t162, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t150, t147 ^ 2 * t162, t147 * t150 * t149, t147 * t157, t149 * t157, qJD(3) ^ 2 / 0.2e1, t155 * qJD(3) - t139 * t159, -t160 * qJD(3) + t139 * t158 (t131 * t147 - t132 * t149) * qJD(1), t131 * qJD(3) + t133 * t159, -t132 * qJD(3) - t133 * t158, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * t140, -t136 * t140, t140 ^ 2 / 0.2e1, t129 * t136 + t153 * t140, t129 * t137 - t161 * t140, -t123 * t140 + t125 * t136, t123 * t137 - t124 * t136, t124 * t140 - t125 * t137, t124 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1;];
T_reg  = t1;
