% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:48
% EndTime: 2019-03-09 02:23:48
% DurationCPUTime: 0.08s
% Computational Cost: add. (164->42), mult. (338->87), div. (0->0), fcn. (192->8), ass. (0->37)
t148 = qJD(1) ^ 2;
t161 = t148 / 0.2e1;
t160 = cos(qJ(5));
t142 = cos(pkin(10));
t153 = -pkin(1) * t142 - pkin(2);
t129 = qJD(3) + (-pkin(7) + t153) * qJD(1);
t145 = sin(qJ(4));
t147 = cos(qJ(4));
t158 = t147 * qJD(2) + t145 * t129;
t123 = qJD(4) * pkin(8) + t158;
t141 = sin(pkin(10));
t152 = -pkin(1) * t141 - qJ(3);
t126 = (pkin(4) * t145 - pkin(8) * t147 - t152) * qJD(1);
t144 = sin(qJ(5));
t159 = t160 * t123 + t144 * t126;
t157 = qJD(1) * t147;
t133 = t152 * qJD(1);
t156 = t133 * qJD(1);
t155 = t145 * qJD(1);
t154 = qJD(1) * qJD(4);
t151 = -t144 * t123 + t160 * t126;
t150 = -t145 * qJD(2) + t147 * t129;
t136 = qJD(5) + t155;
t122 = -qJD(4) * pkin(4) - t150;
t146 = cos(qJ(6));
t143 = sin(qJ(6));
t140 = qJD(2) ^ 2 / 0.2e1;
t135 = qJD(6) + t136;
t132 = t153 * qJD(1) + qJD(3);
t131 = t144 * qJD(4) + t160 * t157;
t130 = -t160 * qJD(4) + t144 * t157;
t120 = -t143 * t130 + t146 * t131;
t119 = t146 * t130 + t143 * t131;
t118 = t130 * pkin(5) + t122;
t117 = -t130 * pkin(9) + t159;
t116 = t136 * pkin(5) - t131 * pkin(9) + t151;
t1 = [t161, 0, 0, t140 + (t141 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t148, t132 * qJD(1), -t156, t140 + t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t147 ^ 2 * t161, -t147 * t148 * t145, t147 * t154, -t145 * t154, qJD(4) ^ 2 / 0.2e1, t150 * qJD(4) - t133 * t155, -t158 * qJD(4) - t147 * t156, t131 ^ 2 / 0.2e1, -t131 * t130, t131 * t136, -t130 * t136, t136 ^ 2 / 0.2e1, t122 * t130 + t151 * t136, t122 * t131 - t159 * t136, t120 ^ 2 / 0.2e1, -t120 * t119, t120 * t135, -t119 * t135, t135 ^ 2 / 0.2e1 (t146 * t116 - t143 * t117) * t135 + t118 * t119 -(t143 * t116 + t146 * t117) * t135 + t118 * t120;];
T_reg  = t1;
