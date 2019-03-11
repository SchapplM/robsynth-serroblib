% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:23
% EndTime: 2019-03-09 02:16:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (286->41), mult. (591->89), div. (0->0), fcn. (368->6), ass. (0->34)
t149 = qJD(1) ^ 2;
t156 = t149 / 0.2e1;
t143 = sin(pkin(9));
t135 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t152 = -pkin(7) * qJD(1) + t135;
t128 = t152 * t143;
t144 = cos(pkin(9));
t129 = t152 * t144;
t146 = sin(qJ(4));
t148 = cos(qJ(4));
t154 = t148 * t128 + t146 * t129;
t121 = qJD(4) * pkin(8) + t154;
t131 = (t143 * t148 + t144 * t146) * qJD(1);
t132 = (-t143 * t146 + t144 * t148) * qJD(1);
t137 = qJD(1) * qJ(2) + qJD(3);
t153 = qJD(1) * t143;
t133 = pkin(3) * t153 + t137;
t122 = t131 * pkin(4) - t132 * pkin(8) + t133;
t145 = sin(qJ(5));
t147 = cos(qJ(5));
t155 = t147 * t121 + t145 * t122;
t151 = -t146 * t128 + t148 * t129;
t150 = -t145 * t121 + t147 * t122;
t120 = -qJD(4) * pkin(4) - t151;
t141 = t144 ^ 2;
t140 = t143 ^ 2;
t138 = -qJD(1) * pkin(1) + qJD(2);
t130 = qJD(5) + t131;
t124 = t145 * qJD(4) + t147 * t132;
t123 = -t147 * qJD(4) + t145 * t132;
t117 = t123 * pkin(5) - t124 * qJ(6) + t120;
t116 = t130 * qJ(6) + t155;
t115 = -t130 * pkin(5) + qJD(6) - t150;
t1 = [t156, 0, 0, t138 * qJD(1), t149 * qJ(2), qJ(2) ^ 2 * t156 + t138 ^ 2 / 0.2e1, t137 * t153, t137 * t144 * qJD(1) (-t140 - t141) * t135 * qJD(1), t137 ^ 2 / 0.2e1 + (t140 / 0.2e1 + t141 / 0.2e1) * t135 ^ 2, t132 ^ 2 / 0.2e1, -t132 * t131, t132 * qJD(4), -t131 * qJD(4), qJD(4) ^ 2 / 0.2e1, t151 * qJD(4) + t133 * t131, -t154 * qJD(4) + t133 * t132, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t130, -t123 * t130, t130 ^ 2 / 0.2e1, t120 * t123 + t150 * t130, t120 * t124 - t155 * t130, -t115 * t130 + t117 * t123, t115 * t124 - t116 * t123, t116 * t130 - t117 * t124, t116 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1;];
T_reg  = t1;
