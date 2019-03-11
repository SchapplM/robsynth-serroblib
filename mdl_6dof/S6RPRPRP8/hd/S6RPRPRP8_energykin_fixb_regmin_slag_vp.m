% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:03
% EndTime: 2019-03-09 03:26:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (305->42), mult. (602->87), div. (0->0), fcn. (362->6), ass. (0->33)
t151 = sin(pkin(9));
t152 = cos(pkin(9));
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t143 = (t151 * t156 + t152 * t154) * qJD(1);
t157 = qJD(1) ^ 2;
t165 = t157 / 0.2e1;
t164 = t157 * qJ(2);
t146 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t160 = -qJ(4) * qJD(1) + t146;
t140 = qJD(3) * pkin(3) + t160 * t156;
t141 = t160 * t154;
t133 = t151 * t140 + t152 * t141;
t130 = qJD(3) * pkin(8) + t133;
t144 = (-t151 * t154 + t152 * t156) * qJD(1);
t145 = qJD(4) + (pkin(3) * t154 + qJ(2)) * qJD(1);
t134 = pkin(4) * t143 - pkin(8) * t144 + t145;
t153 = sin(qJ(5));
t155 = cos(qJ(5));
t163 = t155 * t130 + t153 * t134;
t162 = qJD(3) * t146;
t161 = qJD(1) * qJD(3);
t132 = t140 * t152 - t151 * t141;
t159 = -t130 * t153 + t134 * t155;
t129 = -qJD(3) * pkin(4) - t132;
t148 = -qJD(1) * pkin(1) + qJD(2);
t142 = qJD(5) + t143;
t136 = qJD(3) * t153 + t144 * t155;
t135 = -t155 * qJD(3) + t144 * t153;
t127 = pkin(5) * t135 - qJ(6) * t136 + t129;
t126 = qJ(6) * t142 + t163;
t125 = -pkin(5) * t142 + qJD(6) - t159;
t1 = [t165, 0, 0, t148 * qJD(1), t164, qJ(2) ^ 2 * t165 + t148 ^ 2 / 0.2e1, t156 ^ 2 * t165, -t156 * t157 * t154, t156 * t161, -t154 * t161, qJD(3) ^ 2 / 0.2e1, t154 * t164 + t156 * t162, -t154 * t162 + t156 * t164, -t132 * t144 - t133 * t143, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, t136 ^ 2 / 0.2e1, -t136 * t135, t136 * t142, -t135 * t142, t142 ^ 2 / 0.2e1, t129 * t135 + t159 * t142, t129 * t136 - t163 * t142, -t125 * t142 + t127 * t135, t125 * t136 - t126 * t135, t126 * t142 - t127 * t136, t126 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1;];
T_reg  = t1;
