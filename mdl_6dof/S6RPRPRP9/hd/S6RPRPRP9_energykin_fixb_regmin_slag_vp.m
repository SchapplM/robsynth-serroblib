% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP9
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:27
% EndTime: 2019-03-09 03:29:27
% DurationCPUTime: 0.10s
% Computational Cost: add. (328->45), mult. (640->92), div. (0->0), fcn. (377->6), ass. (0->34)
t157 = qJD(1) ^ 2;
t165 = t157 / 0.2e1;
t164 = t157 * qJ(2);
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t143 = (pkin(3) * t154 - qJ(4) * t156 + qJ(2)) * qJD(1);
t147 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t144 = qJD(3) * qJ(4) + t154 * t147;
t151 = sin(pkin(9));
t152 = cos(pkin(9));
t133 = t152 * t143 - t151 * t144;
t162 = qJD(1) * t156;
t146 = t151 * qJD(3) + t152 * t162;
t160 = t154 * qJD(1);
t130 = pkin(4) * t160 - t146 * pkin(8) + t133;
t134 = t151 * t143 + t152 * t144;
t145 = -t152 * qJD(3) + t151 * t162;
t132 = -t145 * pkin(8) + t134;
t153 = sin(qJ(5));
t155 = cos(qJ(5));
t163 = t153 * t130 + t155 * t132;
t161 = qJD(3) * t147;
t159 = qJD(1) * qJD(3);
t158 = t155 * t130 - t153 * t132;
t142 = -qJD(3) * pkin(3) - t156 * t147 + qJD(4);
t137 = t145 * pkin(4) + t142;
t149 = -qJD(1) * pkin(1) + qJD(2);
t148 = qJD(5) + t160;
t136 = -t153 * t145 + t155 * t146;
t135 = t155 * t145 + t153 * t146;
t128 = t135 * pkin(5) - t136 * qJ(6) + t137;
t127 = t148 * qJ(6) + t163;
t126 = -t148 * pkin(5) + qJD(6) - t158;
t1 = [t165, 0, 0, t149 * qJD(1), t164, qJ(2) ^ 2 * t165 + t149 ^ 2 / 0.2e1, t156 ^ 2 * t165, -t156 * t157 * t154, t156 * t159, -t154 * t159, qJD(3) ^ 2 / 0.2e1, t154 * t164 + t156 * t161, -t154 * t161 + t156 * t164, t133 * t160 + t142 * t145, -t134 * t160 + t142 * t146, -t133 * t146 - t134 * t145, t134 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t136 ^ 2 / 0.2e1, -t136 * t135, t136 * t148, -t135 * t148, t148 ^ 2 / 0.2e1, t137 * t135 + t158 * t148, t137 * t136 - t163 * t148, -t126 * t148 + t128 * t135, t126 * t136 - t127 * t135, t127 * t148 - t128 * t136, t127 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1;];
T_reg  = t1;
