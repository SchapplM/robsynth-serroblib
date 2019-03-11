% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:32
% EndTime: 2019-03-09 08:35:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (172->45), mult. (359->91), div. (0->0), fcn. (159->4), ass. (0->33)
t162 = -pkin(2) - pkin(3);
t150 = qJD(1) ^ 2;
t161 = t150 / 0.2e1;
t149 = cos(qJ(2));
t160 = t149 * t150;
t147 = sin(qJ(2));
t157 = t147 * qJD(1);
t158 = qJD(1) * t149;
t131 = -qJD(1) * pkin(1) - pkin(2) * t158 - qJ(3) * t157;
t127 = pkin(3) * t158 + qJD(4) - t131;
t124 = (pkin(4) * t147 + pkin(8) * t149) * qJD(1) + t127;
t156 = pkin(7) * t157 + qJD(3);
t151 = -qJ(4) * t157 + t156;
t126 = (-pkin(8) + t162) * qJD(2) + t151;
t146 = sin(qJ(5));
t148 = cos(qJ(5));
t159 = t146 * t124 + t148 * t126;
t136 = pkin(7) * t158 + qJD(2) * qJ(3);
t155 = qJD(1) * qJD(2);
t154 = t147 * t155;
t153 = t149 * t155;
t152 = t148 * t124 - t126 * t146;
t130 = qJ(4) * t158 - t136;
t129 = qJD(2) * pkin(4) - t130;
t137 = qJD(5) + t157;
t135 = -qJD(2) * pkin(2) + t156;
t133 = t146 * qJD(2) + t148 * t158;
t132 = -t148 * qJD(2) + t146 * t158;
t128 = t162 * qJD(2) + t151;
t123 = -t132 * pkin(5) + qJD(6) + t129;
t120 = t132 * qJ(6) + t159;
t119 = pkin(5) * t137 + qJ(6) * t133 + t152;
t1 = [t161, 0, 0, t147 ^ 2 * t161, t147 * t160, t154, t153, qJD(2) ^ 2 / 0.2e1, pkin(1) * t160 - pkin(7) * t154, -pkin(1) * t147 * t150 - pkin(7) * t153, -qJD(2) * t135 - t131 * t158 (t135 * t147 + t136 * t149) * qJD(1), qJD(2) * t136 - t131 * t157, t136 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, -qJD(2) * t130 + t127 * t157, qJD(2) * t128 - t127 * t158 (-t128 * t147 + t130 * t149) * qJD(1), t128 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1, t133 ^ 2 / 0.2e1, -t133 * t132, -t133 * t137, t132 * t137, t137 ^ 2 / 0.2e1, -t129 * t132 + t152 * t137, -t129 * t133 - t159 * t137, t119 * t133 + t120 * t132, t120 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1;];
T_reg  = t1;
