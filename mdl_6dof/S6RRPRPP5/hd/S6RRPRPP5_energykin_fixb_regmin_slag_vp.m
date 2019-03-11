% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:16
% EndTime: 2019-03-09 10:06:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (267->47), mult. (529->94), div. (0->0), fcn. (271->4), ass. (0->36)
t168 = -pkin(2) - pkin(8);
t167 = -pkin(4) - pkin(5);
t153 = qJD(1) ^ 2;
t166 = t153 / 0.2e1;
t152 = cos(qJ(2));
t165 = t152 * t153;
t150 = sin(qJ(2));
t157 = -qJ(3) * t150 - pkin(1);
t134 = (t168 * t152 + t157) * qJD(1);
t162 = t150 * qJD(1);
t161 = pkin(7) * t162 + qJD(3);
t135 = pkin(3) * t162 + t168 * qJD(2) + t161;
t149 = sin(qJ(4));
t151 = cos(qJ(4));
t164 = t151 * t134 + t149 * t135;
t163 = qJD(1) * t152;
t141 = -pkin(7) * t163 - qJD(2) * qJ(3);
t160 = qJD(1) * qJD(2);
t144 = qJD(4) + t162;
t129 = t144 * qJ(5) + t164;
t136 = pkin(3) * t163 - t141;
t159 = t150 * t160;
t158 = t152 * t160;
t156 = -t149 * t134 + t151 * t135;
t155 = qJD(5) - t156;
t139 = t151 * qJD(2) - t149 * t163;
t154 = t139 * qJ(5) - t136;
t140 = -qJD(2) * pkin(2) + t161;
t138 = t149 * qJD(2) + t151 * t163;
t137 = (-pkin(2) * t152 + t157) * qJD(1);
t130 = t138 * pkin(4) - t154;
t128 = -t144 * pkin(4) + t155;
t127 = t167 * t138 + qJD(6) + t154;
t126 = t138 * qJ(6) + t129;
t125 = -t139 * qJ(6) + t167 * t144 + t155;
t1 = [t166, 0, 0, t150 ^ 2 * t166, t150 * t165, t159, t158, qJD(2) ^ 2 / 0.2e1, pkin(1) * t165 - pkin(7) * t159, -t153 * pkin(1) * t150 - pkin(7) * t158 (t140 * t150 - t141 * t152) * qJD(1), t140 * qJD(2) + t137 * t163, -t141 * qJD(2) - t137 * t162, t137 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t144, -t138 * t144, t144 ^ 2 / 0.2e1, t136 * t138 + t156 * t144, t136 * t139 - t164 * t144, -t128 * t144 + t130 * t138, t128 * t139 - t129 * t138, t129 * t144 - t130 * t139, t129 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1, -t125 * t144 - t127 * t138, t126 * t144 + t127 * t139, -t125 * t139 + t126 * t138, t126 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1;];
T_reg  = t1;
