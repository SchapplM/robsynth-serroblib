% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:38
% EndTime: 2019-03-09 02:18:38
% DurationCPUTime: 0.12s
% Computational Cost: add. (300->50), mult. (728->103), div. (0->0), fcn. (525->10), ass. (0->41)
t179 = cos(qJ(4));
t161 = sin(pkin(11));
t163 = cos(pkin(11));
t167 = sin(qJ(4));
t152 = (t179 * t161 + t163 * t167) * qJD(1);
t162 = sin(pkin(10));
t155 = (pkin(1) * t162 + qJ(3)) * qJD(1);
t159 = t163 * qJD(2);
t145 = t159 + (-pkin(7) * qJD(1) - t155) * t161;
t148 = t161 * qJD(2) + t163 * t155;
t175 = qJD(1) * t163;
t146 = pkin(7) * t175 + t148;
t173 = t179 * t145 - t167 * t146;
t134 = qJD(4) * pkin(4) - t152 * pkin(8) + t173;
t176 = qJD(1) * t161;
t151 = t167 * t176 - t179 * t175;
t177 = t167 * t145 + t179 * t146;
t135 = -t151 * pkin(8) + t177;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t178 = t166 * t134 + t169 * t135;
t164 = cos(pkin(10));
t174 = -pkin(1) * t164 - pkin(2);
t139 = t169 * t151 + t166 * t152;
t172 = t169 * t134 - t166 * t135;
t150 = qJD(3) + (-pkin(3) * t163 + t174) * qJD(1);
t141 = t151 * pkin(4) + t150;
t170 = qJD(1) ^ 2;
t168 = cos(qJ(6));
t165 = sin(qJ(6));
t160 = qJD(4) + qJD(5);
t154 = t174 * qJD(1) + qJD(3);
t147 = -t161 * t155 + t159;
t140 = -t166 * t151 + t169 * t152;
t138 = qJD(6) + t139;
t137 = t168 * t140 + t165 * t160;
t136 = t165 * t140 - t168 * t160;
t131 = t139 * pkin(5) - t140 * pkin(9) + t141;
t130 = t160 * pkin(9) + t178;
t129 = -t160 * pkin(5) - t172;
t1 = [t170 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t162 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t170, -t154 * t175, t154 * t176 (-t147 * t161 + t148 * t163) * qJD(1), t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * qJD(4), -t151 * qJD(4), qJD(4) ^ 2 / 0.2e1, t173 * qJD(4) + t150 * t151, -t177 * qJD(4) + t150 * t152, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t160, -t139 * t160, t160 ^ 2 / 0.2e1, t141 * t139 + t172 * t160, t141 * t140 - t178 * t160, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * t138, -t136 * t138, t138 ^ 2 / 0.2e1 (-t165 * t130 + t168 * t131) * t138 + t129 * t136 -(t168 * t130 + t165 * t131) * t138 + t129 * t137;];
T_reg  = t1;
