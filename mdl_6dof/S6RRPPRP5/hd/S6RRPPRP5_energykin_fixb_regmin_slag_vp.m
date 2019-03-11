% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:56
% EndTime: 2019-03-09 08:43:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (375->50), mult. (772->102), div. (0->0), fcn. (445->6), ass. (0->39)
t177 = qJD(1) ^ 2;
t189 = t177 / 0.2e1;
t188 = -pkin(2) - qJ(4);
t176 = cos(qJ(2));
t187 = t176 * t177;
t174 = sin(qJ(2));
t179 = -qJ(3) * t174 - pkin(1);
t157 = (t188 * t176 + t179) * qJD(1);
t184 = t174 * qJD(1);
t183 = pkin(7) * t184 + qJD(3);
t158 = pkin(3) * t184 + t188 * qJD(2) + t183;
t171 = sin(pkin(9));
t172 = cos(pkin(9));
t149 = -t171 * t157 + t172 * t158;
t185 = qJD(1) * t176;
t163 = t172 * qJD(2) - t171 * t185;
t146 = pkin(4) * t184 - t163 * pkin(8) + t149;
t150 = t172 * t157 + t171 * t158;
t162 = t171 * qJD(2) + t172 * t185;
t148 = -t162 * pkin(8) + t150;
t173 = sin(qJ(5));
t175 = cos(qJ(5));
t186 = t173 * t146 + t175 * t148;
t165 = -pkin(7) * t185 - qJD(2) * qJ(3);
t182 = qJD(1) * qJD(2);
t181 = t174 * t182;
t180 = t176 * t182;
t160 = pkin(3) * t185 + qJD(4) - t165;
t178 = t175 * t146 - t173 * t148;
t153 = t162 * pkin(4) + t160;
t166 = qJD(5) + t184;
t164 = -qJD(2) * pkin(2) + t183;
t161 = (-pkin(2) * t176 + t179) * qJD(1);
t152 = -t173 * t162 + t175 * t163;
t151 = t175 * t162 + t173 * t163;
t144 = t151 * pkin(5) - t152 * qJ(6) + t153;
t143 = t166 * qJ(6) + t186;
t142 = -t166 * pkin(5) + qJD(6) - t178;
t1 = [t189, 0, 0, t174 ^ 2 * t189, t174 * t187, t181, t180, qJD(2) ^ 2 / 0.2e1, pkin(1) * t187 - pkin(7) * t181, -t177 * pkin(1) * t174 - pkin(7) * t180 (t164 * t174 - t165 * t176) * qJD(1), t164 * qJD(2) + t161 * t185, -t165 * qJD(2) - t161 * t184, t161 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t149 * t184 + t160 * t162, -t150 * t184 + t160 * t163, -t149 * t163 - t150 * t162, t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * t166, -t151 * t166, t166 ^ 2 / 0.2e1, t153 * t151 + t178 * t166, t153 * t152 - t186 * t166, -t142 * t166 + t144 * t151, t142 * t152 - t143 * t151, t143 * t166 - t144 * t152, t143 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
