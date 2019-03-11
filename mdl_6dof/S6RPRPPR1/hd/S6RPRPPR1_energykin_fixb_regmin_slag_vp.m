% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:34
% EndTime: 2019-03-09 02:39:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (356->48), mult. (809->104), div. (0->0), fcn. (535->10), ass. (0->41)
t185 = qJD(1) ^ 2;
t193 = t185 / 0.2e1;
t192 = cos(pkin(11));
t178 = sin(pkin(9));
t170 = (pkin(1) * t178 + pkin(7)) * qJD(1);
t184 = cos(qJ(3));
t175 = t184 * qJD(2);
t182 = sin(qJ(3));
t161 = qJD(3) * pkin(3) + t175 + (-qJ(4) * qJD(1) - t170) * t182;
t189 = qJD(1) * t184;
t191 = t182 * qJD(2) + t184 * t170;
t164 = qJ(4) * t189 + t191;
t177 = sin(pkin(10));
t179 = cos(pkin(10));
t151 = t177 * t161 + t179 * t164;
t149 = qJD(3) * qJ(5) + t151;
t180 = cos(pkin(9));
t187 = -pkin(1) * t180 - pkin(2);
t166 = qJD(4) + (-pkin(3) * t184 + t187) * qJD(1);
t190 = qJD(1) * t182;
t167 = t177 * t190 - t179 * t189;
t168 = (t177 * t184 + t179 * t182) * qJD(1);
t156 = t167 * pkin(4) - t168 * qJ(5) + t166;
t176 = sin(pkin(11));
t145 = t192 * t149 + t176 * t156;
t188 = qJD(1) * qJD(3);
t144 = -t176 * t149 + t192 * t156;
t150 = t179 * t161 - t177 * t164;
t148 = -qJD(3) * pkin(4) + qJD(5) - t150;
t183 = cos(qJ(6));
t181 = sin(qJ(6));
t171 = t187 * qJD(1);
t165 = qJD(6) + t167;
t163 = t176 * qJD(3) + t192 * t168;
t162 = -t192 * qJD(3) + t176 * t168;
t153 = -t181 * t162 + t183 * t163;
t152 = t183 * t162 + t181 * t163;
t146 = t162 * pkin(5) + t148;
t143 = -t162 * pkin(8) + t145;
t142 = t167 * pkin(5) - t163 * pkin(8) + t144;
t1 = [t193, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t178 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t185, t182 ^ 2 * t193, t182 * t185 * t184, t182 * t188, t184 * t188, qJD(3) ^ 2 / 0.2e1, -t171 * t189 + (-t182 * t170 + t175) * qJD(3), -t191 * qJD(3) + t171 * t190, -t150 * t168 - t151 * t167, t151 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t144 * t167 + t148 * t162, -t145 * t167 + t148 * t163, -t144 * t163 - t145 * t162, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t165, -t152 * t165, t165 ^ 2 / 0.2e1 (t183 * t142 - t181 * t143) * t165 + t146 * t152 -(t181 * t142 + t183 * t143) * t165 + t146 * t153;];
T_reg  = t1;
