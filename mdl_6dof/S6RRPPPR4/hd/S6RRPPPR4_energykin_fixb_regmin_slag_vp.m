% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:43
% EndTime: 2019-03-09 08:19:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (271->51), mult. (572->102), div. (0->0), fcn. (306->6), ass. (0->40)
t192 = -pkin(4) - pkin(5);
t179 = qJD(1) ^ 2;
t191 = t179 / 0.2e1;
t190 = -pkin(2) - qJ(4);
t178 = cos(qJ(2));
t189 = t178 * t179;
t176 = sin(qJ(2));
t182 = -qJ(3) * t176 - pkin(1);
t156 = (t190 * t178 + t182) * qJD(1);
t187 = t176 * qJD(1);
t186 = pkin(7) * t187 + qJD(3);
t157 = pkin(3) * t187 + t190 * qJD(2) + t186;
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t149 = t174 * t156 + t173 * t157;
t188 = qJD(1) * t178;
t164 = -pkin(7) * t188 - qJD(2) * qJ(3);
t185 = qJD(1) * qJD(2);
t147 = qJ(5) * t187 + t149;
t184 = t176 * t185;
t183 = t178 * t185;
t148 = -t173 * t156 + t174 * t157;
t159 = pkin(3) * t188 + qJD(4) - t164;
t181 = qJD(5) - t148;
t162 = t174 * qJD(2) - t173 * t188;
t180 = t162 * qJ(5) - t159;
t177 = cos(qJ(6));
t175 = sin(qJ(6));
t166 = -qJD(6) + t187;
t163 = -qJD(2) * pkin(2) + t186;
t161 = t173 * qJD(2) + t174 * t188;
t160 = (-pkin(2) * t178 + t182) * qJD(1);
t152 = t175 * t161 + t177 * t162;
t151 = -t177 * t161 + t175 * t162;
t150 = t161 * pkin(4) - t180;
t146 = -pkin(4) * t187 + t181;
t145 = t192 * t161 + t180;
t144 = t161 * pkin(8) + t147;
t143 = -t162 * pkin(8) + t192 * t187 + t181;
t1 = [t191, 0, 0, t176 ^ 2 * t191, t176 * t189, t184, t183, qJD(2) ^ 2 / 0.2e1, pkin(1) * t189 - pkin(7) * t184, -t179 * pkin(1) * t176 - pkin(7) * t183 (t163 * t176 - t164 * t178) * qJD(1), t163 * qJD(2) + t160 * t188, -t164 * qJD(2) - t160 * t187, t160 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t148 * t187 + t159 * t161, -t149 * t187 + t159 * t162, -t148 * t162 - t149 * t161, t149 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, -t146 * t187 + t150 * t161, t146 * t162 - t147 * t161, t147 * t187 - t150 * t162, t147 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t152 ^ 2 / 0.2e1, -t152 * t151, -t152 * t166, t151 * t166, t166 ^ 2 / 0.2e1 -(t177 * t143 - t175 * t144) * t166 + t145 * t151 (t175 * t143 + t177 * t144) * t166 + t145 * t152;];
T_reg  = t1;
