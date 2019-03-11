% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:39
% EndTime: 2019-03-09 17:19:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (301->47), mult. (627->95), div. (0->0), fcn. (393->6), ass. (0->38)
t177 = qJD(1) ^ 2;
t190 = t177 / 0.2e1;
t189 = cos(qJ(5));
t176 = cos(qJ(2));
t188 = t176 * t177;
t173 = sin(qJ(3));
t175 = cos(qJ(3));
t174 = sin(qJ(2));
t185 = qJD(1) * t174;
t160 = t173 * qJD(2) + t175 * t185;
t184 = t176 * qJD(1);
t167 = -qJD(3) + t184;
t156 = (-pkin(2) * t176 - pkin(8) * t174 - pkin(1)) * qJD(1);
t164 = pkin(7) * t184 + qJD(2) * pkin(8);
t179 = t175 * t156 - t173 * t164;
t178 = qJD(4) - t179;
t145 = -t160 * pkin(9) + (pkin(3) + pkin(4)) * t167 + t178;
t186 = t173 * t156 + t175 * t164;
t150 = -t167 * qJ(4) + t186;
t159 = -t175 * qJD(2) + t173 * t185;
t147 = t159 * pkin(9) + t150;
t172 = sin(qJ(5));
t187 = t172 * t145 + t189 * t147;
t163 = -qJD(2) * pkin(2) + pkin(7) * t185;
t183 = qJD(1) * qJD(2);
t182 = t174 * t183;
t181 = t176 * t183;
t180 = t189 * t145 - t172 * t147;
t151 = t159 * pkin(3) - t160 * qJ(4) + t163;
t148 = -t159 * pkin(4) - t151;
t166 = qJD(5) + t167;
t153 = t172 * t159 + t189 * t160;
t152 = -t189 * t159 + t172 * t160;
t149 = t167 * pkin(3) + t178;
t142 = t152 * pkin(5) + qJD(6) + t148;
t141 = -t152 * qJ(6) + t187;
t140 = t166 * pkin(5) - t153 * qJ(6) + t180;
t1 = [t190, 0, 0, t174 ^ 2 * t190, t174 * t188, t182, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t188 - pkin(7) * t182, -t177 * pkin(1) * t174 - pkin(7) * t181, t160 ^ 2 / 0.2e1, -t160 * t159, -t160 * t167, t159 * t167, t167 ^ 2 / 0.2e1, t163 * t159 - t179 * t167, t163 * t160 + t186 * t167, t149 * t167 + t151 * t159, t149 * t160 - t150 * t159, -t150 * t167 - t151 * t160, t150 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t153 ^ 2 / 0.2e1, -t153 * t152, t166 * t153, -t166 * t152, t166 ^ 2 / 0.2e1, t148 * t152 + t180 * t166, t148 * t153 - t187 * t166, -t140 * t153 - t141 * t152, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1;];
T_reg  = t1;
