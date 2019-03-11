% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:35
% EndTime: 2019-03-09 12:48:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (275->47), mult. (577->97), div. (0->0), fcn. (341->6), ass. (0->40)
t190 = -pkin(2) - pkin(8);
t175 = qJD(1) ^ 2;
t189 = t175 / 0.2e1;
t188 = cos(qJ(5));
t174 = cos(qJ(2));
t187 = t174 * t175;
t171 = sin(qJ(4));
t173 = cos(qJ(4));
t184 = qJD(1) * t174;
t161 = t173 * qJD(2) - t171 * t184;
t172 = sin(qJ(2));
t183 = t172 * qJD(1);
t165 = qJD(4) + t183;
t178 = -qJ(3) * t172 - pkin(1);
t155 = (t190 * t174 + t178) * qJD(1);
t182 = pkin(7) * t183 + qJD(3);
t156 = pkin(3) * t183 + t190 * qJD(2) + t182;
t176 = -t171 * t155 + t173 * t156;
t145 = t165 * pkin(4) - t161 * pkin(9) + t176;
t160 = t171 * qJD(2) + t173 * t184;
t185 = t173 * t155 + t171 * t156;
t148 = -t160 * pkin(9) + t185;
t170 = sin(qJ(5));
t186 = t170 * t145 + t188 * t148;
t163 = -pkin(7) * t184 - qJD(2) * qJ(3);
t181 = qJD(1) * qJD(2);
t158 = pkin(3) * t184 - t163;
t180 = t172 * t181;
t179 = t174 * t181;
t177 = t188 * t145 - t170 * t148;
t151 = t160 * pkin(4) + t158;
t164 = qJD(5) + t165;
t162 = -qJD(2) * pkin(2) + t182;
t159 = (-pkin(2) * t174 + t178) * qJD(1);
t150 = -t170 * t160 + t161 * t188;
t149 = t160 * t188 + t170 * t161;
t146 = t149 * pkin(5) + qJD(6) + t151;
t142 = -t149 * qJ(6) + t186;
t141 = t164 * pkin(5) - t150 * qJ(6) + t177;
t1 = [t189, 0, 0, t172 ^ 2 * t189, t172 * t187, t180, t179, qJD(2) ^ 2 / 0.2e1, pkin(1) * t187 - pkin(7) * t180, -t175 * pkin(1) * t172 - pkin(7) * t179 (t162 * t172 - t163 * t174) * qJD(1), t162 * qJD(2) + t159 * t184, -t163 * qJD(2) - t159 * t183, t159 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t165, -t160 * t165, t165 ^ 2 / 0.2e1, t158 * t160 + t165 * t176, t158 * t161 - t165 * t185, t150 ^ 2 / 0.2e1, -t150 * t149, t164 * t150, -t164 * t149, t164 ^ 2 / 0.2e1, t151 * t149 + t164 * t177, t151 * t150 - t164 * t186, -t141 * t150 - t142 * t149, t142 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1;];
T_reg  = t1;
