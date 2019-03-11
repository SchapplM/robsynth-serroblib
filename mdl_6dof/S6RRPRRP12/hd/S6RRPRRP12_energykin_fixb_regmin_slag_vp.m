% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP12
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:31
% EndTime: 2019-03-09 12:53:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (359->49), mult. (718->101), div. (0->0), fcn. (427->6), ass. (0->40)
t191 = -pkin(2) - pkin(8);
t177 = qJD(1) ^ 2;
t190 = t177 / 0.2e1;
t176 = cos(qJ(2));
t189 = t176 * t177;
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t186 = qJD(1) * t176;
t162 = t175 * qJD(2) - t172 * t186;
t173 = sin(qJ(2));
t185 = t173 * qJD(1);
t166 = qJD(4) + t185;
t180 = -qJ(3) * t173 - pkin(1);
t156 = (t191 * t176 + t180) * qJD(1);
t184 = pkin(7) * t185 + qJD(3);
t157 = pkin(3) * t185 + t191 * qJD(2) + t184;
t179 = -t172 * t156 + t175 * t157;
t147 = t166 * pkin(4) - t162 * pkin(9) + t179;
t161 = t172 * qJD(2) + t175 * t186;
t187 = t175 * t156 + t172 * t157;
t149 = -t161 * pkin(9) + t187;
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t188 = t171 * t147 + t174 * t149;
t164 = -pkin(7) * t186 - qJD(2) * qJ(3);
t183 = qJD(1) * qJD(2);
t159 = pkin(3) * t186 - t164;
t182 = t173 * t183;
t181 = t176 * t183;
t178 = t174 * t147 - t171 * t149;
t152 = t161 * pkin(4) + t159;
t165 = qJD(5) + t166;
t163 = -qJD(2) * pkin(2) + t184;
t160 = (-pkin(2) * t176 + t180) * qJD(1);
t151 = -t171 * t161 + t174 * t162;
t150 = t174 * t161 + t171 * t162;
t145 = t150 * pkin(5) - t151 * qJ(6) + t152;
t144 = t165 * qJ(6) + t188;
t143 = -t165 * pkin(5) + qJD(6) - t178;
t1 = [t190, 0, 0, t173 ^ 2 * t190, t173 * t189, t182, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t189 - pkin(7) * t182, -t177 * pkin(1) * t173 - pkin(7) * t181 (t163 * t173 - t164 * t176) * qJD(1), t163 * qJD(2) + t160 * t186, -t164 * qJD(2) - t160 * t185, t160 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t166, -t161 * t166, t166 ^ 2 / 0.2e1, t159 * t161 + t179 * t166, t159 * t162 - t187 * t166, t151 ^ 2 / 0.2e1, -t151 * t150, t165 * t151, -t165 * t150, t165 ^ 2 / 0.2e1, t152 * t150 + t178 * t165, t152 * t151 - t188 * t165, -t143 * t165 + t145 * t150, t143 * t151 - t144 * t150, t144 * t165 - t145 * t151, t144 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg  = t1;
