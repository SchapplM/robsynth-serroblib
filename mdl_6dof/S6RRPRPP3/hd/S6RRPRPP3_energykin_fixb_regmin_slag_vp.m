% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:08
% EndTime: 2019-03-09 09:57:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (458->50), mult. (1009->100), div. (0->0), fcn. (666->6), ass. (0->39)
t174 = qJD(1) ^ 2;
t186 = t174 / 0.2e1;
t185 = pkin(4) + qJ(6);
t173 = cos(qJ(2));
t184 = t173 * t174;
t171 = sin(qJ(2));
t158 = (-pkin(2) * t173 - qJ(3) * t171 - pkin(1)) * qJD(1);
t181 = t173 * qJD(1);
t163 = pkin(7) * t181 + qJD(2) * qJ(3);
t168 = sin(pkin(9));
t169 = cos(pkin(9));
t152 = t169 * t158 - t168 * t163;
t182 = qJD(1) * t171;
t160 = t168 * qJD(2) + t169 * t182;
t146 = -pkin(3) * t181 - t160 * pkin(8) + t152;
t153 = t168 * t158 + t169 * t163;
t159 = -t169 * qJD(2) + t168 * t182;
t149 = -t159 * pkin(8) + t153;
t170 = sin(qJ(4));
t172 = cos(qJ(4));
t183 = t170 * t146 + t172 * t149;
t180 = qJD(1) * qJD(2);
t179 = t171 * t180;
t178 = t173 * t180;
t177 = t172 * t146 - t170 * t149;
t162 = -qJD(2) * pkin(2) + pkin(7) * t182 + qJD(3);
t164 = -qJD(4) + t181;
t143 = t164 * qJ(5) - t183;
t176 = qJD(5) - t177;
t154 = t159 * pkin(3) + t162;
t151 = -t170 * t159 + t172 * t160;
t175 = -t151 * qJ(5) + t154;
t150 = t172 * t159 + t170 * t160;
t144 = t150 * pkin(4) + t175;
t142 = t164 * pkin(4) + t176;
t141 = t185 * t150 + t175;
t140 = -t150 * pkin(5) + qJD(6) - t143;
t139 = t151 * pkin(5) + t185 * t164 + t176;
t1 = [t186, 0, 0, t171 ^ 2 * t186, t171 * t184, t179, t178, qJD(2) ^ 2 / 0.2e1, pkin(1) * t184 - pkin(7) * t179, -t174 * pkin(1) * t171 - pkin(7) * t178, -t152 * t181 + t162 * t159, t153 * t181 + t162 * t160, -t152 * t160 - t153 * t159, t153 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t151 ^ 2 / 0.2e1, -t151 * t150, -t151 * t164, t150 * t164, t164 ^ 2 / 0.2e1, t154 * t150 - t177 * t164, t154 * t151 + t183 * t164, t142 * t151 + t143 * t150, -t142 * t164 - t144 * t150, t143 * t164 - t144 * t151, t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t139 * t151 - t140 * t150, -t140 * t164 - t141 * t151, t139 * t164 + t141 * t150, t141 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1;];
T_reg  = t1;
