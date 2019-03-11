% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:38
% EndTime: 2019-03-09 08:15:38
% DurationCPUTime: 0.09s
% Computational Cost: add. (250->51), mult. (506->103), div. (0->0), fcn. (249->6), ass. (0->37)
t178 = -pkin(2) - pkin(3);
t168 = qJD(1) ^ 2;
t177 = t168 / 0.2e1;
t167 = cos(qJ(2));
t176 = t167 * t168;
t165 = sin(qJ(2));
t174 = t165 * qJD(1);
t175 = qJD(1) * t167;
t148 = -qJD(1) * pkin(1) - pkin(2) * t175 - qJ(3) * t174;
t143 = pkin(3) * t175 + qJD(4) - t148;
t140 = (pkin(4) * t165 + qJ(5) * t167) * qJD(1) + t143;
t173 = pkin(7) * t174 + qJD(3);
t169 = -qJ(4) * t174 + t173;
t142 = (-qJ(5) + t178) * qJD(2) + t169;
t160 = sin(pkin(9));
t161 = cos(pkin(9));
t134 = t160 * t140 + t161 * t142;
t152 = pkin(7) * t175 + qJD(2) * qJ(3);
t172 = qJD(1) * qJD(2);
t171 = t165 * t172;
t170 = t167 * t172;
t133 = t161 * t140 - t160 * t142;
t147 = qJ(4) * t175 - t152;
t144 = qJD(2) * pkin(4) + qJD(5) - t147;
t166 = cos(qJ(6));
t164 = sin(qJ(6));
t153 = qJD(6) + t174;
t151 = -qJD(2) * pkin(2) + t173;
t150 = t160 * qJD(2) + t161 * t175;
t149 = -t161 * qJD(2) + t160 * t175;
t145 = t178 * qJD(2) + t169;
t139 = t164 * t149 - t166 * t150;
t138 = -t166 * t149 - t164 * t150;
t137 = -t149 * pkin(5) + t144;
t132 = t149 * pkin(8) + t134;
t131 = pkin(5) * t174 + t150 * pkin(8) + t133;
t1 = [t177, 0, 0, t165 ^ 2 * t177, t165 * t176, t171, t170, qJD(2) ^ 2 / 0.2e1, pkin(1) * t176 - pkin(7) * t171, -t168 * pkin(1) * t165 - pkin(7) * t170, -t151 * qJD(2) - t148 * t175 (t151 * t165 + t152 * t167) * qJD(1), t152 * qJD(2) - t148 * t174, t152 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1, -t147 * qJD(2) + t143 * t174, t145 * qJD(2) - t143 * t175 (-t145 * t165 + t147 * t167) * qJD(1), t145 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t133 * t174 - t144 * t149, -t134 * t174 - t144 * t150, t133 * t150 + t134 * t149, t134 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t153, -t138 * t153, t153 ^ 2 / 0.2e1 (t166 * t131 - t164 * t132) * t153 + t137 * t138 -(t164 * t131 + t166 * t132) * t153 + t137 * t139;];
T_reg  = t1;
