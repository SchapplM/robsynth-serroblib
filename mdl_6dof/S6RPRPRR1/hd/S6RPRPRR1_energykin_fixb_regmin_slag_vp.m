% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:19
% EndTime: 2019-03-09 03:35:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (313->47), mult. (740->104), div. (0->0), fcn. (519->10), ass. (0->41)
t178 = qJD(1) ^ 2;
t186 = t178 / 0.2e1;
t169 = sin(pkin(10));
t162 = (pkin(1) * t169 + pkin(7)) * qJD(1);
t177 = cos(qJ(3));
t166 = t177 * qJD(2);
t174 = sin(qJ(3));
t155 = qJD(3) * pkin(3) + t166 + (-qJ(4) * qJD(1) - t162) * t174;
t183 = qJD(1) * t177;
t184 = qJD(2) * t174 + t162 * t177;
t156 = qJ(4) * t183 + t184;
t168 = sin(pkin(11));
t170 = cos(pkin(11));
t144 = t155 * t170 - t156 * t168;
t160 = (t168 * t177 + t170 * t174) * qJD(1);
t142 = qJD(3) * pkin(4) - pkin(8) * t160 + t144;
t145 = t155 * t168 + t156 * t170;
t159 = (-t168 * t174 + t170 * t177) * qJD(1);
t143 = pkin(8) * t159 + t145;
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t185 = t142 * t173 + t143 * t176;
t182 = qJD(1) * qJD(3);
t171 = cos(pkin(10));
t181 = -pkin(1) * t171 - pkin(2);
t149 = -t159 * t176 + t160 * t173;
t180 = t142 * t176 - t143 * t173;
t158 = qJD(4) + (-pkin(3) * t177 + t181) * qJD(1);
t151 = -t159 * pkin(4) + t158;
t175 = cos(qJ(6));
t172 = sin(qJ(6));
t167 = qJD(3) + qJD(5);
t163 = t181 * qJD(1);
t150 = t159 * t173 + t160 * t176;
t148 = qJD(6) + t149;
t147 = t150 * t175 + t167 * t172;
t146 = t150 * t172 - t167 * t175;
t139 = t149 * pkin(5) - t150 * pkin(9) + t151;
t138 = pkin(9) * t167 + t185;
t137 = -pkin(5) * t167 - t180;
t1 = [t186, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t169 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t178, t174 ^ 2 * t186, t174 * t178 * t177, t174 * t182, t177 * t182, qJD(3) ^ 2 / 0.2e1, -t163 * t183 + (-t162 * t174 + t166) * qJD(3), qJD(1) * t163 * t174 - qJD(3) * t184, -t144 * t160 + t145 * t159, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t167, -t149 * t167, t167 ^ 2 / 0.2e1, t149 * t151 + t167 * t180, t150 * t151 - t167 * t185, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * t148, -t146 * t148, t148 ^ 2 / 0.2e1 (-t138 * t172 + t139 * t175) * t148 + t137 * t146 -(t138 * t175 + t139 * t172) * t148 + t137 * t147;];
T_reg  = t1;
