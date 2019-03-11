% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:55
% EndTime: 2019-03-08 20:20:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (173->38), mult. (346->80), div. (0->0), fcn. (211->8), ass. (0->36)
t166 = qJD(2) ^ 2;
t182 = t166 / 0.2e1;
t181 = qJD(1) ^ 2 / 0.2e1;
t165 = cos(qJ(2));
t178 = qJD(1) * sin(pkin(6));
t169 = -t165 * t178 + qJD(3);
t148 = (-pkin(2) - pkin(8)) * qJD(2) + t169;
t161 = sin(qJ(4));
t164 = cos(qJ(4));
t159 = cos(pkin(6));
t177 = qJD(1) * t159;
t179 = t161 * t148 + t164 * t177;
t144 = qJD(4) * pkin(9) + t179;
t162 = sin(qJ(2));
t172 = t162 * t178;
t146 = t172 + (pkin(4) * t161 - pkin(9) * t164 + qJ(3)) * qJD(2);
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t180 = t163 * t144 + t160 * t146;
t176 = qJD(2) * t164;
t152 = qJD(2) * qJ(3) + t172;
t175 = t152 * qJD(2);
t174 = t161 * qJD(2);
t173 = qJD(2) * qJD(4);
t171 = qJD(2) * t178;
t170 = t164 * t148 - t161 * t177;
t168 = -t160 * t144 + t163 * t146;
t143 = -qJD(4) * pkin(4) - t170;
t156 = qJD(5) + t174;
t151 = t160 * qJD(4) + t163 * t176;
t150 = -t163 * qJD(4) + t160 * t176;
t149 = -qJD(2) * pkin(2) + t169;
t141 = t150 * pkin(5) - t151 * qJ(6) + t143;
t140 = t156 * qJ(6) + t180;
t139 = -t156 * pkin(5) + qJD(6) - t168;
t1 = [t181, t182, t165 * t171, -t162 * t171, t149 * qJD(2), t175, t159 ^ 2 * t181 + t152 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t164 ^ 2 * t182, -t164 * t166 * t161, t164 * t173, -t161 * t173, qJD(4) ^ 2 / 0.2e1, t170 * qJD(4) + t152 * t174, -t179 * qJD(4) + t164 * t175, t151 ^ 2 / 0.2e1, -t151 * t150, t151 * t156, -t150 * t156, t156 ^ 2 / 0.2e1, t143 * t150 + t168 * t156, t143 * t151 - t180 * t156, -t139 * t156 + t141 * t150, t139 * t151 - t140 * t150, t140 * t156 - t141 * t151, t140 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1;];
T_reg  = t1;
