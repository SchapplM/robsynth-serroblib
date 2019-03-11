% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:02
% EndTime: 2019-03-09 03:03:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (234->42), mult. (542->92), div. (0->0), fcn. (341->8), ass. (0->37)
t166 = qJD(1) ^ 2;
t176 = t166 / 0.2e1;
t175 = cos(qJ(5));
t160 = sin(pkin(9));
t153 = (pkin(1) * t160 + pkin(7)) * qJD(1);
t165 = cos(qJ(3));
t158 = t165 * qJD(2);
t164 = sin(qJ(3));
t144 = qJD(3) * pkin(3) + t158 + (-qJ(4) * qJD(1) - t153) * t164;
t171 = qJD(1) * t165;
t173 = t164 * qJD(2) + t165 * t153;
t147 = qJ(4) * t171 + t173;
t159 = sin(pkin(10));
t161 = cos(pkin(10));
t137 = t159 * t144 + t161 * t147;
t135 = qJD(3) * pkin(8) + t137;
t162 = cos(pkin(9));
t169 = -pkin(1) * t162 - pkin(2);
t149 = qJD(4) + (-pkin(3) * t165 + t169) * qJD(1);
t172 = qJD(1) * t164;
t150 = -t159 * t172 + t161 * t171;
t151 = (t159 * t165 + t161 * t164) * qJD(1);
t140 = -t150 * pkin(4) - t151 * pkin(8) + t149;
t163 = sin(qJ(5));
t174 = t175 * t135 + t163 * t140;
t170 = qJD(1) * qJD(3);
t168 = -t163 * t135 + t175 * t140;
t136 = t161 * t144 - t159 * t147;
t134 = -qJD(3) * pkin(4) - t136;
t154 = t169 * qJD(1);
t148 = qJD(5) - t150;
t146 = t163 * qJD(3) + t175 * t151;
t145 = -t175 * qJD(3) + t163 * t151;
t132 = t145 * pkin(5) + qJD(6) + t134;
t131 = -t145 * qJ(6) + t174;
t130 = t148 * pkin(5) - t146 * qJ(6) + t168;
t1 = [t176, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t160 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t166, t164 ^ 2 * t176, t164 * t166 * t165, t164 * t170, t165 * t170, qJD(3) ^ 2 / 0.2e1, -t154 * t171 + (-t164 * t153 + t158) * qJD(3), -t173 * qJD(3) + t154 * t172, -t136 * t151 + t137 * t150, t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t148, -t145 * t148, t148 ^ 2 / 0.2e1, t134 * t145 + t168 * t148, t134 * t146 - t174 * t148, -t130 * t146 - t131 * t145, t131 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1;];
T_reg  = t1;
