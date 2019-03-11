% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:46
% EndTime: 2019-03-08 20:48:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (163->41), mult. (346->87), div. (0->0), fcn. (231->10), ass. (0->41)
t177 = qJD(2) ^ 2;
t194 = t177 / 0.2e1;
t193 = qJD(1) ^ 2 / 0.2e1;
t192 = cos(qJ(5));
t176 = cos(qJ(2));
t189 = qJD(1) * sin(pkin(6));
t179 = -t176 * t189 + qJD(3);
t157 = (-pkin(2) - pkin(8)) * qJD(2) + t179;
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t169 = cos(pkin(6));
t188 = qJD(1) * t169;
t190 = t172 * t157 + t175 * t188;
t151 = qJD(4) * pkin(9) + t190;
t173 = sin(qJ(2));
t183 = t173 * t189;
t154 = t183 + (pkin(4) * t172 - pkin(9) * t175 + qJ(3)) * qJD(2);
t171 = sin(qJ(5));
t191 = t192 * t151 + t171 * t154;
t187 = qJD(2) * t175;
t161 = qJD(2) * qJ(3) + t183;
t186 = t161 * qJD(2);
t185 = t172 * qJD(2);
t184 = qJD(2) * qJD(4);
t182 = qJD(2) * t189;
t181 = -t171 * t151 + t192 * t154;
t180 = t175 * t157 - t172 * t188;
t166 = qJD(5) + t185;
t150 = -qJD(4) * pkin(4) - t180;
t174 = cos(qJ(6));
t170 = sin(qJ(6));
t163 = qJD(6) + t166;
t160 = t171 * qJD(4) + t187 * t192;
t159 = -qJD(4) * t192 + t171 * t187;
t158 = -qJD(2) * pkin(2) + t179;
t149 = -t170 * t159 + t174 * t160;
t148 = t174 * t159 + t170 * t160;
t146 = t159 * pkin(5) + t150;
t145 = -t159 * pkin(10) + t191;
t144 = t166 * pkin(5) - t160 * pkin(10) + t181;
t1 = [t193, t194, t176 * t182, -t173 * t182, t158 * qJD(2), t186, t169 ^ 2 * t193 + t161 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t175 ^ 2 * t194, -t175 * t177 * t172, t175 * t184, -t172 * t184, qJD(4) ^ 2 / 0.2e1, qJD(4) * t180 + t161 * t185, -qJD(4) * t190 + t175 * t186, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t166, -t159 * t166, t166 ^ 2 / 0.2e1, t150 * t159 + t166 * t181, t150 * t160 - t166 * t191, t149 ^ 2 / 0.2e1, -t149 * t148, t149 * t163, -t148 * t163, t163 ^ 2 / 0.2e1 (t174 * t144 - t170 * t145) * t163 + t146 * t148 -(t170 * t144 + t174 * t145) * t163 + t146 * t149;];
T_reg  = t1;
