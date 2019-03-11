% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:14
% EndTime: 2019-03-08 19:16:14
% DurationCPUTime: 0.10s
% Computational Cost: add. (191->42), mult. (461->89), div. (0->0), fcn. (345->12), ass. (0->39)
t188 = cos(qJ(2));
t195 = qJD(1) * sin(pkin(6));
t170 = qJD(2) * pkin(2) + t188 * t195;
t179 = sin(pkin(11));
t182 = cos(pkin(11));
t185 = sin(qJ(2));
t192 = t185 * t195;
t163 = t179 * t170 + t182 * t192;
t161 = qJD(2) * qJ(4) + t163;
t176 = cos(pkin(6)) * qJD(1) + qJD(3);
t181 = cos(pkin(12));
t172 = t181 * t176;
t178 = sin(pkin(12));
t154 = t172 + (-pkin(8) * qJD(2) - t161) * t178;
t157 = t181 * t161 + t178 * t176;
t193 = qJD(2) * t181;
t155 = pkin(8) * t193 + t157;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t196 = t184 * t154 + t187 * t155;
t194 = qJD(2) * t178;
t191 = qJD(2) * t195;
t162 = t182 * t170 - t179 * t192;
t167 = t184 * t194 - t187 * t193;
t190 = qJD(4) - t162;
t189 = t187 * t154 - t184 * t155;
t158 = (-pkin(4) * t181 - pkin(3)) * qJD(2) + t190;
t186 = cos(qJ(6));
t183 = sin(qJ(6));
t168 = (t178 * t187 + t181 * t184) * qJD(2);
t166 = qJD(6) + t167;
t165 = t183 * qJD(5) + t186 * t168;
t164 = -t186 * qJD(5) + t183 * t168;
t160 = -qJD(2) * pkin(3) + t190;
t156 = -t178 * t161 + t172;
t151 = t167 * pkin(5) - t168 * pkin(9) + t158;
t150 = qJD(5) * pkin(9) + t196;
t149 = -qJD(5) * pkin(5) - t189;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t188 * t191, -t185 * t191, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1, -t160 * t193, t160 * t194 (-t156 * t178 + t157 * t181) * qJD(2), t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * qJD(5), -t167 * qJD(5), qJD(5) ^ 2 / 0.2e1, t189 * qJD(5) + t158 * t167, -t196 * qJD(5) + t158 * t168, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t166, -t164 * t166, t166 ^ 2 / 0.2e1 (-t183 * t150 + t186 * t151) * t166 + t149 * t164 -(t186 * t150 + t183 * t151) * t166 + t149 * t165;];
T_reg  = t1;
