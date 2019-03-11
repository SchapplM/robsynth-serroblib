% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:21
% EndTime: 2019-03-09 00:45:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (316->48), mult. (694->104), div. (0->0), fcn. (533->12), ass. (0->46)
t199 = qJD(2) ^ 2;
t213 = t199 / 0.2e1;
t212 = cos(qJ(5));
t187 = qJD(3) + qJD(4);
t194 = sin(qJ(2));
t208 = qJD(1) * sin(pkin(6));
t181 = qJD(2) * pkin(8) + t194 * t208;
t197 = cos(qJ(3));
t207 = qJD(1) * cos(pkin(6));
t184 = t197 * t207;
t193 = sin(qJ(3));
t170 = qJD(3) * pkin(3) + t184 + (-pkin(9) * qJD(2) - t181) * t193;
t205 = qJD(2) * t197;
t209 = t197 * t181 + t193 * t207;
t172 = pkin(9) * t205 + t209;
t192 = sin(qJ(4));
t196 = cos(qJ(4));
t210 = t192 * t170 + t196 * t172;
t161 = t187 * pkin(10) + t210;
t198 = cos(qJ(2));
t203 = t198 * t208;
t176 = -t203 + (-pkin(3) * t197 - pkin(2)) * qJD(2);
t206 = qJD(2) * t193;
t178 = t192 * t206 - t196 * t205;
t179 = (t192 * t197 + t193 * t196) * qJD(2);
t166 = t178 * pkin(4) - t179 * pkin(10) + t176;
t191 = sin(qJ(5));
t211 = t212 * t161 + t191 * t166;
t204 = qJD(2) * qJD(3);
t202 = qJD(2) * t208;
t201 = -t191 * t161 + t212 * t166;
t200 = t196 * t170 - t192 * t172;
t160 = -t187 * pkin(4) - t200;
t177 = qJD(5) + t178;
t195 = cos(qJ(6));
t190 = sin(qJ(6));
t182 = -qJD(2) * pkin(2) - t203;
t175 = qJD(6) + t177;
t174 = t212 * t179 + t191 * t187;
t173 = t191 * t179 - t212 * t187;
t163 = -t190 * t173 + t195 * t174;
t162 = t195 * t173 + t190 * t174;
t158 = t173 * pkin(5) + t160;
t157 = -t173 * pkin(11) + t211;
t156 = t177 * pkin(5) - t174 * pkin(11) + t201;
t1 = [qJD(1) ^ 2 / 0.2e1, t213, t198 * t202, -t194 * t202, t193 ^ 2 * t213, t193 * t199 * t197, t193 * t204, t197 * t204, qJD(3) ^ 2 / 0.2e1 (-t193 * t181 + t184) * qJD(3) - t182 * t205, -t209 * qJD(3) + t182 * t206, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t187, -t178 * t187, t187 ^ 2 / 0.2e1, t176 * t178 + t200 * t187, t176 * t179 - t210 * t187, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t177, -t173 * t177, t177 ^ 2 / 0.2e1, t160 * t173 + t201 * t177, t160 * t174 - t211 * t177, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t175, -t162 * t175, t175 ^ 2 / 0.2e1 (t195 * t156 - t190 * t157) * t175 + t158 * t162 -(t190 * t156 + t195 * t157) * t175 + t158 * t163;];
T_reg  = t1;
