% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:58
% EndTime: 2019-03-09 09:20:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (303->53), mult. (722->112), div. (0->0), fcn. (486->8), ass. (0->43)
t202 = -pkin(2) - pkin(3);
t181 = sin(pkin(6));
t189 = qJD(1) ^ 2;
t201 = t181 ^ 2 * t189;
t182 = cos(pkin(6));
t197 = t182 * qJD(1);
t178 = qJD(2) + t197;
t185 = sin(qJ(2));
t198 = qJD(1) * t181;
t193 = t185 * t198;
t173 = pkin(8) * t193;
t188 = cos(qJ(2));
t190 = qJD(3) + t173 + (-pkin(1) * t182 * t188 - qJ(4) * t181 * t185) * qJD(1);
t154 = (-pkin(9) + t202) * t178 + t190;
t194 = t188 * t198;
t164 = -pkin(1) * t198 - pkin(2) * t194 - qJ(3) * t193;
t161 = pkin(3) * t194 + qJD(4) - t164;
t155 = (pkin(4) * t185 + pkin(9) * t188) * t198 + t161;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t200 = t154 * t187 + t155 * t184;
t196 = pkin(1) * t197;
t199 = pkin(8) * t194 + t185 * t196;
t195 = t188 * t201;
t163 = qJ(3) * t178 + t199;
t192 = -t154 * t184 + t155 * t187;
t191 = t188 * t196 - t173;
t167 = -t178 * t187 + t184 * t194;
t160 = qJ(4) * t194 - t163;
t157 = pkin(4) * t178 - t160;
t186 = cos(qJ(6));
t183 = sin(qJ(6));
t170 = qJD(5) + t193;
t168 = -t178 * t184 - t187 * t194;
t165 = -qJD(6) + t167;
t162 = -t178 * pkin(2) + qJD(3) - t191;
t159 = t168 * t186 + t170 * t183;
t158 = t168 * t183 - t170 * t186;
t156 = t178 * t202 + t190;
t151 = -pkin(5) * t167 - pkin(10) * t168 + t157;
t150 = pkin(10) * t170 + t200;
t149 = -pkin(5) * t170 - t192;
t1 = [t189 / 0.2e1, 0, 0, t185 ^ 2 * t201 / 0.2e1, t185 * t195, t178 * t193, t178 * t194, t178 ^ 2 / 0.2e1, pkin(1) * t195 + t178 * t191, -pkin(1) * t185 * t201 - t178 * t199, -t162 * t178 - t164 * t194 (t162 * t185 + t163 * t188) * t198, t163 * t178 - t164 * t193, t163 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, -t160 * t178 + t161 * t193, t156 * t178 - t161 * t194 (-t156 * t185 + t160 * t188) * t198, t156 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, t168 * t167, t168 * t170, t167 * t170, t170 ^ 2 / 0.2e1, -t157 * t167 + t170 * t192, t157 * t168 - t170 * t200, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t165, t158 * t165, t165 ^ 2 / 0.2e1 -(-t150 * t183 + t151 * t186) * t165 + t149 * t158 (t150 * t186 + t151 * t183) * t165 + t149 * t159;];
T_reg  = t1;
