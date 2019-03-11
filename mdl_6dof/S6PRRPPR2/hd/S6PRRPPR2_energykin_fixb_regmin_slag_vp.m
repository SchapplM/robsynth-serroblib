% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:30
% EndTime: 2019-03-08 21:06:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (230->45), mult. (556->94), div. (0->0), fcn. (393->10), ass. (0->41)
t201 = pkin(4) + pkin(9);
t189 = qJD(2) ^ 2;
t200 = t189 / 0.2e1;
t185 = sin(qJ(2));
t198 = qJD(1) * sin(pkin(6));
t174 = qJD(2) * pkin(8) + t185 * t198;
t187 = cos(qJ(3));
t197 = qJD(1) * cos(pkin(6));
t178 = t187 * t197;
t184 = sin(qJ(3));
t164 = qJD(3) * pkin(3) + t178 + (-qJ(4) * qJD(2) - t174) * t184;
t195 = qJD(2) * t187;
t199 = t187 * t174 + t184 * t197;
t165 = qJ(4) * t195 + t199;
t179 = sin(pkin(11));
t181 = cos(pkin(11));
t159 = t179 * t164 + t181 * t165;
t196 = qJD(2) * t184;
t194 = qJD(2) * qJD(3);
t188 = cos(qJ(2));
t193 = t188 * t198;
t192 = qJD(2) * t198;
t158 = t181 * t164 - t179 * t165;
t156 = -qJD(3) * qJ(5) - t159;
t191 = qJD(5) - t158;
t172 = (t179 * t187 + t181 * t184) * qJD(2);
t169 = -t193 + qJD(4) + (-pkin(3) * t187 - pkin(2)) * qJD(2);
t190 = -t172 * qJ(5) + t169;
t186 = cos(qJ(6));
t183 = sin(qJ(6));
t175 = -qJD(2) * pkin(2) - t193;
t171 = t179 * t196 - t181 * t195;
t170 = qJD(6) + t172;
t167 = t186 * qJD(3) + t183 * t171;
t166 = t183 * qJD(3) - t186 * t171;
t160 = t171 * pkin(4) + t190;
t157 = t201 * t171 + t190;
t155 = -qJD(3) * pkin(4) + t191;
t154 = -t171 * pkin(5) - t156;
t153 = t172 * pkin(5) - t201 * qJD(3) + t191;
t1 = [qJD(1) ^ 2 / 0.2e1, t200, t188 * t192, -t185 * t192, t184 ^ 2 * t200, t184 * t189 * t187, t184 * t194, t187 * t194, qJD(3) ^ 2 / 0.2e1 (-t184 * t174 + t178) * qJD(3) - t175 * t195, -t199 * qJD(3) + t175 * t196, -t158 * t172 - t159 * t171, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t155 * t172 + t156 * t171, t155 * qJD(3) - t160 * t171, -t156 * qJD(3) - t160 * t172, t160 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t170, -t166 * t170, t170 ^ 2 / 0.2e1 (t186 * t153 - t183 * t157) * t170 + t154 * t166 -(t183 * t153 + t186 * t157) * t170 + t154 * t167;];
T_reg  = t1;
