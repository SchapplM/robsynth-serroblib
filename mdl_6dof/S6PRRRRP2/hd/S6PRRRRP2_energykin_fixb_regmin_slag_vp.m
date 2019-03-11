% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:19
% EndTime: 2019-03-09 00:04:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (328->45), mult. (719->97), div. (0->0), fcn. (527->10), ass. (0->41)
t189 = qJD(2) ^ 2;
t202 = t189 / 0.2e1;
t178 = qJD(3) + qJD(4);
t184 = sin(qJ(2));
t198 = qJD(1) * sin(pkin(6));
t172 = qJD(2) * pkin(8) + t184 * t198;
t187 = cos(qJ(3));
t197 = qJD(1) * cos(pkin(6));
t175 = t187 * t197;
t183 = sin(qJ(3));
t163 = qJD(3) * pkin(3) + t175 + (-pkin(9) * qJD(2) - t172) * t183;
t195 = qJD(2) * t187;
t199 = t187 * t172 + t183 * t197;
t164 = pkin(9) * t195 + t199;
t182 = sin(qJ(4));
t186 = cos(qJ(4));
t200 = t182 * t163 + t186 * t164;
t157 = t178 * pkin(10) + t200;
t188 = cos(qJ(2));
t193 = t188 * t198;
t167 = -t193 + (-pkin(3) * t187 - pkin(2)) * qJD(2);
t196 = qJD(2) * t183;
t169 = t182 * t196 - t186 * t195;
t170 = (t182 * t187 + t183 * t186) * qJD(2);
t159 = t169 * pkin(4) - t170 * pkin(10) + t167;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t201 = t185 * t157 + t181 * t159;
t194 = qJD(2) * qJD(3);
t192 = qJD(2) * t198;
t191 = t186 * t163 - t182 * t164;
t190 = -t181 * t157 + t185 * t159;
t156 = -t178 * pkin(4) - t191;
t173 = -qJD(2) * pkin(2) - t193;
t168 = qJD(5) + t169;
t166 = t185 * t170 + t181 * t178;
t165 = t181 * t170 - t185 * t178;
t154 = t165 * pkin(5) - t166 * qJ(6) + t156;
t153 = t168 * qJ(6) + t201;
t152 = -t168 * pkin(5) + qJD(6) - t190;
t1 = [qJD(1) ^ 2 / 0.2e1, t202, t188 * t192, -t184 * t192, t183 ^ 2 * t202, t183 * t189 * t187, t183 * t194, t187 * t194, qJD(3) ^ 2 / 0.2e1 (-t183 * t172 + t175) * qJD(3) - t173 * t195, -t199 * qJD(3) + t173 * t196, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t178, -t169 * t178, t178 ^ 2 / 0.2e1, t167 * t169 + t191 * t178, t167 * t170 - t200 * t178, t166 ^ 2 / 0.2e1, -t166 * t165, t166 * t168, -t165 * t168, t168 ^ 2 / 0.2e1, t156 * t165 + t190 * t168, t156 * t166 - t201 * t168, -t152 * t168 + t154 * t165, t152 * t166 - t153 * t165, t153 * t168 - t154 * t166, t153 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1;];
T_reg  = t1;
