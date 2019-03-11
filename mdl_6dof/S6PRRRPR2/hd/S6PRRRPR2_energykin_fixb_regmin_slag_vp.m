% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:08:57
% EndTime: 2019-03-08 23:08:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (378->49), mult. (833->105), div. (0->0), fcn. (625->12), ass. (0->45)
t200 = qJD(2) ^ 2;
t212 = t200 / 0.2e1;
t211 = cos(pkin(12));
t188 = qJD(3) + qJD(4);
t195 = sin(qJ(2));
t208 = qJD(1) * sin(pkin(6));
t182 = qJD(2) * pkin(8) + t195 * t208;
t198 = cos(qJ(3));
t207 = qJD(1) * cos(pkin(6));
t185 = t198 * t207;
t194 = sin(qJ(3));
t172 = qJD(3) * pkin(3) + t185 + (-pkin(9) * qJD(2) - t182) * t194;
t205 = qJD(2) * t198;
t209 = t198 * t182 + t194 * t207;
t174 = pkin(9) * t205 + t209;
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t210 = t193 * t172 + t197 * t174;
t163 = t188 * qJ(5) + t210;
t199 = cos(qJ(2));
t203 = t199 * t208;
t177 = -t203 + (-pkin(3) * t198 - pkin(2)) * qJD(2);
t206 = qJD(2) * t194;
t179 = t193 * t206 - t197 * t205;
t180 = (t193 * t198 + t194 * t197) * qJD(2);
t168 = t179 * pkin(4) - t180 * qJ(5) + t177;
t189 = sin(pkin(12));
t159 = t211 * t163 + t189 * t168;
t204 = qJD(2) * qJD(3);
t202 = qJD(2) * t208;
t158 = -t189 * t163 + t211 * t168;
t201 = t197 * t172 - t193 * t174;
t162 = -t188 * pkin(4) + qJD(5) - t201;
t196 = cos(qJ(6));
t192 = sin(qJ(6));
t183 = -qJD(2) * pkin(2) - t203;
t178 = qJD(6) + t179;
t176 = t211 * t180 + t189 * t188;
t175 = t189 * t180 - t211 * t188;
t165 = -t192 * t175 + t196 * t176;
t164 = t196 * t175 + t192 * t176;
t160 = t175 * pkin(5) + t162;
t157 = -t175 * pkin(10) + t159;
t156 = t179 * pkin(5) - t176 * pkin(10) + t158;
t1 = [qJD(1) ^ 2 / 0.2e1, t212, t199 * t202, -t195 * t202, t194 ^ 2 * t212, t194 * t200 * t198, t194 * t204, t198 * t204, qJD(3) ^ 2 / 0.2e1 (-t194 * t182 + t185) * qJD(3) - t183 * t205, -t209 * qJD(3) + t183 * t206, t180 ^ 2 / 0.2e1, -t180 * t179, t180 * t188, -t179 * t188, t188 ^ 2 / 0.2e1, t177 * t179 + t201 * t188, t177 * t180 - t210 * t188, t158 * t179 + t162 * t175, -t159 * t179 + t162 * t176, -t158 * t176 - t159 * t175, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t178, -t164 * t178, t178 ^ 2 / 0.2e1 (t196 * t156 - t192 * t157) * t178 + t160 * t164 -(t192 * t156 + t196 * t157) * t178 + t160 * t165;];
T_reg  = t1;
