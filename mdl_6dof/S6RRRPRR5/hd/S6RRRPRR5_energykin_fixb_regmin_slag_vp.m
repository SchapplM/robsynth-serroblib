% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:09
% EndTime: 2019-03-09 18:23:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (368->51), mult. (787->105), div. (0->0), fcn. (548->8), ass. (0->46)
t206 = pkin(3) + pkin(9);
t205 = -pkin(8) - pkin(7);
t190 = qJD(1) ^ 2;
t204 = t190 / 0.2e1;
t203 = cos(qJ(5));
t189 = cos(qJ(2));
t202 = t189 * t190;
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t188 = cos(qJ(3));
t174 = (t185 * t189 + t186 * t188) * qJD(1);
t182 = qJD(2) + qJD(3);
t199 = qJD(1) * t186;
t177 = qJD(2) * pkin(2) + t205 * t199;
t198 = qJD(1) * t189;
t178 = t205 * t198;
t193 = t188 * t177 + t185 * t178;
t192 = qJD(4) - t193;
t158 = t174 * pkin(4) - t206 * t182 + t192;
t173 = t185 * t199 - t188 * t198;
t179 = (-pkin(2) * t189 - pkin(1)) * qJD(1);
t191 = -t174 * qJ(4) + t179;
t159 = t206 * t173 + t191;
t184 = sin(qJ(5));
t201 = t184 * t158 + t203 * t159;
t200 = t185 * t177 - t188 * t178;
t197 = qJD(1) * qJD(2);
t165 = -t182 * qJ(4) - t200;
t196 = t186 * t197;
t195 = t189 * t197;
t194 = t203 * t158 - t184 * t159;
t162 = -t173 * pkin(4) - t165;
t172 = qJD(5) + t174;
t187 = cos(qJ(6));
t183 = sin(qJ(6));
t170 = qJD(6) + t172;
t168 = t184 * t173 + t203 * t182;
t167 = -t203 * t173 + t184 * t182;
t164 = -t182 * pkin(3) + t192;
t163 = t173 * pkin(3) + t191;
t161 = -t183 * t167 + t187 * t168;
t160 = t187 * t167 + t183 * t168;
t154 = t167 * pkin(5) + t162;
t153 = -t167 * pkin(10) + t201;
t152 = t172 * pkin(5) - t168 * pkin(10) + t194;
t1 = [t204, 0, 0, t186 ^ 2 * t204, t186 * t202, t196, t195, qJD(2) ^ 2 / 0.2e1, pkin(1) * t202 - pkin(7) * t196, -t190 * pkin(1) * t186 - pkin(7) * t195, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t182, -t173 * t182, t182 ^ 2 / 0.2e1, t179 * t173 + t193 * t182, t179 * t174 - t200 * t182, t164 * t174 + t165 * t173, -t163 * t173 + t164 * t182, -t163 * t174 - t165 * t182, t163 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t172, -t167 * t172, t172 ^ 2 / 0.2e1, t162 * t167 + t194 * t172, t162 * t168 - t201 * t172, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t170, -t160 * t170, t170 ^ 2 / 0.2e1 (t187 * t152 - t183 * t153) * t170 + t154 * t160 -(t183 * t152 + t187 * t153) * t170 + t154 * t161;];
T_reg  = t1;
