% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:05
% EndTime: 2019-03-09 18:04:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (596->52), mult. (1470->109), div. (0->0), fcn. (1144->10), ass. (0->49)
t206 = -pkin(8) - pkin(7);
t193 = qJD(1) ^ 2;
t205 = t193 / 0.2e1;
t204 = cos(qJ(3));
t192 = cos(qJ(2));
t203 = t192 * t193;
t188 = sin(qJ(3));
t189 = sin(qJ(2));
t175 = (t188 * t192 + t204 * t189) * qJD(1);
t183 = qJD(2) + qJD(3);
t200 = qJD(1) * t189;
t177 = qJD(2) * pkin(2) + t206 * t200;
t199 = qJD(1) * t192;
t178 = t206 * t199;
t195 = t204 * t177 + t188 * t178;
t165 = t183 * pkin(3) - t175 * qJ(4) + t195;
t174 = t188 * t200 - t204 * t199;
t201 = t188 * t177 - t204 * t178;
t167 = -t174 * qJ(4) + t201;
t184 = sin(pkin(11));
t185 = cos(pkin(11));
t155 = t185 * t165 - t184 * t167;
t170 = -t184 * t174 + t185 * t175;
t152 = t183 * pkin(4) - t170 * pkin(9) + t155;
t156 = t184 * t165 + t185 * t167;
t169 = -t185 * t174 - t184 * t175;
t154 = t169 * pkin(9) + t156;
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t202 = t187 * t152 + t191 * t154;
t198 = qJD(1) * qJD(2);
t197 = t189 * t198;
t196 = t192 * t198;
t160 = -t191 * t169 + t187 * t170;
t179 = (-pkin(2) * t192 - pkin(1)) * qJD(1);
t194 = t191 * t152 - t187 * t154;
t171 = t174 * pkin(3) + qJD(4) + t179;
t162 = -t169 * pkin(4) + t171;
t190 = cos(qJ(6));
t186 = sin(qJ(6));
t182 = qJD(5) + t183;
t161 = t187 * t169 + t191 * t170;
t159 = qJD(6) + t160;
t158 = t190 * t161 + t186 * t182;
t157 = t186 * t161 - t190 * t182;
t150 = t160 * pkin(5) - t161 * pkin(10) + t162;
t149 = t182 * pkin(10) + t202;
t148 = -t182 * pkin(5) - t194;
t1 = [t205, 0, 0, t189 ^ 2 * t205, t189 * t203, t197, t196, qJD(2) ^ 2 / 0.2e1, pkin(1) * t203 - pkin(7) * t197, -t193 * pkin(1) * t189 - pkin(7) * t196, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t183, -t174 * t183, t183 ^ 2 / 0.2e1, t179 * t174 + t195 * t183, t179 * t175 - t201 * t183, -t155 * t170 + t156 * t169, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t182, -t160 * t182, t182 ^ 2 / 0.2e1, t162 * t160 + t194 * t182, t162 * t161 - t202 * t182, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t159, -t157 * t159, t159 ^ 2 / 0.2e1 (-t186 * t149 + t190 * t150) * t159 + t148 * t157 -(t190 * t149 + t186 * t150) * t159 + t148 * t158;];
T_reg  = t1;
