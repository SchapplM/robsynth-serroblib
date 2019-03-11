% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:29
% EndTime: 2019-03-09 16:09:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (413->54), mult. (970->105), div. (0->0), fcn. (702->8), ass. (0->43)
t213 = -pkin(4) - pkin(10);
t190 = sin(pkin(6));
t198 = qJD(1) ^ 2;
t212 = t190 ^ 2 * t198;
t207 = cos(pkin(6)) * qJD(1);
t188 = qJD(2) + t207;
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t208 = qJD(1) * t190;
t203 = t197 * t208;
t206 = pkin(1) * t207;
t209 = pkin(8) * t203 + t194 * t206;
t171 = t188 * pkin(9) + t209;
t172 = (-pkin(2) * t197 - pkin(9) * t194 - pkin(1)) * t208;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t211 = t196 * t171 + t193 * t172;
t204 = t194 * t208;
t210 = -pkin(8) * t204 + t197 * t206;
t205 = t197 * t212;
t181 = -qJD(3) + t203;
t164 = -t181 * qJ(4) + t211;
t170 = -t188 * pkin(2) - t210;
t202 = -t193 * t171 + t196 * t172;
t201 = qJD(4) - t202;
t176 = -t196 * t188 + t193 * t204;
t177 = t193 * t188 + t196 * t204;
t162 = t176 * pkin(3) - t177 * qJ(4) + t170;
t161 = -t176 * qJ(5) - t164;
t200 = qJD(5) - t162;
t199 = -t177 * qJ(5) + t201;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t175 = qJD(6) + t177;
t166 = t195 * t176 + t192 * t181;
t165 = t192 * t176 - t195 * t181;
t163 = t181 * pkin(3) + t201;
t160 = -t176 * pkin(4) + t200;
t159 = -t181 * pkin(5) - t161;
t158 = (pkin(3) + pkin(4)) * t181 + t199;
t157 = (pkin(3) - t213) * t181 + t199;
t156 = t177 * pkin(5) + t213 * t176 + t200;
t1 = [t198 / 0.2e1, 0, 0, t194 ^ 2 * t212 / 0.2e1, t194 * t205, t188 * t204, t188 * t203, t188 ^ 2 / 0.2e1, pkin(1) * t205 + t210 * t188, -pkin(1) * t194 * t212 - t209 * t188, t177 ^ 2 / 0.2e1, -t177 * t176, -t177 * t181, t176 * t181, t181 ^ 2 / 0.2e1, t170 * t176 - t202 * t181, t170 * t177 + t211 * t181, t162 * t176 + t163 * t181, t163 * t177 - t164 * t176, -t162 * t177 - t164 * t181, t164 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t160 * t177 + t161 * t181, -t158 * t181 + t160 * t176, -t158 * t177 - t161 * t176, t158 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1, t166 ^ 2 / 0.2e1, -t166 * t165, t166 * t175, -t165 * t175, t175 ^ 2 / 0.2e1 (t195 * t156 - t192 * t157) * t175 + t159 * t165 -(t192 * t156 + t195 * t157) * t175 + t159 * t166;];
T_reg  = t1;
