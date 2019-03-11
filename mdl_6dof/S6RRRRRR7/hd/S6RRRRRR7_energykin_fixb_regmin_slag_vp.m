% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:05
% EndTime: 2019-03-10 04:41:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (820->57), mult. (1841->119), div. (0->0), fcn. (1515->12), ass. (0->53)
t251 = cos(qJ(5));
t224 = sin(pkin(6));
t235 = qJD(1) ^ 2;
t250 = t224 ^ 2 * t235;
t244 = cos(pkin(6)) * qJD(1);
t222 = qJD(2) + t244;
t229 = sin(qJ(3));
t233 = cos(qJ(3));
t230 = sin(qJ(2));
t245 = qJD(1) * t224;
t241 = t230 * t245;
t214 = t229 * t222 + t233 * t241;
t234 = cos(qJ(2));
t240 = t234 * t245;
t217 = -qJD(3) + t240;
t228 = sin(qJ(4));
t232 = cos(qJ(4));
t203 = t232 * t214 - t228 * t217;
t213 = -t233 * t222 + t229 * t241;
t212 = qJD(4) + t213;
t243 = pkin(1) * t244;
t236 = -pkin(8) * t241 + t234 * t243;
t207 = -t222 * pkin(2) - t236;
t197 = t213 * pkin(3) - t214 * pkin(10) + t207;
t246 = pkin(8) * t240 + t230 * t243;
t208 = t222 * pkin(9) + t246;
t211 = (-pkin(2) * t234 - pkin(9) * t230 - pkin(1)) * t245;
t247 = t233 * t208 + t229 * t211;
t200 = -t217 * pkin(10) + t247;
t238 = t232 * t197 - t228 * t200;
t185 = t212 * pkin(4) - t203 * pkin(11) + t238;
t202 = t228 * t214 + t232 * t217;
t248 = t228 * t197 + t232 * t200;
t188 = -t202 * pkin(11) + t248;
t227 = sin(qJ(5));
t249 = t227 * t185 + t251 * t188;
t242 = t234 * t250;
t239 = t251 * t185 - t227 * t188;
t237 = -t229 * t208 + t233 * t211;
t199 = t217 * pkin(3) - t237;
t210 = qJD(5) + t212;
t191 = t202 * pkin(4) + t199;
t231 = cos(qJ(6));
t226 = sin(qJ(6));
t209 = qJD(6) + t210;
t194 = -t227 * t202 + t251 * t203;
t193 = t251 * t202 + t227 * t203;
t190 = -t226 * t193 + t231 * t194;
t189 = t231 * t193 + t226 * t194;
t186 = t193 * pkin(5) + t191;
t182 = -t193 * pkin(12) + t249;
t181 = t210 * pkin(5) - t194 * pkin(12) + t239;
t1 = [t235 / 0.2e1, 0, 0, t230 ^ 2 * t250 / 0.2e1, t230 * t242, t222 * t241, t222 * t240, t222 ^ 2 / 0.2e1, pkin(1) * t242 + t236 * t222, -pkin(1) * t230 * t250 - t246 * t222, t214 ^ 2 / 0.2e1, -t214 * t213, -t214 * t217, t213 * t217, t217 ^ 2 / 0.2e1, t207 * t213 - t237 * t217, t207 * t214 + t247 * t217, t203 ^ 2 / 0.2e1, -t203 * t202, t203 * t212, -t202 * t212, t212 ^ 2 / 0.2e1, t199 * t202 + t238 * t212, t199 * t203 - t248 * t212, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t210, -t193 * t210, t210 ^ 2 / 0.2e1, t191 * t193 + t210 * t239, t191 * t194 - t249 * t210, t190 ^ 2 / 0.2e1, -t190 * t189, t190 * t209, -t189 * t209, t209 ^ 2 / 0.2e1 (t231 * t181 - t226 * t182) * t209 + t186 * t189 -(t226 * t181 + t231 * t182) * t209 + t186 * t190;];
T_reg  = t1;
