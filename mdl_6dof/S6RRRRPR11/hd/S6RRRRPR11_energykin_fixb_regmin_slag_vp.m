% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:33
% EndTime: 2019-03-09 23:23:33
% DurationCPUTime: 0.19s
% Computational Cost: add. (841->56), mult. (1910->116), div. (0->0), fcn. (1541->12), ass. (0->52)
t249 = cos(qJ(4));
t225 = sin(pkin(6));
t235 = qJD(1) ^ 2;
t248 = t225 ^ 2 * t235;
t243 = cos(pkin(6)) * qJD(1);
t222 = qJD(2) + t243;
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t231 = sin(qJ(2));
t244 = qJD(1) * t225;
t240 = t231 * t244;
t214 = t230 * t222 + t233 * t240;
t234 = cos(qJ(2));
t239 = t234 * t244;
t217 = -qJD(3) + t239;
t229 = sin(qJ(4));
t204 = t249 * t214 - t229 * t217;
t213 = -t233 * t222 + t230 * t240;
t212 = qJD(4) + t213;
t242 = pkin(1) * t243;
t236 = -pkin(8) * t240 + t234 * t242;
t208 = -t222 * pkin(2) - t236;
t199 = t213 * pkin(3) - t214 * pkin(10) + t208;
t245 = pkin(8) * t239 + t231 * t242;
t209 = t222 * pkin(9) + t245;
t211 = (-pkin(2) * t234 - pkin(9) * t231 - pkin(1)) * t244;
t246 = t233 * t209 + t230 * t211;
t202 = -t217 * pkin(10) + t246;
t238 = t249 * t199 - t229 * t202;
t187 = t212 * pkin(4) - t204 * qJ(5) + t238;
t203 = t229 * t214 + t249 * t217;
t247 = t229 * t199 + t249 * t202;
t190 = -t203 * qJ(5) + t247;
t224 = sin(pkin(12));
t226 = cos(pkin(12));
t184 = t224 * t187 + t226 * t190;
t241 = t234 * t248;
t183 = t226 * t187 - t224 * t190;
t237 = -t230 * t209 + t233 * t211;
t201 = t217 * pkin(3) - t237;
t193 = t203 * pkin(4) + qJD(5) + t201;
t232 = cos(qJ(6));
t228 = sin(qJ(6));
t210 = qJD(6) + t212;
t196 = -t224 * t203 + t226 * t204;
t195 = -t226 * t203 - t224 * t204;
t192 = t228 * t195 + t232 * t196;
t191 = -t232 * t195 + t228 * t196;
t188 = -t195 * pkin(5) + t193;
t182 = t195 * pkin(11) + t184;
t181 = t212 * pkin(5) - t196 * pkin(11) + t183;
t1 = [t235 / 0.2e1, 0, 0, t231 ^ 2 * t248 / 0.2e1, t231 * t241, t222 * t240, t222 * t239, t222 ^ 2 / 0.2e1, pkin(1) * t241 + t236 * t222, -pkin(1) * t231 * t248 - t245 * t222, t214 ^ 2 / 0.2e1, -t214 * t213, -t214 * t217, t213 * t217, t217 ^ 2 / 0.2e1, t208 * t213 - t237 * t217, t208 * t214 + t246 * t217, t204 ^ 2 / 0.2e1, -t204 * t203, t204 * t212, -t203 * t212, t212 ^ 2 / 0.2e1, t201 * t203 + t238 * t212, t201 * t204 - t247 * t212, -t183 * t196 + t184 * t195, t184 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t210, -t191 * t210, t210 ^ 2 / 0.2e1 (t232 * t181 - t228 * t182) * t210 + t188 * t191 -(t228 * t181 + t232 * t182) * t210 + t188 * t192;];
T_reg  = t1;
