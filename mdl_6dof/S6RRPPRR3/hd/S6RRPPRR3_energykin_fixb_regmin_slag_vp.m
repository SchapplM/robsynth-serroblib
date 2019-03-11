% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:36
% EndTime: 2019-03-09 08:58:37
% DurationCPUTime: 0.19s
% Computational Cost: add. (722->58), mult. (2019->118), div. (0->0), fcn. (1621->12), ass. (0->51)
t258 = cos(pkin(12));
t238 = sin(pkin(6));
t247 = qJD(1) ^ 2;
t257 = t238 ^ 2 * t247;
t246 = cos(qJ(2));
t253 = cos(pkin(6)) * qJD(1);
t252 = pkin(1) * t253;
t233 = t246 * t252;
t234 = qJD(2) + t253;
t243 = sin(qJ(2));
t254 = qJD(1) * t238;
t250 = t243 * t254;
t220 = t234 * pkin(2) + t233 + (-pkin(8) - qJ(3)) * t250;
t249 = t246 * t254;
t255 = pkin(8) * t249 + t243 * t252;
t223 = qJ(3) * t249 + t255;
t237 = sin(pkin(11));
t239 = cos(pkin(11));
t212 = t220 * t237 + t223 * t239;
t210 = qJ(4) * t234 + t212;
t226 = t237 * t250 - t239 * t249;
t227 = (t237 * t246 + t239 * t243) * t254;
t228 = qJD(3) + (-pkin(2) * t246 - pkin(1)) * t254;
t215 = t226 * pkin(3) - t227 * qJ(4) + t228;
t236 = sin(pkin(12));
t200 = -t210 * t236 + t215 * t258;
t219 = t227 * t258 + t234 * t236;
t197 = pkin(4) * t226 - pkin(9) * t219 + t200;
t201 = t210 * t258 + t215 * t236;
t218 = t227 * t236 - t234 * t258;
t199 = -pkin(9) * t218 + t201;
t242 = sin(qJ(5));
t245 = cos(qJ(5));
t256 = t197 * t242 + t199 * t245;
t251 = t246 * t257;
t207 = t218 * t245 + t219 * t242;
t211 = t220 * t239 - t223 * t237;
t248 = t197 * t245 - t199 * t242;
t209 = -pkin(3) * t234 + qJD(4) - t211;
t202 = pkin(4) * t218 + t209;
t244 = cos(qJ(6));
t241 = sin(qJ(6));
t225 = qJD(5) + t226;
t208 = -t218 * t242 + t219 * t245;
t205 = qJD(6) + t207;
t204 = t208 * t244 + t225 * t241;
t203 = t208 * t241 - t225 * t244;
t195 = pkin(5) * t207 - pkin(10) * t208 + t202;
t194 = pkin(10) * t225 + t256;
t193 = -pkin(5) * t225 - t248;
t1 = [t247 / 0.2e1, 0, 0, t243 ^ 2 * t257 / 0.2e1, t243 * t251, t234 * t250, t234 * t249, t234 ^ 2 / 0.2e1, pkin(1) * t251 + (-pkin(8) * t250 + t233) * t234, -pkin(1) * t243 * t257 - t234 * t255, -t211 * t227 - t212 * t226, t212 ^ 2 / 0.2e1 + t211 ^ 2 / 0.2e1 + t228 ^ 2 / 0.2e1, t200 * t226 + t209 * t218, -t201 * t226 + t209 * t219, -t200 * t219 - t201 * t218, t201 ^ 2 / 0.2e1 + t200 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1, t208 ^ 2 / 0.2e1, -t208 * t207, t208 * t225, -t207 * t225, t225 ^ 2 / 0.2e1, t202 * t207 + t225 * t248, t202 * t208 - t225 * t256, t204 ^ 2 / 0.2e1, -t204 * t203, t204 * t205, -t203 * t205, t205 ^ 2 / 0.2e1 (-t194 * t241 + t195 * t244) * t205 + t193 * t203 -(t194 * t244 + t195 * t241) * t205 + t193 * t204;];
T_reg  = t1;
