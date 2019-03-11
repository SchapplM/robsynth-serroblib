% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:37
% EndTime: 2019-03-08 22:27:37
% DurationCPUTime: 0.21s
% Computational Cost: add. (497->52), mult. (1225->113), div. (0->0), fcn. (1006->14), ass. (0->50)
t248 = cos(qJ(2));
t260 = qJD(1) * sin(pkin(6));
t228 = qJD(2) * pkin(2) + t248 * t260;
t237 = sin(pkin(7));
t239 = cos(pkin(7));
t259 = qJD(1) * cos(pkin(6));
t250 = t228 * t239 + t237 * t259;
t264 = cos(pkin(13));
t249 = qJD(2) ^ 2;
t262 = t237 ^ 2 * t249;
t234 = qJD(2) * t239 + qJD(3);
t244 = sin(qJ(2));
t258 = qJD(2) * t237;
t226 = pkin(9) * t258 + t244 * t260;
t243 = sin(qJ(3));
t247 = cos(qJ(3));
t257 = t247 * t226 + t243 * t250;
t212 = qJ(4) * t234 + t257;
t233 = t239 * t259;
t218 = t233 + (-t228 + (-pkin(3) * t247 - qJ(4) * t243) * qJD(2)) * t237;
t236 = sin(pkin(13));
t205 = -t212 * t236 + t218 * t264;
t255 = t243 * t258;
t222 = t234 * t236 + t255 * t264;
t253 = t247 * t258;
t202 = -pkin(4) * t253 - pkin(10) * t222 + t205;
t206 = t212 * t264 + t218 * t236;
t221 = -t234 * t264 + t236 * t255;
t204 = -pkin(10) * t221 + t206;
t242 = sin(qJ(5));
t246 = cos(qJ(5));
t261 = t202 * t242 + t204 * t246;
t254 = (-t237 * t228 + t233) * t258;
t252 = qJD(2) * t260;
t214 = t221 * t246 + t222 * t242;
t251 = t202 * t246 - t204 * t242;
t224 = t243 * t226;
t211 = -t234 * pkin(3) - t247 * t250 + qJD(4) + t224;
t207 = t221 * pkin(4) + t211;
t245 = cos(qJ(6));
t241 = sin(qJ(6));
t231 = -qJD(5) + t253;
t215 = -t221 * t242 + t222 * t246;
t213 = qJD(6) + t214;
t209 = t215 * t245 - t231 * t241;
t208 = t215 * t241 + t231 * t245;
t200 = t214 * pkin(5) - t215 * pkin(11) + t207;
t199 = -pkin(11) * t231 + t261;
t198 = pkin(5) * t231 - t251;
t1 = [qJD(1) ^ 2 / 0.2e1, t249 / 0.2e1, t248 * t252, -t244 * t252, t243 ^ 2 * t262 / 0.2e1, t243 * t247 * t262, t234 * t255, t234 * t253, t234 ^ 2 / 0.2e1, -t224 * t234 + (t234 * t250 - t254) * t247, -t234 * t257 + t243 * t254, -t205 * t253 + t211 * t221, t206 * t253 + t211 * t222, -t205 * t222 - t206 * t221, t206 ^ 2 / 0.2e1 + t205 ^ 2 / 0.2e1 + t211 ^ 2 / 0.2e1, t215 ^ 2 / 0.2e1, -t215 * t214, -t215 * t231, t214 * t231, t231 ^ 2 / 0.2e1, t207 * t214 - t231 * t251, t207 * t215 + t231 * t261, t209 ^ 2 / 0.2e1, -t209 * t208, t209 * t213, -t208 * t213, t213 ^ 2 / 0.2e1 (-t199 * t241 + t200 * t245) * t213 + t198 * t208 -(t199 * t245 + t200 * t241) * t213 + t198 * t209;];
T_reg  = t1;
