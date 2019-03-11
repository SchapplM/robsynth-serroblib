% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:16
% EndTime: 2019-03-09 10:26:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (696->56), mult. (1930->114), div. (0->0), fcn. (1547->12), ass. (0->51)
t262 = cos(qJ(4));
t242 = sin(pkin(6));
t251 = qJD(1) ^ 2;
t261 = t242 ^ 2 * t251;
t241 = sin(pkin(11));
t244 = cos(pkin(11));
t248 = sin(qJ(2));
t250 = cos(qJ(2));
t258 = qJD(1) * t242;
t231 = (t241 * t250 + t244 * t248) * t258;
t257 = cos(pkin(6)) * qJD(1);
t238 = qJD(2) + t257;
t247 = sin(qJ(4));
t223 = t262 * t231 + t247 * t238;
t253 = t250 * t258;
t254 = t248 * t258;
t230 = -t241 * t254 + t244 * t253;
t229 = qJD(4) - t230;
t256 = pkin(1) * t257;
t237 = t250 * t256;
t224 = t238 * pkin(2) + t237 + (-pkin(8) - qJ(3)) * t254;
t259 = pkin(8) * t253 + t248 * t256;
t227 = qJ(3) * t253 + t259;
t216 = t241 * t224 + t244 * t227;
t214 = t238 * pkin(9) + t216;
t232 = qJD(3) + (-pkin(2) * t250 - pkin(1)) * t258;
t219 = -t230 * pkin(3) - t231 * pkin(9) + t232;
t252 = -t247 * t214 + t262 * t219;
t203 = t229 * pkin(4) - t223 * qJ(5) + t252;
t222 = t247 * t231 - t262 * t238;
t260 = t262 * t214 + t247 * t219;
t205 = -t222 * qJ(5) + t260;
t240 = sin(pkin(12));
t243 = cos(pkin(12));
t200 = t240 * t203 + t243 * t205;
t255 = t250 * t261;
t211 = -t243 * t222 - t240 * t223;
t215 = t244 * t224 - t241 * t227;
t199 = t243 * t203 - t240 * t205;
t213 = -t238 * pkin(3) - t215;
t206 = t222 * pkin(4) + qJD(5) + t213;
t249 = cos(qJ(6));
t246 = sin(qJ(6));
t212 = -t240 * t222 + t243 * t223;
t209 = qJD(6) - t211;
t208 = t249 * t212 + t246 * t229;
t207 = t246 * t212 - t249 * t229;
t201 = -t211 * pkin(5) - t212 * pkin(10) + t206;
t198 = t229 * pkin(10) + t200;
t197 = -t229 * pkin(5) - t199;
t1 = [t251 / 0.2e1, 0, 0, t248 ^ 2 * t261 / 0.2e1, t248 * t255, t238 * t254, t238 * t253, t238 ^ 2 / 0.2e1, pkin(1) * t255 + (-pkin(8) * t254 + t237) * t238, -pkin(1) * t248 * t261 - t259 * t238, -t215 * t231 + t216 * t230, t216 ^ 2 / 0.2e1 + t215 ^ 2 / 0.2e1 + t232 ^ 2 / 0.2e1, t223 ^ 2 / 0.2e1, -t223 * t222, t223 * t229, -t222 * t229, t229 ^ 2 / 0.2e1, t213 * t222 + t252 * t229, t213 * t223 - t260 * t229, -t199 * t212 + t200 * t211, t200 ^ 2 / 0.2e1 + t199 ^ 2 / 0.2e1 + t206 ^ 2 / 0.2e1, t208 ^ 2 / 0.2e1, -t208 * t207, t208 * t209, -t207 * t209, t209 ^ 2 / 0.2e1 (-t246 * t198 + t249 * t201) * t209 + t197 * t207 -(t249 * t198 + t246 * t201) * t209 + t197 * t208;];
T_reg  = t1;
