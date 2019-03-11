% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:12
% EndTime: 2019-03-09 06:41:12
% DurationCPUTime: 0.22s
% Computational Cost: add. (843->65), mult. (2726->131), div. (0->0), fcn. (2265->12), ass. (0->55)
t259 = sin(pkin(12));
t285 = cos(pkin(6));
t273 = qJD(1) * t285;
t270 = pkin(1) * t273;
t262 = cos(pkin(12));
t261 = sin(pkin(6));
t278 = qJD(1) * t261;
t275 = t262 * t278;
t251 = qJ(2) * t275 + t259 * t270;
t263 = cos(pkin(7));
t260 = sin(pkin(7));
t274 = t285 * t260;
t282 = t261 * t262;
t239 = (t263 * t282 + t274) * qJD(1) * pkin(9) + t251;
t256 = t262 * t270;
t283 = t259 * t261;
t241 = t256 + (t285 * pkin(2) + (-pkin(9) * t263 - qJ(2)) * t283) * qJD(1);
t266 = sin(qJ(3));
t268 = cos(qJ(3));
t246 = qJD(2) + (-pkin(9) * t259 * t260 - pkin(2) * t262 - pkin(1)) * t278;
t284 = t246 * t260;
t287 = -t266 * t239 + (t241 * t263 + t284) * t268;
t286 = cos(qJ(5));
t281 = t263 * t266;
t276 = t259 * t278;
t242 = t266 * t276 + (-t260 * t273 - t263 * t275) * t268;
t240 = qJD(4) + t242;
t231 = -t260 * t241 + t263 * t246;
t243 = (t266 * t274 + (t259 * t268 + t262 * t281) * t261) * qJD(1);
t224 = t242 * pkin(3) - t243 * pkin(10) + t231;
t248 = t260 * t275 - t263 * t273 - qJD(3);
t277 = t268 * t239 + t241 * t281 + t266 * t284;
t228 = -t248 * pkin(10) + t277;
t265 = sin(qJ(4));
t267 = cos(qJ(4));
t279 = t265 * t224 + t267 * t228;
t219 = t240 * pkin(11) + t279;
t227 = t248 * pkin(3) - t287;
t233 = t265 * t243 + t267 * t248;
t234 = t267 * t243 - t265 * t248;
t222 = t233 * pkin(4) - t234 * pkin(11) + t227;
t264 = sin(qJ(5));
t280 = t286 * t219 + t264 * t222;
t272 = -t264 * t219 + t286 * t222;
t271 = t267 * t224 - t265 * t228;
t218 = -t240 * pkin(4) - t271;
t257 = -pkin(1) * t278 + qJD(2);
t250 = -qJ(2) * t276 + t256;
t232 = qJD(5) + t233;
t230 = t286 * t234 + t264 * t240;
t229 = t264 * t234 - t286 * t240;
t216 = t229 * pkin(5) + qJD(6) + t218;
t215 = -t229 * qJ(6) + t280;
t214 = t232 * pkin(5) - t230 * qJ(6) + t272;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t285 * t250 - t257 * t282) * qJD(1) (-t285 * t251 + t257 * t283) * qJD(1) (-t250 * t259 + t251 * t262) * t278, t251 ^ 2 / 0.2e1 + t250 ^ 2 / 0.2e1 + t257 ^ 2 / 0.2e1, t243 ^ 2 / 0.2e1, -t243 * t242, -t243 * t248, t242 * t248, t248 ^ 2 / 0.2e1, t231 * t242 - t287 * t248, t231 * t243 + t248 * t277, t234 ^ 2 / 0.2e1, -t234 * t233, t234 * t240, -t233 * t240, t240 ^ 2 / 0.2e1, t227 * t233 + t240 * t271, t227 * t234 - t279 * t240, t230 ^ 2 / 0.2e1, -t230 * t229, t230 * t232, -t229 * t232, t232 ^ 2 / 0.2e1, t218 * t229 + t232 * t272, t218 * t230 - t280 * t232, -t214 * t230 - t215 * t229, t215 ^ 2 / 0.2e1 + t214 ^ 2 / 0.2e1 + t216 ^ 2 / 0.2e1;];
T_reg  = t1;
