% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR15_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:50:11
% EndTime: 2019-03-10 00:50:12
% DurationCPUTime: 0.24s
% Computational Cost: add. (1014->63), mult. (2702->128), div. (0->0), fcn. (2233->12), ass. (0->55)
t281 = cos(pkin(6)) * qJD(1);
t256 = qJD(2) + t281;
t258 = sin(pkin(7));
t260 = cos(pkin(7));
t269 = cos(qJ(2));
t259 = sin(pkin(6));
t282 = qJD(1) * t259;
t276 = t269 * t282;
t291 = t256 * t258 + t260 * t276;
t265 = sin(qJ(2));
t280 = pkin(1) * t281;
t283 = pkin(9) * t276 + t265 * t280;
t239 = t291 * pkin(10) + t283;
t255 = t269 * t280;
t277 = t265 * t282;
t241 = t256 * pkin(2) + t255 + (-pkin(10) * t260 - pkin(9)) * t277;
t247 = (-pkin(10) * t258 * t265 - pkin(2) * t269 - pkin(1)) * t282;
t264 = sin(qJ(3));
t268 = cos(qJ(3));
t290 = -t264 * t239 + (t241 * t260 + t247 * t258) * t268;
t289 = pkin(4) + pkin(12);
t270 = qJD(1) ^ 2;
t287 = t259 ^ 2 * t270;
t286 = t258 * t264;
t285 = t260 * t264;
t231 = -t258 * t241 + t260 * t247;
t242 = t264 * t277 - t291 * t268;
t243 = t256 * t286 + (t265 * t268 + t269 * t285) * t282;
t224 = t242 * pkin(3) - t243 * pkin(11) + t231;
t248 = -t260 * t256 + t258 * t276 - qJD(3);
t278 = t268 * t239 + t241 * t285 + t247 * t286;
t228 = -t248 * pkin(11) + t278;
t263 = sin(qJ(4));
t267 = cos(qJ(4));
t284 = t263 * t224 + t267 * t228;
t279 = t269 * t287;
t275 = t267 * t224 - t263 * t228;
t240 = qJD(4) + t242;
t221 = -t240 * qJ(5) - t284;
t273 = qJD(5) - t275;
t235 = t267 * t243 - t263 * t248;
t227 = t248 * pkin(3) - t290;
t271 = -t235 * qJ(5) + t227;
t266 = cos(qJ(6));
t262 = sin(qJ(6));
t234 = t263 * t243 + t267 * t248;
t233 = qJD(6) + t235;
t230 = t262 * t234 + t266 * t240;
t229 = -t266 * t234 + t262 * t240;
t222 = t234 * pkin(4) + t271;
t220 = -t240 * pkin(4) + t273;
t219 = t289 * t234 + t271;
t218 = -t234 * pkin(5) - t221;
t217 = t235 * pkin(5) - t289 * t240 + t273;
t1 = [t270 / 0.2e1, 0, 0, t265 ^ 2 * t287 / 0.2e1, t265 * t279, t256 * t277, t256 * t276, t256 ^ 2 / 0.2e1, pkin(1) * t279 + (-pkin(9) * t277 + t255) * t256, -pkin(1) * t265 * t287 - t283 * t256, t243 ^ 2 / 0.2e1, -t243 * t242, -t243 * t248, t242 * t248, t248 ^ 2 / 0.2e1, t231 * t242 - t290 * t248, t231 * t243 + t278 * t248, t235 ^ 2 / 0.2e1, -t235 * t234, t235 * t240, -t234 * t240, t240 ^ 2 / 0.2e1, t227 * t234 + t275 * t240, t227 * t235 - t284 * t240, t220 * t235 + t221 * t234, t220 * t240 - t222 * t234, -t221 * t240 - t222 * t235, t222 ^ 2 / 0.2e1 + t221 ^ 2 / 0.2e1 + t220 ^ 2 / 0.2e1, t230 ^ 2 / 0.2e1, -t230 * t229, t230 * t233, -t229 * t233, t233 ^ 2 / 0.2e1 (t266 * t217 - t262 * t219) * t233 + t218 * t229 -(t262 * t217 + t266 * t219) * t233 + t218 * t230;];
T_reg  = t1;
