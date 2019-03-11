% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:36
% EndTime: 2019-03-10 02:53:37
% DurationCPUTime: 0.24s
% Computational Cost: add. (1020->60), mult. (2710->124), div. (0->0), fcn. (2254->12), ass. (0->54)
t278 = cos(pkin(6)) * qJD(1);
t255 = qJD(2) + t278;
t257 = sin(pkin(7));
t259 = cos(pkin(7));
t267 = cos(qJ(2));
t258 = sin(pkin(6));
t279 = qJD(1) * t258;
t273 = t267 * t279;
t289 = t255 * t257 + t259 * t273;
t264 = sin(qJ(2));
t277 = pkin(1) * t278;
t280 = pkin(9) * t273 + t264 * t277;
t238 = t289 * pkin(10) + t280;
t254 = t267 * t277;
t274 = t264 * t279;
t240 = t255 * pkin(2) + t254 + (-pkin(10) * t259 - pkin(9)) * t274;
t246 = (-pkin(10) * t257 * t264 - pkin(2) * t267 - pkin(1)) * t279;
t263 = sin(qJ(3));
t266 = cos(qJ(3));
t288 = -t263 * t238 + (t240 * t259 + t246 * t257) * t266;
t287 = cos(qJ(5));
t268 = qJD(1) ^ 2;
t285 = t258 ^ 2 * t268;
t284 = t257 * t263;
t283 = t259 * t263;
t241 = t263 * t274 - t289 * t266;
t239 = qJD(4) + t241;
t230 = -t257 * t240 + t259 * t246;
t242 = t255 * t284 + (t264 * t266 + t267 * t283) * t279;
t223 = t241 * pkin(3) - t242 * pkin(11) + t230;
t247 = -t259 * t255 + t257 * t273 - qJD(3);
t275 = t266 * t238 + t240 * t283 + t246 * t284;
t227 = -t247 * pkin(11) + t275;
t262 = sin(qJ(4));
t265 = cos(qJ(4));
t281 = t262 * t223 + t265 * t227;
t218 = t239 * pkin(12) + t281;
t226 = t247 * pkin(3) - t288;
t232 = t262 * t242 + t265 * t247;
t233 = t265 * t242 - t262 * t247;
t221 = t232 * pkin(4) - t233 * pkin(12) + t226;
t261 = sin(qJ(5));
t282 = t287 * t218 + t261 * t221;
t276 = t267 * t285;
t272 = -t261 * t218 + t287 * t221;
t271 = t265 * t223 - t262 * t227;
t217 = -t239 * pkin(4) - t271;
t231 = qJD(5) + t232;
t229 = t287 * t233 + t261 * t239;
t228 = t261 * t233 - t287 * t239;
t215 = t228 * pkin(5) + qJD(6) + t217;
t214 = -t228 * qJ(6) + t282;
t213 = t231 * pkin(5) - t229 * qJ(6) + t272;
t1 = [t268 / 0.2e1, 0, 0, t264 ^ 2 * t285 / 0.2e1, t264 * t276, t255 * t274, t255 * t273, t255 ^ 2 / 0.2e1, pkin(1) * t276 + (-pkin(9) * t274 + t254) * t255, -pkin(1) * t264 * t285 - t280 * t255, t242 ^ 2 / 0.2e1, -t242 * t241, -t242 * t247, t241 * t247, t247 ^ 2 / 0.2e1, t230 * t241 - t288 * t247, t230 * t242 + t275 * t247, t233 ^ 2 / 0.2e1, -t233 * t232, t233 * t239, -t232 * t239, t239 ^ 2 / 0.2e1, t226 * t232 + t271 * t239, t226 * t233 - t281 * t239, t229 ^ 2 / 0.2e1, -t229 * t228, t229 * t231, -t228 * t231, t231 ^ 2 / 0.2e1, t217 * t228 + t272 * t231, t217 * t229 - t282 * t231, -t213 * t229 - t214 * t228, t214 ^ 2 / 0.2e1 + t213 ^ 2 / 0.2e1 + t215 ^ 2 / 0.2e1;];
T_reg  = t1;
