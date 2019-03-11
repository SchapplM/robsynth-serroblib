% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:40
% EndTime: 2019-03-08 19:11:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (343->45), mult. (876->105), div. (0->0), fcn. (762->16), ass. (0->52)
t238 = cos(pkin(6)) * qJD(1) + qJD(2);
t243 = sin(pkin(7));
t247 = cos(pkin(7));
t245 = cos(pkin(14));
t244 = sin(pkin(6));
t270 = qJD(1) * t244;
t265 = t245 * t270;
t276 = t238 * t243 + t247 * t265;
t251 = sin(qJ(3));
t255 = cos(qJ(3));
t241 = sin(pkin(14));
t266 = t241 * t270;
t258 = -t251 * t266 + t276 * t255;
t219 = qJD(3) * pkin(3) + t258;
t225 = t247 * t238 - t243 * t265;
t242 = sin(pkin(8));
t246 = cos(pkin(8));
t259 = t219 * t246 + t225 * t242;
t256 = qJD(3) ^ 2;
t272 = t242 ^ 2 * t256;
t237 = t246 * qJD(3) + qJD(4);
t267 = t276 * t251 + t255 * t266;
t269 = qJD(3) * t242;
t218 = pkin(10) * t269 + t267;
t250 = sin(qJ(4));
t254 = cos(qJ(4));
t268 = t254 * t218 + t259 * t250;
t211 = t237 * pkin(11) + t268;
t223 = t246 * t225;
t213 = t223 + (-t219 + (-pkin(4) * t254 - pkin(11) * t250) * qJD(3)) * t242;
t249 = sin(qJ(5));
t253 = cos(qJ(5));
t271 = t253 * t211 + t249 * t213;
t264 = t250 * t269;
t263 = (-t242 * t219 + t223) * t269;
t262 = t254 * t269;
t260 = -t249 * t211 + t253 * t213;
t226 = -t253 * t237 + t249 * t264;
t216 = t250 * t218;
t210 = -t237 * pkin(4) - t259 * t254 + t216;
t257 = qJD(1) ^ 2;
t252 = cos(qJ(6));
t248 = sin(qJ(6));
t235 = -qJD(5) + t262;
t227 = t249 * t237 + t253 * t264;
t224 = qJD(6) + t226;
t221 = t252 * t227 - t248 * t235;
t220 = t248 * t227 + t252 * t235;
t208 = t226 * pkin(5) - t227 * pkin(12) + t210;
t207 = -t235 * pkin(12) + t271;
t206 = t235 * pkin(5) - t260;
t1 = [t257 / 0.2e1, t238 ^ 2 / 0.2e1 + (t241 ^ 2 / 0.2e1 + t245 ^ 2 / 0.2e1) * t257 * t244 ^ 2, t256 / 0.2e1, t258 * qJD(3), -t267 * qJD(3), t250 ^ 2 * t272 / 0.2e1, t254 * t250 * t272, t237 * t264, t237 * t262, t237 ^ 2 / 0.2e1, -t216 * t237 + (t259 * t237 - t263) * t254, -t268 * t237 + t250 * t263, t227 ^ 2 / 0.2e1, -t227 * t226, -t227 * t235, t226 * t235, t235 ^ 2 / 0.2e1, t210 * t226 - t260 * t235, t210 * t227 + t271 * t235, t221 ^ 2 / 0.2e1, -t221 * t220, t221 * t224, -t220 * t224, t224 ^ 2 / 0.2e1 (-t248 * t207 + t252 * t208) * t224 + t206 * t220 -(t252 * t207 + t248 * t208) * t224 + t206 * t221;];
T_reg  = t1;
