% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR3
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
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:42
% EndTime: 2019-03-08 22:07:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (412->50), mult. (1134->109), div. (0->0), fcn. (936->14), ass. (0->51)
t253 = cos(qJ(2));
t264 = qJD(1) * sin(pkin(6));
t232 = qJD(2) * pkin(2) + t253 * t264;
t241 = sin(pkin(7));
t244 = cos(pkin(7));
t263 = qJD(1) * cos(pkin(6));
t269 = t232 * t244 + t241 * t263;
t254 = qJD(2) ^ 2;
t267 = t241 ^ 2 * t254;
t249 = sin(qJ(2));
t262 = qJD(2) * t241;
t231 = pkin(9) * t262 + t249 * t264;
t238 = t244 * qJD(2) + qJD(3);
t248 = sin(qJ(3));
t252 = cos(qJ(3));
t265 = t269 * t252;
t215 = t238 * pkin(3) + (-qJ(4) * t262 - t231) * t248 + t265;
t261 = qJD(2) * t252;
t257 = t241 * t261;
t260 = t252 * t231 + t269 * t248;
t218 = qJ(4) * t257 + t260;
t240 = sin(pkin(13));
t243 = cos(pkin(13));
t209 = t240 * t215 + t243 * t218;
t207 = t238 * pkin(10) + t209;
t237 = t244 * t263;
t222 = qJD(4) + t237 + (-pkin(3) * t261 - t232) * t241;
t258 = t248 * t262;
t226 = -t240 * t258 + t243 * t257;
t227 = (t240 * t252 + t243 * t248) * t262;
t211 = -t226 * pkin(4) - t227 * pkin(10) + t222;
t247 = sin(qJ(5));
t251 = cos(qJ(5));
t266 = t251 * t207 + t247 * t211;
t256 = qJD(2) * t264;
t208 = t243 * t215 - t240 * t218;
t220 = t247 * t227 - t251 * t238;
t255 = -t247 * t207 + t251 * t211;
t206 = -t238 * pkin(4) - t208;
t250 = cos(qJ(6));
t246 = sin(qJ(6));
t225 = qJD(5) - t226;
t224 = -t241 * t232 + t237;
t221 = t251 * t227 + t247 * t238;
t219 = qJD(6) + t220;
t214 = t250 * t221 + t246 * t225;
t213 = t246 * t221 - t250 * t225;
t204 = t220 * pkin(5) - t221 * pkin(11) + t206;
t203 = t225 * pkin(11) + t266;
t202 = -t225 * pkin(5) - t255;
t1 = [qJD(1) ^ 2 / 0.2e1, t254 / 0.2e1, t253 * t256, -t249 * t256, t248 ^ 2 * t267 / 0.2e1, t248 * t252 * t267, t238 * t258, t238 * t257, t238 ^ 2 / 0.2e1 (-t248 * t231 + t265) * t238 - t224 * t257, t224 * t258 - t260 * t238, -t208 * t227 + t209 * t226, t209 ^ 2 / 0.2e1 + t208 ^ 2 / 0.2e1 + t222 ^ 2 / 0.2e1, t221 ^ 2 / 0.2e1, -t221 * t220, t221 * t225, -t220 * t225, t225 ^ 2 / 0.2e1, t206 * t220 + t255 * t225, t206 * t221 - t266 * t225, t214 ^ 2 / 0.2e1, -t214 * t213, t214 * t219, -t213 * t219, t219 ^ 2 / 0.2e1 (-t246 * t203 + t250 * t204) * t219 + t202 * t213 -(t250 * t203 + t246 * t204) * t219 + t202 * t214;];
T_reg  = t1;
