% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:28
% EndTime: 2019-03-09 15:18:28
% DurationCPUTime: 0.21s
% Computational Cost: add. (863->55), mult. (1899->110), div. (0->0), fcn. (1361->8), ass. (0->46)
t240 = cos(qJ(2));
t252 = t240 * qJD(1);
t231 = qJD(3) - t252;
t238 = sin(qJ(2));
t223 = (-pkin(2) * t240 - pkin(9) * t238 - pkin(1)) * qJD(1);
t230 = pkin(8) * t252 + qJD(2) * pkin(9);
t237 = sin(qJ(3));
t239 = cos(qJ(3));
t246 = t223 * t239 - t230 * t237;
t253 = qJD(1) * t238;
t226 = qJD(2) * t237 + t239 * t253;
t258 = qJ(4) * t226;
t259 = cos(pkin(6));
t213 = pkin(3) * t231 - t258 * t259 + t246;
t225 = -qJD(2) * t239 + t237 * t253;
t229 = -qJD(2) * pkin(2) + pkin(8) * t253;
t235 = sin(pkin(6));
t218 = pkin(3) * t225 - t235 * t258 + t229;
t263 = t213 * t259 + t218 * t235;
t262 = t225 * t259 - t231 * t235;
t254 = t223 * t237 + t230 * t239;
t212 = -qJ(4) * t262 + t254;
t234 = sin(pkin(10));
t236 = cos(pkin(10));
t206 = -t212 * t234 + t236 * t263;
t241 = qJD(1) ^ 2;
t261 = t241 / 0.2e1;
t260 = pkin(4) + qJ(6);
t255 = t240 * t241;
t251 = qJD(1) * qJD(2);
t207 = t236 * t212 + t234 * t263;
t250 = t238 * t251;
t249 = t240 * t251;
t208 = -t235 * t213 + t218 * t259 + qJD(4);
t219 = -t225 * t235 - t231 * t259;
t205 = qJ(5) * t219 - t207;
t215 = t236 * t226 - t234 * t262;
t243 = -t215 * qJ(5) + t208;
t242 = qJD(5) - t206;
t214 = t234 * t226 + t236 * t262;
t204 = t219 * pkin(4) + t242;
t203 = pkin(4) * t214 + t243;
t202 = -pkin(5) * t214 + qJD(6) - t205;
t201 = t214 * t260 + t243;
t200 = t215 * pkin(5) + t219 * t260 + t242;
t1 = [t261, 0, 0, t238 ^ 2 * t261, t238 * t255, t250, t249, qJD(2) ^ 2 / 0.2e1, pkin(1) * t255 - pkin(8) * t250, -pkin(1) * t238 * t241 - pkin(8) * t249, t226 ^ 2 / 0.2e1, -t226 * t225, t226 * t231, -t225 * t231, t231 ^ 2 / 0.2e1, t225 * t229 + t231 * t246, t226 * t229 - t231 * t254, -t206 * t219 + t208 * t214, t207 * t219 + t208 * t215, -t206 * t215 - t207 * t214, t207 ^ 2 / 0.2e1 + t206 ^ 2 / 0.2e1 + t208 ^ 2 / 0.2e1, t204 * t215 + t205 * t214, -t203 * t214 - t204 * t219, -t203 * t215 + t205 * t219, t203 ^ 2 / 0.2e1 + t205 ^ 2 / 0.2e1 + t204 ^ 2 / 0.2e1, t200 * t215 - t202 * t214, -t201 * t215 - t202 * t219, t200 * t219 + t201 * t214, t201 ^ 2 / 0.2e1 + t200 ^ 2 / 0.2e1 + t202 ^ 2 / 0.2e1;];
T_reg  = t1;
