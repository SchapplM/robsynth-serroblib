% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:01
% EndTime: 2019-03-09 01:10:01
% DurationCPUTime: 0.17s
% Computational Cost: add. (453->51), mult. (1069->112), div. (0->0), fcn. (884->14), ass. (0->51)
t244 = cos(qJ(2));
t257 = qJD(1) * sin(pkin(6));
t224 = qJD(2) * pkin(2) + t244 * t257;
t232 = sin(pkin(7));
t234 = cos(pkin(7));
t256 = qJD(1) * cos(pkin(6));
t246 = t224 * t234 + t232 * t256;
t262 = cos(qJ(5));
t245 = qJD(2) ^ 2;
t260 = t232 ^ 2 * t245;
t243 = cos(qJ(3));
t255 = qJD(2) * t232;
t250 = t243 * t255;
t227 = -qJD(4) + t250;
t230 = t234 * qJD(2) + qJD(3);
t240 = sin(qJ(2));
t222 = pkin(9) * t255 + t240 * t257;
t239 = sin(qJ(3));
t254 = t243 * t222 + t246 * t239;
t208 = t230 * pkin(10) + t254;
t229 = t234 * t256;
t211 = t229 + (-t224 + (-pkin(3) * t243 - pkin(10) * t239) * qJD(2)) * t232;
t238 = sin(qJ(4));
t242 = cos(qJ(4));
t258 = t242 * t208 + t238 * t211;
t199 = -t227 * pkin(11) + t258;
t220 = t239 * t222;
t207 = -t230 * pkin(3) - t246 * t243 + t220;
t252 = t239 * t255;
t217 = -t242 * t230 + t238 * t252;
t218 = t238 * t230 + t242 * t252;
t202 = t217 * pkin(4) - t218 * pkin(11) + t207;
t237 = sin(qJ(5));
t259 = t262 * t199 + t237 * t202;
t251 = (-t232 * t224 + t229) * t255;
t249 = qJD(2) * t257;
t248 = -t237 * t199 + t262 * t202;
t247 = -t238 * t208 + t242 * t211;
t198 = t227 * pkin(4) - t247;
t216 = qJD(5) + t217;
t241 = cos(qJ(6));
t236 = sin(qJ(6));
t214 = qJD(6) + t216;
t213 = t262 * t218 - t237 * t227;
t212 = t237 * t218 + t262 * t227;
t204 = -t236 * t212 + t241 * t213;
t203 = t241 * t212 + t236 * t213;
t196 = t212 * pkin(5) + t198;
t195 = -t212 * pkin(12) + t259;
t194 = t216 * pkin(5) - t213 * pkin(12) + t248;
t1 = [qJD(1) ^ 2 / 0.2e1, t245 / 0.2e1, t244 * t249, -t240 * t249, t239 ^ 2 * t260 / 0.2e1, t239 * t243 * t260, t230 * t252, t230 * t250, t230 ^ 2 / 0.2e1, -t220 * t230 + (t246 * t230 - t251) * t243, -t254 * t230 + t239 * t251, t218 ^ 2 / 0.2e1, -t218 * t217, -t218 * t227, t217 * t227, t227 ^ 2 / 0.2e1, t207 * t217 - t247 * t227, t207 * t218 + t258 * t227, t213 ^ 2 / 0.2e1, -t213 * t212, t213 * t216, -t212 * t216, t216 ^ 2 / 0.2e1, t198 * t212 + t248 * t216, t198 * t213 - t259 * t216, t204 ^ 2 / 0.2e1, -t204 * t203, t204 * t214, -t203 * t214, t214 ^ 2 / 0.2e1 (t241 * t194 - t236 * t195) * t214 + t196 * t203 -(t236 * t194 + t241 * t195) * t214 + t196 * t204;];
T_reg  = t1;
