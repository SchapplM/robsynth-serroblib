% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:20
% EndTime: 2019-03-08 23:53:20
% DurationCPUTime: 0.17s
% Computational Cost: add. (354->49), mult. (860->105), div. (0->0), fcn. (679->12), ass. (0->47)
t228 = cos(qJ(2));
t242 = qJD(1) * sin(pkin(6));
t209 = qJD(2) * pkin(2) + t228 * t242;
t217 = sin(pkin(7));
t219 = cos(pkin(7));
t241 = qJD(1) * cos(pkin(6));
t231 = t209 * t219 + t217 * t241;
t246 = pkin(4) + pkin(11);
t229 = qJD(2) ^ 2;
t244 = t217 ^ 2 * t229;
t215 = t219 * qJD(2) + qJD(3);
t224 = sin(qJ(2));
t240 = qJD(2) * t217;
t208 = pkin(9) * t240 + t224 * t242;
t223 = sin(qJ(3));
t227 = cos(qJ(3));
t239 = t227 * t208 + t231 * t223;
t195 = t215 * pkin(10) + t239;
t214 = t219 * t241;
t197 = t214 + (-t209 + (-pkin(3) * t227 - pkin(10) * t223) * qJD(2)) * t217;
t222 = sin(qJ(4));
t226 = cos(qJ(4));
t243 = t226 * t195 + t222 * t197;
t237 = t223 * t240;
t236 = (-t217 * t209 + t214) * t240;
t235 = t227 * t240;
t234 = qJD(2) * t242;
t233 = -t222 * t195 + t226 * t197;
t212 = -qJD(4) + t235;
t190 = t212 * qJ(5) - t243;
t232 = qJD(5) - t233;
t204 = t222 * t215 + t226 * t237;
t206 = t223 * t208;
t194 = -t215 * pkin(3) - t231 * t227 + t206;
t230 = -t204 * qJ(5) + t194;
t225 = cos(qJ(6));
t221 = sin(qJ(6));
t203 = -t226 * t215 + t222 * t237;
t202 = qJD(6) + t204;
t199 = t221 * t203 - t225 * t212;
t198 = -t225 * t203 - t221 * t212;
t191 = t203 * pkin(4) + t230;
t189 = t212 * pkin(4) + t232;
t188 = t246 * t203 + t230;
t187 = -t203 * pkin(5) - t190;
t186 = t204 * pkin(5) + t246 * t212 + t232;
t1 = [qJD(1) ^ 2 / 0.2e1, t229 / 0.2e1, t228 * t234, -t224 * t234, t223 ^ 2 * t244 / 0.2e1, t223 * t227 * t244, t215 * t237, t215 * t235, t215 ^ 2 / 0.2e1, -t206 * t215 + (t231 * t215 - t236) * t227, -t239 * t215 + t223 * t236, t204 ^ 2 / 0.2e1, -t204 * t203, -t204 * t212, t203 * t212, t212 ^ 2 / 0.2e1, t194 * t203 - t233 * t212, t194 * t204 + t243 * t212, t189 * t204 + t190 * t203, -t189 * t212 - t191 * t203, t190 * t212 - t191 * t204, t191 ^ 2 / 0.2e1 + t190 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1, t199 ^ 2 / 0.2e1, -t199 * t198, t199 * t202, -t198 * t202, t202 ^ 2 / 0.2e1 (t225 * t186 - t221 * t188) * t202 + t187 * t198 -(t221 * t186 + t225 * t188) * t202 + t187 * t199;];
T_reg  = t1;
