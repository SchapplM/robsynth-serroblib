% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:52
% EndTime: 2019-03-09 00:25:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (360->46), mult. (873->101), div. (0->0), fcn. (699->12), ass. (0->46)
t226 = cos(qJ(2));
t239 = qJD(1) * sin(pkin(6));
t208 = qJD(2) * pkin(2) + t226 * t239;
t216 = sin(pkin(7));
t218 = cos(pkin(7));
t238 = qJD(1) * cos(pkin(6));
t228 = t208 * t218 + t216 * t238;
t244 = cos(qJ(5));
t227 = qJD(2) ^ 2;
t242 = t216 ^ 2 * t227;
t225 = cos(qJ(3));
t237 = qJD(2) * t216;
t232 = t225 * t237;
t211 = -qJD(4) + t232;
t214 = t218 * qJD(2) + qJD(3);
t223 = sin(qJ(2));
t206 = pkin(9) * t237 + t223 * t239;
t222 = sin(qJ(3));
t236 = t225 * t206 + t228 * t222;
t194 = t214 * pkin(10) + t236;
t213 = t218 * t238;
t196 = t213 + (-t208 + (-pkin(3) * t225 - pkin(10) * t222) * qJD(2)) * t216;
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t240 = t224 * t194 + t221 * t196;
t187 = -t211 * pkin(11) + t240;
t204 = t222 * t206;
t193 = -t214 * pkin(3) - t228 * t225 + t204;
t234 = t222 * t237;
t201 = -t224 * t214 + t221 * t234;
t202 = t221 * t214 + t224 * t234;
t190 = t201 * pkin(4) - t202 * pkin(11) + t193;
t220 = sin(qJ(5));
t241 = t244 * t187 + t220 * t190;
t233 = (-t216 * t208 + t213) * t237;
t231 = qJD(2) * t239;
t230 = -t220 * t187 + t244 * t190;
t229 = -t221 * t194 + t224 * t196;
t186 = t211 * pkin(4) - t229;
t200 = qJD(5) + t201;
t198 = t244 * t202 - t220 * t211;
t197 = t220 * t202 + t244 * t211;
t184 = t197 * pkin(5) + qJD(6) + t186;
t183 = -t197 * qJ(6) + t241;
t182 = t200 * pkin(5) - t198 * qJ(6) + t230;
t1 = [qJD(1) ^ 2 / 0.2e1, t227 / 0.2e1, t226 * t231, -t223 * t231, t222 ^ 2 * t242 / 0.2e1, t222 * t225 * t242, t214 * t234, t214 * t232, t214 ^ 2 / 0.2e1, -t204 * t214 + (t228 * t214 - t233) * t225, -t236 * t214 + t222 * t233, t202 ^ 2 / 0.2e1, -t202 * t201, -t202 * t211, t201 * t211, t211 ^ 2 / 0.2e1, t193 * t201 - t229 * t211, t193 * t202 + t240 * t211, t198 ^ 2 / 0.2e1, -t198 * t197, t198 * t200, -t197 * t200, t200 ^ 2 / 0.2e1, t186 * t197 + t230 * t200, t186 * t198 - t241 * t200, -t182 * t198 - t183 * t197, t183 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t184 ^ 2 / 0.2e1;];
T_reg  = t1;
