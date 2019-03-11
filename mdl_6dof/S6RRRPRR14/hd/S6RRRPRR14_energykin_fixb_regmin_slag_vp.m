% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:20:35
% EndTime: 2019-03-09 20:20:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (517->55), mult. (1187->112), div. (0->0), fcn. (913->10), ass. (0->49)
t239 = pkin(3) + pkin(10);
t238 = cos(qJ(5));
t213 = sin(pkin(6));
t222 = qJD(1) ^ 2;
t237 = t213 ^ 2 * t222;
t232 = cos(pkin(6)) * qJD(1);
t211 = qJD(2) + t232;
t217 = sin(qJ(3));
t220 = cos(qJ(3));
t218 = sin(qJ(2));
t233 = qJD(1) * t213;
t229 = t218 * t233;
t203 = t217 * t211 + t220 * t229;
t221 = cos(qJ(2));
t228 = t221 * t233;
t206 = -qJD(3) + t228;
t231 = pkin(1) * t232;
t234 = pkin(8) * t228 + t218 * t231;
t197 = t211 * pkin(9) + t234;
t200 = (-pkin(2) * t221 - pkin(9) * t218 - pkin(1)) * t233;
t226 = -t217 * t197 + t220 * t200;
t225 = qJD(4) - t226;
t181 = t203 * pkin(4) + t239 * t206 + t225;
t202 = -t220 * t211 + t217 * t229;
t224 = -pkin(8) * t229 + t221 * t231;
t196 = -t211 * pkin(2) - t224;
t223 = -t203 * qJ(4) + t196;
t183 = t239 * t202 + t223;
t216 = sin(qJ(5));
t236 = t216 * t181 + t238 * t183;
t235 = t220 * t197 + t217 * t200;
t230 = t221 * t237;
t189 = t206 * qJ(4) - t235;
t227 = t238 * t181 - t216 * t183;
t184 = -t202 * pkin(4) - t189;
t201 = qJD(5) + t203;
t219 = cos(qJ(6));
t215 = sin(qJ(6));
t199 = qJD(6) + t201;
t192 = t216 * t202 - t238 * t206;
t191 = -t238 * t202 - t216 * t206;
t188 = t206 * pkin(3) + t225;
t187 = t202 * pkin(3) + t223;
t186 = -t215 * t191 + t219 * t192;
t185 = t219 * t191 + t215 * t192;
t178 = t191 * pkin(5) + t184;
t177 = -t191 * pkin(11) + t236;
t176 = t201 * pkin(5) - t192 * pkin(11) + t227;
t1 = [t222 / 0.2e1, 0, 0, t218 ^ 2 * t237 / 0.2e1, t218 * t230, t211 * t229, t211 * t228, t211 ^ 2 / 0.2e1, pkin(1) * t230 + t224 * t211, -pkin(1) * t218 * t237 - t234 * t211, t203 ^ 2 / 0.2e1, -t203 * t202, -t203 * t206, t202 * t206, t206 ^ 2 / 0.2e1, t196 * t202 - t226 * t206, t196 * t203 + t235 * t206, t188 * t203 + t189 * t202, -t187 * t202 - t188 * t206, -t187 * t203 + t189 * t206, t187 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t188 ^ 2 / 0.2e1, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t201, -t191 * t201, t201 ^ 2 / 0.2e1, t184 * t191 + t227 * t201, t184 * t192 - t236 * t201, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t199, -t185 * t199, t199 ^ 2 / 0.2e1 (t219 * t176 - t215 * t177) * t199 + t178 * t185 -(t215 * t176 + t219 * t177) * t199 + t178 * t186;];
T_reg  = t1;
