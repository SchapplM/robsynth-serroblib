% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:46
% EndTime: 2019-03-09 11:21:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (494->56), mult. (1160->115), div. (0->0), fcn. (850->10), ass. (0->49)
t235 = -pkin(2) - pkin(9);
t213 = sin(pkin(6));
t222 = qJD(1) ^ 2;
t234 = t213 ^ 2 * t222;
t215 = cos(pkin(6));
t230 = t215 * qJD(1);
t210 = qJD(2) + t230;
t217 = sin(qJ(4));
t220 = cos(qJ(4));
t221 = cos(qJ(2));
t231 = qJD(1) * t213;
t227 = t221 * t231;
t202 = t220 * t210 - t217 * t227;
t218 = sin(qJ(2));
t226 = t218 * t231;
t204 = qJD(4) + t226;
t206 = pkin(8) * t226;
t193 = qJD(3) + t206 + t235 * t210 + (-pkin(1) * t215 * t221 + pkin(3) * t213 * t218) * qJD(1);
t225 = -qJ(3) * t218 - pkin(1);
t196 = (t235 * t221 + t225) * t231;
t224 = t220 * t193 - t217 * t196;
t182 = t204 * pkin(4) - t202 * qJ(5) + t224;
t201 = t217 * t210 + t220 * t227;
t233 = t217 * t193 + t220 * t196;
t184 = -t201 * qJ(5) + t233;
t212 = sin(pkin(11));
t214 = cos(pkin(11));
t179 = t212 * t182 + t214 * t184;
t229 = pkin(1) * t230;
t232 = pkin(8) * t227 + t218 * t229;
t228 = t221 * t234;
t198 = -t210 * qJ(3) - t232;
t189 = -t214 * t201 - t212 * t202;
t195 = pkin(3) * t227 - t198;
t178 = t214 * t182 - t212 * t184;
t223 = t221 * t229 - t206;
t187 = t201 * pkin(4) + qJD(5) + t195;
t219 = cos(qJ(6));
t216 = sin(qJ(6));
t200 = (-pkin(2) * t221 + t225) * t231;
t197 = -t210 * pkin(2) + qJD(3) - t223;
t190 = -t212 * t201 + t214 * t202;
t188 = qJD(6) - t189;
t186 = t219 * t190 + t216 * t204;
t185 = t216 * t190 - t219 * t204;
t180 = -t189 * pkin(5) - t190 * pkin(10) + t187;
t177 = t204 * pkin(10) + t179;
t176 = -t204 * pkin(5) - t178;
t1 = [t222 / 0.2e1, 0, 0, t218 ^ 2 * t234 / 0.2e1, t218 * t228, t210 * t226, t210 * t227, t210 ^ 2 / 0.2e1, pkin(1) * t228 + t223 * t210, -pkin(1) * t218 * t234 - t232 * t210 (t197 * t218 - t198 * t221) * t231, t197 * t210 + t200 * t227, -t198 * t210 - t200 * t226, t200 ^ 2 / 0.2e1 + t198 ^ 2 / 0.2e1 + t197 ^ 2 / 0.2e1, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t204, -t201 * t204, t204 ^ 2 / 0.2e1, t195 * t201 + t224 * t204, t195 * t202 - t233 * t204, -t178 * t190 + t179 * t189, t179 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t188, -t185 * t188, t188 ^ 2 / 0.2e1 (-t216 * t177 + t219 * t180) * t188 + t176 * t185 -(t219 * t177 + t216 * t180) * t188 + t176 * t186;];
T_reg  = t1;
