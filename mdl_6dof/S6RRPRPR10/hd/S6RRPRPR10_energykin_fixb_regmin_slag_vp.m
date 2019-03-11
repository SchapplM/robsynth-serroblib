% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:24
% EndTime: 2019-03-09 11:09:24
% DurationCPUTime: 0.17s
% Computational Cost: add. (636->56), mult. (1552->113), div. (0->0), fcn. (1201->10), ass. (0->48)
t232 = pkin(4) + pkin(10);
t231 = cos(pkin(11));
t209 = sin(pkin(6));
t217 = qJD(1) ^ 2;
t230 = t209 ^ 2 * t217;
t226 = cos(pkin(6)) * qJD(1);
t206 = qJD(2) + t226;
t213 = sin(qJ(2));
t216 = cos(qJ(2));
t227 = qJD(1) * t209;
t222 = t216 * t227;
t225 = pkin(1) * t226;
t228 = pkin(8) * t222 + t213 * t225;
t195 = t206 * qJ(3) + t228;
t197 = (-pkin(2) * t216 - qJ(3) * t213 - pkin(1)) * t227;
t208 = sin(pkin(11));
t184 = -t208 * t195 + t231 * t197;
t223 = t213 * t227;
t199 = t208 * t206 + t231 * t223;
t178 = -pkin(3) * t222 - t199 * pkin(9) + t184;
t185 = t231 * t195 + t208 * t197;
t198 = -t231 * t206 + t208 * t223;
t181 = -t198 * pkin(9) + t185;
t212 = sin(qJ(4));
t215 = cos(qJ(4));
t229 = t212 * t178 + t215 * t181;
t224 = t216 * t230;
t221 = t215 * t178 - t212 * t181;
t201 = -qJD(4) + t222;
t175 = t201 * qJ(5) - t229;
t220 = qJD(5) - t221;
t190 = -t212 * t198 + t215 * t199;
t219 = -pkin(8) * t223 + t216 * t225;
t192 = -t206 * pkin(2) + qJD(3) - t219;
t188 = t198 * pkin(3) + t192;
t218 = -t190 * qJ(5) + t188;
t214 = cos(qJ(6));
t211 = sin(qJ(6));
t189 = t215 * t198 + t212 * t199;
t187 = qJD(6) + t190;
t183 = t211 * t189 - t214 * t201;
t182 = -t214 * t189 - t211 * t201;
t176 = t189 * pkin(4) + t218;
t174 = t201 * pkin(4) + t220;
t173 = t232 * t189 + t218;
t172 = -t189 * pkin(5) - t175;
t171 = t190 * pkin(5) + t232 * t201 + t220;
t1 = [t217 / 0.2e1, 0, 0, t213 ^ 2 * t230 / 0.2e1, t213 * t224, t206 * t223, t206 * t222, t206 ^ 2 / 0.2e1, pkin(1) * t224 + t219 * t206, -pkin(1) * t213 * t230 - t228 * t206, -t184 * t222 + t192 * t198, t185 * t222 + t192 * t199, -t184 * t199 - t185 * t198, t185 ^ 2 / 0.2e1 + t184 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1, t190 ^ 2 / 0.2e1, -t190 * t189, -t190 * t201, t189 * t201, t201 ^ 2 / 0.2e1, t188 * t189 - t221 * t201, t188 * t190 + t229 * t201, t174 * t190 + t175 * t189, -t174 * t201 - t176 * t189, t175 * t201 - t176 * t190, t176 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1, t183 ^ 2 / 0.2e1, -t183 * t182, t183 * t187, -t182 * t187, t187 ^ 2 / 0.2e1 (t214 * t171 - t211 * t173) * t187 + t172 * t182 -(t211 * t171 + t214 * t173) * t187 + t172 * t183;];
T_reg  = t1;
