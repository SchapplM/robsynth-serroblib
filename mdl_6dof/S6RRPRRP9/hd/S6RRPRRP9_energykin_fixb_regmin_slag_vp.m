% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:33:57
% EndTime: 2019-03-09 12:33:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (637->53), mult. (1547->109), div. (0->0), fcn. (1211->10), ass. (0->47)
t233 = cos(qJ(5));
t232 = cos(pkin(11));
t211 = sin(pkin(6));
t218 = qJD(1) ^ 2;
t231 = t211 ^ 2 * t218;
t217 = cos(qJ(2));
t227 = qJD(1) * t211;
t222 = t217 * t227;
t203 = -qJD(4) + t222;
t226 = cos(pkin(6)) * qJD(1);
t208 = qJD(2) + t226;
t215 = sin(qJ(2));
t225 = pkin(1) * t226;
t228 = pkin(8) * t222 + t215 * t225;
t196 = t208 * qJ(3) + t228;
t198 = (-pkin(2) * t217 - qJ(3) * t215 - pkin(1)) * t227;
t210 = sin(pkin(11));
t186 = -t210 * t196 + t232 * t198;
t223 = t215 * t227;
t200 = t210 * t208 + t232 * t223;
t180 = -pkin(3) * t222 - t200 * pkin(9) + t186;
t187 = t232 * t196 + t210 * t198;
t199 = -t232 * t208 + t210 * t223;
t183 = -t199 * pkin(9) + t187;
t214 = sin(qJ(4));
t216 = cos(qJ(4));
t229 = t214 * t180 + t216 * t183;
t175 = -t203 * pkin(10) + t229;
t219 = -pkin(8) * t223 + t217 * t225;
t193 = -t208 * pkin(2) + qJD(3) - t219;
t189 = t199 * pkin(3) + t193;
t190 = t216 * t199 + t214 * t200;
t191 = -t214 * t199 + t216 * t200;
t178 = t190 * pkin(4) - t191 * pkin(10) + t189;
t213 = sin(qJ(5));
t230 = t233 * t175 + t213 * t178;
t224 = t217 * t231;
t221 = -t213 * t175 + t233 * t178;
t220 = t216 * t180 - t214 * t183;
t174 = t203 * pkin(4) - t220;
t188 = qJD(5) + t190;
t185 = t233 * t191 - t213 * t203;
t184 = t213 * t191 + t233 * t203;
t172 = t184 * pkin(5) + qJD(6) + t174;
t171 = -t184 * qJ(6) + t230;
t170 = t188 * pkin(5) - t185 * qJ(6) + t221;
t1 = [t218 / 0.2e1, 0, 0, t215 ^ 2 * t231 / 0.2e1, t215 * t224, t208 * t223, t208 * t222, t208 ^ 2 / 0.2e1, pkin(1) * t224 + t219 * t208, -pkin(1) * t215 * t231 - t228 * t208, -t186 * t222 + t193 * t199, t187 * t222 + t193 * t200, -t186 * t200 - t187 * t199, t187 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t191 ^ 2 / 0.2e1, -t191 * t190, -t191 * t203, t190 * t203, t203 ^ 2 / 0.2e1, t189 * t190 - t220 * t203, t189 * t191 + t229 * t203, t185 ^ 2 / 0.2e1, -t185 * t184, t185 * t188, -t184 * t188, t188 ^ 2 / 0.2e1, t174 * t184 + t221 * t188, t174 * t185 - t230 * t188, -t170 * t185 - t171 * t184, t171 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1;];
T_reg  = t1;
