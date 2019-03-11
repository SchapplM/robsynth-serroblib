% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:41:42
% EndTime: 2019-03-09 14:41:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (481->57), mult. (1112->118), div. (0->0), fcn. (836->10), ass. (0->50)
t234 = -pkin(2) - pkin(9);
t209 = sin(pkin(6));
t219 = qJD(1) ^ 2;
t233 = t209 ^ 2 * t219;
t210 = cos(pkin(6));
t228 = t210 * qJD(1);
t207 = qJD(2) + t228;
t213 = sin(qJ(4));
t217 = cos(qJ(4));
t218 = cos(qJ(2));
t229 = qJD(1) * t209;
t225 = t218 * t229;
t198 = t217 * t207 - t213 * t225;
t214 = sin(qJ(2));
t224 = t214 * t229;
t201 = qJD(4) + t224;
t203 = pkin(8) * t224;
t189 = qJD(3) + t203 + t234 * t207 + (-pkin(1) * t210 * t218 + pkin(3) * t209 * t214) * qJD(1);
t223 = -qJ(3) * t214 - pkin(1);
t192 = (t234 * t218 + t223) * t229;
t222 = t217 * t189 - t213 * t192;
t178 = t201 * pkin(4) - t198 * pkin(10) + t222;
t197 = t213 * t207 + t217 * t225;
t231 = t213 * t189 + t217 * t192;
t180 = -t197 * pkin(10) + t231;
t212 = sin(qJ(5));
t216 = cos(qJ(5));
t232 = t212 * t178 + t216 * t180;
t227 = pkin(1) * t228;
t230 = pkin(8) * t225 + t214 * t227;
t226 = t218 * t233;
t194 = -t207 * qJ(3) - t230;
t185 = t216 * t197 + t212 * t198;
t191 = pkin(3) * t225 - t194;
t221 = t216 * t178 - t212 * t180;
t220 = t218 * t227 - t203;
t183 = t197 * pkin(4) + t191;
t215 = cos(qJ(6));
t211 = sin(qJ(6));
t200 = qJD(5) + t201;
t196 = (-pkin(2) * t218 + t223) * t229;
t193 = -t207 * pkin(2) + qJD(3) - t220;
t186 = -t212 * t197 + t216 * t198;
t184 = qJD(6) + t185;
t182 = t215 * t186 + t211 * t200;
t181 = t211 * t186 - t215 * t200;
t176 = t185 * pkin(5) - t186 * pkin(11) + t183;
t175 = t200 * pkin(11) + t232;
t174 = -t200 * pkin(5) - t221;
t1 = [t219 / 0.2e1, 0, 0, t214 ^ 2 * t233 / 0.2e1, t214 * t226, t207 * t224, t207 * t225, t207 ^ 2 / 0.2e1, pkin(1) * t226 + t220 * t207, -pkin(1) * t214 * t233 - t230 * t207 (t193 * t214 - t194 * t218) * t229, t193 * t207 + t196 * t225, -t194 * t207 - t196 * t224, t196 ^ 2 / 0.2e1 + t194 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t198 ^ 2 / 0.2e1, -t198 * t197, t198 * t201, -t197 * t201, t201 ^ 2 / 0.2e1, t191 * t197 + t222 * t201, t191 * t198 - t231 * t201, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t200, -t185 * t200, t200 ^ 2 / 0.2e1, t183 * t185 + t221 * t200, t183 * t186 - t232 * t200, t182 ^ 2 / 0.2e1, -t182 * t181, t182 * t184, -t181 * t184, t184 ^ 2 / 0.2e1 (-t211 * t175 + t215 * t176) * t184 + t174 * t181 -(t215 * t175 + t211 * t176) * t184 + t174 * t182;];
T_reg  = t1;
