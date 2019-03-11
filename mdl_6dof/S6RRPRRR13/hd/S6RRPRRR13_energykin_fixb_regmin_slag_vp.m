% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR13
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:18
% EndTime: 2019-03-09 14:53:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (460->57), mult. (1065->118), div. (0->0), fcn. (799->10), ass. (0->50)
t232 = -pkin(2) - pkin(9);
t231 = cos(qJ(5));
t207 = sin(pkin(6));
t216 = qJD(1) ^ 2;
t230 = t207 ^ 2 * t216;
t212 = sin(qJ(2));
t226 = qJD(1) * t207;
t222 = t212 * t226;
t199 = qJD(4) + t222;
t201 = pkin(8) * t222;
t208 = cos(pkin(6));
t225 = qJD(1) * t208;
t205 = qJD(2) + t225;
t215 = cos(qJ(2));
t183 = qJD(3) + t201 + t232 * t205 + (-pkin(1) * t208 * t215 + pkin(3) * t207 * t212) * qJD(1);
t220 = -qJ(3) * t212 - pkin(1);
t190 = (t232 * t215 + t220) * t226;
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t228 = t211 * t183 + t214 * t190;
t176 = pkin(10) * t199 + t228;
t221 = t215 * t226;
t224 = pkin(1) * t225;
t227 = pkin(8) * t221 + t212 * t224;
t192 = -t205 * qJ(3) - t227;
t189 = pkin(3) * t221 - t192;
t196 = t205 * t211 + t214 * t221;
t197 = t205 * t214 - t211 * t221;
t181 = pkin(4) * t196 - pkin(10) * t197 + t189;
t210 = sin(qJ(5));
t229 = t231 * t176 + t210 * t181;
t223 = t215 * t230;
t219 = -t176 * t210 + t231 * t181;
t218 = t183 * t214 - t211 * t190;
t217 = t215 * t224 - t201;
t175 = -pkin(4) * t199 - t218;
t195 = qJD(5) + t196;
t213 = cos(qJ(6));
t209 = sin(qJ(6));
t194 = (-pkin(2) * t215 + t220) * t226;
t193 = qJD(6) + t195;
t191 = -t205 * pkin(2) + qJD(3) - t217;
t186 = t231 * t197 + t210 * t199;
t185 = t210 * t197 - t231 * t199;
t180 = -t185 * t209 + t186 * t213;
t179 = t213 * t185 + t186 * t209;
t173 = pkin(5) * t185 + t175;
t172 = -pkin(11) * t185 + t229;
t171 = pkin(5) * t195 - pkin(11) * t186 + t219;
t1 = [t216 / 0.2e1, 0, 0, t212 ^ 2 * t230 / 0.2e1, t212 * t223, t205 * t222, t205 * t221, t205 ^ 2 / 0.2e1, pkin(1) * t223 + t217 * t205, -pkin(1) * t212 * t230 - t227 * t205 (t191 * t212 - t192 * t215) * t226, t191 * t205 + t194 * t221, -t192 * t205 - t194 * t222, t194 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1, t197 ^ 2 / 0.2e1, -t197 * t196, t197 * t199, -t196 * t199, t199 ^ 2 / 0.2e1, t189 * t196 + t218 * t199, t189 * t197 - t228 * t199, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t195, -t185 * t195, t195 ^ 2 / 0.2e1, t175 * t185 + t219 * t195, t175 * t186 - t229 * t195, t180 ^ 2 / 0.2e1, -t180 * t179, t180 * t193, -t179 * t193, t193 ^ 2 / 0.2e1 (t171 * t213 - t172 * t209) * t193 + t173 * t179 -(t171 * t209 + t172 * t213) * t193 + t173 * t180;];
T_reg  = t1;
