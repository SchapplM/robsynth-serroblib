% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:29
% EndTime: 2019-03-09 18:29:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (616->53), mult. (1375->110), div. (0->0), fcn. (1049->10), ass. (0->48)
t206 = qJD(1) ^ 2;
t219 = t206 / 0.2e1;
t218 = cos(qJ(3));
t217 = cos(qJ(5));
t205 = cos(qJ(2));
t216 = t205 * t206;
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t213 = qJD(1) * t203;
t186 = t202 * qJD(2) + t218 * t213;
t212 = t205 * qJD(1);
t194 = -qJD(3) + t212;
t184 = (-pkin(2) * t205 - pkin(8) * t203 - pkin(1)) * qJD(1);
t191 = pkin(7) * t212 + qJD(2) * pkin(8);
t207 = t218 * t184 - t202 * t191;
t175 = -t194 * pkin(3) - t186 * qJ(4) + t207;
t185 = -t218 * qJD(2) + t202 * t213;
t214 = t202 * t184 + t218 * t191;
t177 = -t185 * qJ(4) + t214;
t198 = sin(pkin(11));
t199 = cos(pkin(11));
t167 = t199 * t175 - t198 * t177;
t180 = -t198 * t185 + t199 * t186;
t163 = -t194 * pkin(4) - t180 * pkin(9) + t167;
t168 = t198 * t175 + t199 * t177;
t179 = -t199 * t185 - t198 * t186;
t165 = t179 * pkin(9) + t168;
t201 = sin(qJ(5));
t215 = t201 * t163 + t217 * t165;
t211 = qJD(1) * qJD(2);
t210 = t203 * t211;
t209 = t205 * t211;
t190 = -qJD(2) * pkin(2) + pkin(7) * t213;
t208 = t217 * t163 - t201 * t165;
t192 = -qJD(5) + t194;
t181 = t185 * pkin(3) + qJD(4) + t190;
t172 = -t179 * pkin(4) + t181;
t204 = cos(qJ(6));
t200 = sin(qJ(6));
t188 = -qJD(6) + t192;
t171 = t201 * t179 + t217 * t180;
t170 = -t217 * t179 + t201 * t180;
t166 = t170 * pkin(5) + t172;
t160 = -t200 * t170 + t204 * t171;
t159 = t204 * t170 + t200 * t171;
t158 = -t170 * pkin(10) + t215;
t157 = -t192 * pkin(5) - t171 * pkin(10) + t208;
t1 = [t219, 0, 0, t203 ^ 2 * t219, t203 * t216, t210, t209, qJD(2) ^ 2 / 0.2e1, pkin(1) * t216 - pkin(7) * t210, -t206 * pkin(1) * t203 - pkin(7) * t209, t186 ^ 2 / 0.2e1, -t186 * t185, -t186 * t194, t185 * t194, t194 ^ 2 / 0.2e1, t190 * t185 - t207 * t194, t190 * t186 + t214 * t194, -t167 * t180 + t168 * t179, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t192, t170 * t192, t192 ^ 2 / 0.2e1, t172 * t170 - t208 * t192, t172 * t171 + t215 * t192, t160 ^ 2 / 0.2e1, -t160 * t159, -t160 * t188, t159 * t188, t188 ^ 2 / 0.2e1 -(t204 * t157 - t200 * t158) * t188 + t166 * t159 (t200 * t157 + t204 * t158) * t188 + t166 * t160;];
T_reg  = t1;
