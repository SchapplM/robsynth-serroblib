% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:29
% EndTime: 2019-03-09 14:05:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (610->55), mult. (1398->114), div. (0->0), fcn. (1069->10), ass. (0->48)
t205 = qJD(1) ^ 2;
t218 = t205 / 0.2e1;
t217 = cos(qJ(5));
t216 = cos(pkin(11));
t204 = cos(qJ(2));
t215 = t204 * t205;
t197 = sin(pkin(11));
t201 = sin(qJ(2));
t212 = qJD(1) * t201;
t184 = -t216 * qJD(2) + t197 * t212;
t185 = t197 * qJD(2) + t216 * t212;
t200 = sin(qJ(4));
t203 = cos(qJ(4));
t176 = -t200 * t184 + t203 * t185;
t211 = t204 * qJD(1);
t193 = -qJD(4) + t211;
t183 = (-pkin(2) * t204 - qJ(3) * t201 - pkin(1)) * qJD(1);
t190 = pkin(7) * t211 + qJD(2) * qJ(3);
t177 = t216 * t183 - t197 * t190;
t171 = -pkin(3) * t211 - t185 * pkin(8) + t177;
t178 = t197 * t183 + t216 * t190;
t173 = -t184 * pkin(8) + t178;
t206 = t203 * t171 - t200 * t173;
t161 = -t193 * pkin(4) - t176 * pkin(9) + t206;
t175 = t203 * t184 + t200 * t185;
t213 = t200 * t171 + t203 * t173;
t163 = -t175 * pkin(9) + t213;
t199 = sin(qJ(5));
t214 = t199 * t161 + t217 * t163;
t210 = qJD(1) * qJD(2);
t209 = t201 * t210;
t208 = t204 * t210;
t207 = t217 * t161 - t199 * t163;
t187 = -qJD(2) * pkin(2) + pkin(7) * t212 + qJD(3);
t191 = -qJD(5) + t193;
t179 = t184 * pkin(3) + t187;
t168 = t175 * pkin(4) + t179;
t202 = cos(qJ(6));
t198 = sin(qJ(6));
t188 = -qJD(6) + t191;
t167 = -t199 * t175 + t217 * t176;
t166 = t217 * t175 + t199 * t176;
t164 = t166 * pkin(5) + t168;
t158 = -t198 * t166 + t202 * t167;
t157 = t202 * t166 + t198 * t167;
t156 = -t166 * pkin(10) + t214;
t155 = -t191 * pkin(5) - t167 * pkin(10) + t207;
t1 = [t218, 0, 0, t201 ^ 2 * t218, t201 * t215, t209, t208, qJD(2) ^ 2 / 0.2e1, pkin(1) * t215 - pkin(7) * t209, -t205 * pkin(1) * t201 - pkin(7) * t208, -t177 * t211 + t187 * t184, t178 * t211 + t187 * t185, -t177 * t185 - t178 * t184, t178 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t176 ^ 2 / 0.2e1, -t176 * t175, -t176 * t193, t175 * t193, t193 ^ 2 / 0.2e1, t179 * t175 - t206 * t193, t179 * t176 + t213 * t193, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t191, t166 * t191, t191 ^ 2 / 0.2e1, t168 * t166 - t207 * t191, t168 * t167 + t214 * t191, t158 ^ 2 / 0.2e1, -t158 * t157, -t158 * t188, t157 * t188, t188 ^ 2 / 0.2e1 -(t202 * t155 - t198 * t156) * t188 + t164 * t157 (t198 * t155 + t202 * t156) * t188 + t164 * t158;];
T_reg  = t1;
