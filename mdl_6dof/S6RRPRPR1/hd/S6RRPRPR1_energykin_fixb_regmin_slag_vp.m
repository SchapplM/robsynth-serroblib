% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:56
% EndTime: 2019-03-09 10:09:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (640->53), mult. (1544->110), div. (0->0), fcn. (1162->10), ass. (0->46)
t213 = qJD(1) * (pkin(7) + qJ(3));
t203 = qJD(1) ^ 2;
t212 = t203 / 0.2e1;
t210 = cos(pkin(11));
t202 = cos(qJ(2));
t209 = t202 * t203;
t193 = qJD(2) + qJD(4);
t199 = sin(qJ(2));
t189 = qJD(2) * pkin(2) - t199 * t213;
t190 = t202 * t213;
t195 = sin(pkin(10));
t196 = cos(pkin(10));
t180 = t196 * t189 - t195 * t190;
t187 = (t195 * t202 + t196 * t199) * qJD(1);
t172 = qJD(2) * pkin(3) - t187 * pkin(8) + t180;
t181 = t195 * t189 + t196 * t190;
t186 = (-t195 * t199 + t196 * t202) * qJD(1);
t174 = t186 * pkin(8) + t181;
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t208 = t198 * t172 + t201 * t174;
t163 = t193 * qJ(5) + t208;
t178 = -t201 * t186 + t198 * t187;
t179 = t198 * t186 + t201 * t187;
t191 = qJD(3) + (-pkin(2) * t202 - pkin(1)) * qJD(1);
t182 = -t186 * pkin(3) + t191;
t166 = t178 * pkin(4) - t179 * qJ(5) + t182;
t194 = sin(pkin(11));
t159 = t210 * t163 + t194 * t166;
t207 = qJD(1) * qJD(2);
t206 = t199 * t207;
t205 = t202 * t207;
t158 = -t194 * t163 + t210 * t166;
t204 = t201 * t172 - t198 * t174;
t162 = -t193 * pkin(4) + qJD(5) - t204;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t177 = qJD(6) + t178;
t176 = t210 * t179 + t194 * t193;
t175 = t194 * t179 - t210 * t193;
t168 = -t197 * t175 + t200 * t176;
t167 = t200 * t175 + t197 * t176;
t160 = t175 * pkin(5) + t162;
t157 = -t175 * pkin(9) + t159;
t156 = t178 * pkin(5) - t176 * pkin(9) + t158;
t1 = [t212, 0, 0, t199 ^ 2 * t212, t199 * t209, t206, t205, qJD(2) ^ 2 / 0.2e1, pkin(1) * t209 - pkin(7) * t206, -t203 * pkin(1) * t199 - pkin(7) * t205, -t180 * t187 + t181 * t186, t181 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t193, -t178 * t193, t193 ^ 2 / 0.2e1, t182 * t178 + t204 * t193, t182 * t179 - t208 * t193, t158 * t178 + t162 * t175, -t159 * t178 + t162 * t176, -t158 * t176 - t159 * t175, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t177, -t167 * t177, t177 ^ 2 / 0.2e1 (t200 * t156 - t197 * t157) * t177 + t160 * t167 -(t197 * t156 + t200 * t157) * t177 + t160 * t168;];
T_reg  = t1;
