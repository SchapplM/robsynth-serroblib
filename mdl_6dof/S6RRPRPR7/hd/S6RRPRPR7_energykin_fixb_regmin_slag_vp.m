% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:42
% EndTime: 2019-03-09 10:47:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (334->51), mult. (728->106), div. (0->0), fcn. (466->8), ass. (0->42)
t202 = qJD(1) ^ 2;
t212 = t202 / 0.2e1;
t201 = cos(qJ(2));
t211 = t201 * t202;
t197 = sin(qJ(4));
t198 = sin(qJ(2));
t200 = cos(qJ(4));
t179 = (-t197 * t201 + t198 * t200) * qJD(1);
t190 = qJD(2) - qJD(4);
t209 = qJD(1) * t198;
t207 = pkin(7) * t209 + qJD(3);
t175 = -pkin(8) * t209 + (-pkin(2) - pkin(3)) * qJD(2) + t207;
t208 = qJD(1) * t201;
t183 = pkin(7) * t208 + qJD(2) * qJ(3);
t180 = -pkin(8) * t208 + t183;
t203 = t200 * t175 - t197 * t180;
t163 = -t190 * pkin(4) - t179 * qJ(5) + t203;
t178 = (t197 * t198 + t200 * t201) * qJD(1);
t210 = t197 * t175 + t200 * t180;
t165 = -t178 * qJ(5) + t210;
t193 = sin(pkin(10));
t194 = cos(pkin(10));
t160 = t193 * t163 + t194 * t165;
t206 = qJD(1) * qJD(2);
t181 = -qJD(1) * pkin(1) - pkin(2) * t208 - qJ(3) * t209;
t205 = t198 * t206;
t204 = t201 * t206;
t170 = -t194 * t178 - t193 * t179;
t174 = pkin(3) * t208 - t181;
t159 = t194 * t163 - t193 * t165;
t169 = t178 * pkin(4) + qJD(5) + t174;
t199 = cos(qJ(6));
t196 = sin(qJ(6));
t182 = -qJD(2) * pkin(2) + t207;
t171 = -t193 * t178 + t194 * t179;
t168 = qJD(6) - t170;
t167 = t199 * t171 - t196 * t190;
t166 = t196 * t171 + t199 * t190;
t161 = -t170 * pkin(5) - t171 * pkin(9) + t169;
t158 = -t190 * pkin(9) + t160;
t157 = t190 * pkin(5) - t159;
t1 = [t212, 0, 0, t198 ^ 2 * t212, t198 * t211, t205, t204, qJD(2) ^ 2 / 0.2e1, pkin(1) * t211 - pkin(7) * t205, -t202 * pkin(1) * t198 - pkin(7) * t204, -t182 * qJD(2) - t181 * t208 (t182 * t198 + t183 * t201) * qJD(1), t183 * qJD(2) - t181 * t209, t183 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t179 ^ 2 / 0.2e1, -t179 * t178, -t179 * t190, t178 * t190, t190 ^ 2 / 0.2e1, t174 * t178 - t203 * t190, t174 * t179 + t210 * t190, -t159 * t171 + t160 * t170, t160 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t168, -t166 * t168, t168 ^ 2 / 0.2e1 (-t196 * t158 + t199 * t161) * t168 + t157 * t166 -(t199 * t158 + t196 * t161) * t168 + t157 * t167;];
T_reg  = t1;
