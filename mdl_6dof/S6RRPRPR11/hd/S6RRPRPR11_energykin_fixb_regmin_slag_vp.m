% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR11
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:02
% EndTime: 2019-03-09 11:15:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (355->51), mult. (744->105), div. (0->0), fcn. (469->8), ass. (0->44)
t209 = -pkin(2) - pkin(8);
t197 = qJD(1) ^ 2;
t208 = t197 / 0.2e1;
t196 = cos(qJ(2));
t207 = t196 * t197;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t205 = qJD(1) * t196;
t180 = t195 * qJD(2) - t192 * t205;
t193 = sin(qJ(2));
t204 = t193 * qJD(1);
t184 = qJD(4) + t204;
t199 = -qJ(3) * t193 - pkin(1);
t175 = (t209 * t196 + t199) * qJD(1);
t203 = pkin(7) * t204 + qJD(3);
t176 = pkin(3) * t204 + t209 * qJD(2) + t203;
t198 = -t192 * t175 + t195 * t176;
t164 = t184 * pkin(4) - t180 * qJ(5) + t198;
t179 = t192 * qJD(2) + t195 * t205;
t206 = t195 * t175 + t192 * t176;
t167 = -t179 * qJ(5) + t206;
t189 = sin(pkin(10));
t190 = cos(pkin(10));
t159 = t189 * t164 + t190 * t167;
t182 = -pkin(7) * t205 - qJD(2) * qJ(3);
t202 = qJD(1) * qJD(2);
t177 = pkin(3) * t205 - t182;
t201 = t193 * t202;
t200 = t196 * t202;
t158 = t190 * t164 - t189 * t167;
t171 = t179 * pkin(4) + qJD(5) + t177;
t194 = cos(qJ(6));
t191 = sin(qJ(6));
t183 = qJD(6) + t184;
t181 = -qJD(2) * pkin(2) + t203;
t178 = (-pkin(2) * t196 + t199) * qJD(1);
t170 = -t189 * t179 + t190 * t180;
t169 = -t190 * t179 - t189 * t180;
t165 = -t169 * pkin(5) + t171;
t161 = t191 * t169 + t194 * t170;
t160 = -t194 * t169 + t191 * t170;
t157 = t169 * pkin(9) + t159;
t156 = t184 * pkin(5) - t170 * pkin(9) + t158;
t1 = [t208, 0, 0, t193 ^ 2 * t208, t193 * t207, t201, t200, qJD(2) ^ 2 / 0.2e1, pkin(1) * t207 - pkin(7) * t201, -t197 * pkin(1) * t193 - pkin(7) * t200 (t181 * t193 - t182 * t196) * qJD(1), t181 * qJD(2) + t178 * t205, -t182 * qJD(2) - t178 * t204, t178 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t180 ^ 2 / 0.2e1, -t180 * t179, t180 * t184, -t179 * t184, t184 ^ 2 / 0.2e1, t177 * t179 + t198 * t184, t177 * t180 - t206 * t184, -t158 * t170 + t159 * t169, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t183, -t160 * t183, t183 ^ 2 / 0.2e1 (t194 * t156 - t191 * t157) * t183 + t165 * t160 -(t191 * t156 + t194 * t157) * t183 + t165 * t161;];
T_reg  = t1;
