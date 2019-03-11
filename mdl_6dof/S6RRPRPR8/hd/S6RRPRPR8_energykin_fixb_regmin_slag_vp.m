% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:19
% EndTime: 2019-03-09 10:53:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (440->53), mult. (987->107), div. (0->0), fcn. (682->8), ass. (0->44)
t211 = pkin(4) + pkin(5);
t197 = qJD(1) ^ 2;
t210 = t197 / 0.2e1;
t196 = cos(qJ(2));
t209 = t196 * t197;
t193 = sin(qJ(2));
t177 = (-pkin(2) * t196 - qJ(3) * t193 - pkin(1)) * qJD(1);
t206 = t196 * qJD(1);
t182 = pkin(7) * t206 + qJD(2) * qJ(3);
t189 = sin(pkin(10));
t190 = cos(pkin(10));
t171 = t177 * t190 - t182 * t189;
t207 = qJD(1) * t193;
t179 = qJD(2) * t189 + t190 * t207;
t164 = -pkin(3) * t206 - pkin(8) * t179 + t171;
t172 = t177 * t189 + t182 * t190;
t178 = -qJD(2) * t190 + t189 * t207;
t167 = -pkin(8) * t178 + t172;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t208 = t164 * t192 + t167 * t195;
t205 = qJD(1) * qJD(2);
t185 = -qJD(4) + t206;
t159 = -qJ(5) * t185 + t208;
t204 = t193 * t205;
t203 = t196 * t205;
t202 = t195 * t164 - t167 * t192;
t201 = qJD(2) * pkin(2) - pkin(7) * t207 - qJD(3);
t200 = qJD(5) - t202;
t199 = -t178 * pkin(3) + t201;
t170 = -t178 * t192 + t179 * t195;
t198 = t170 * qJ(5) + t199;
t194 = cos(qJ(6));
t191 = sin(qJ(6));
t184 = qJD(6) + t185;
t169 = t178 * t195 + t179 * t192;
t162 = t169 * t191 + t170 * t194;
t161 = -t169 * t194 + t170 * t191;
t160 = pkin(4) * t169 - t198;
t158 = pkin(4) * t185 + t200;
t157 = -t169 * t211 + t198;
t156 = pkin(9) * t169 + t159;
t155 = -t170 * pkin(9) + t185 * t211 + t200;
t1 = [t210, 0, 0, t193 ^ 2 * t210, t193 * t209, t204, t203, qJD(2) ^ 2 / 0.2e1, pkin(1) * t209 - pkin(7) * t204, -pkin(1) * t193 * t197 - pkin(7) * t203, -t171 * t206 - t178 * t201, t172 * t206 - t179 * t201, -t171 * t179 - t172 * t178, t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1 + t201 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, -t170 * t185, t169 * t185, t185 ^ 2 / 0.2e1, -t169 * t199 - t185 * t202, -t170 * t199 + t185 * t208, t158 * t185 + t160 * t169, t158 * t170 - t159 * t169, -t159 * t185 - t160 * t170, t159 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t184, -t161 * t184, t184 ^ 2 / 0.2e1 (t155 * t194 - t156 * t191) * t184 + t157 * t161 -(t155 * t191 + t156 * t194) * t184 + t157 * t162;];
T_reg  = t1;
