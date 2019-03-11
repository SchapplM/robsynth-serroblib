% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:58
% EndTime: 2019-03-09 09:10:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (297->51), mult. (718->108), div. (0->0), fcn. (484->8), ass. (0->39)
t190 = sin(pkin(6));
t198 = qJD(1) ^ 2;
t209 = t190 ^ 2 * t198;
t197 = cos(qJ(2));
t205 = qJD(1) * t190;
t200 = t197 * t205;
t194 = sin(qJ(2));
t201 = t194 * t205;
t169 = -pkin(1) * t205 - pkin(2) * t200 - qJ(3) * t201;
t166 = pkin(3) * t200 + qJD(4) - t169;
t159 = (pkin(4) * t197 - pkin(9) * t194) * t205 + t166;
t204 = cos(pkin(6)) * qJD(1);
t187 = qJD(2) + t204;
t203 = pkin(1) * t204;
t206 = pkin(8) * t200 + t194 * t203;
t168 = t187 * qJ(3) + t206;
t165 = -qJ(4) * t200 + t168;
t162 = -t187 * pkin(9) + t165;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t208 = t193 * t159 + t196 * t162;
t207 = -pkin(8) * t201 + t197 * t203;
t202 = t197 * t209;
t167 = -t187 * pkin(2) + qJD(3) - t207;
t199 = t196 * t159 - t193 * t162;
t171 = t196 * t187 + t193 * t201;
t161 = -t187 * pkin(3) - qJ(4) * t201 + t167;
t158 = pkin(4) * t187 - t161;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t175 = qJD(5) + t200;
t172 = -t193 * t187 + t196 * t201;
t170 = qJD(6) + t171;
t164 = t195 * t172 + t192 * t175;
t163 = t172 * t192 - t195 * t175;
t156 = pkin(5) * t171 - pkin(10) * t172 + t158;
t155 = pkin(10) * t175 + t208;
t154 = -t175 * pkin(5) - t199;
t1 = [t198 / 0.2e1, 0, 0, t194 ^ 2 * t209 / 0.2e1, t194 * t202, t187 * t201, t187 * t200, t187 ^ 2 / 0.2e1, pkin(1) * t202 + t207 * t187, -pkin(1) * t194 * t209 - t206 * t187, -t167 * t187 - t169 * t200 (t167 * t194 + t168 * t197) * t205, t168 * t187 - t169 * t201, t168 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1, -t161 * t187 + t166 * t200, t165 * t187 + t166 * t201 (-t161 * t194 - t165 * t197) * t205, t165 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t175, -t171 * t175, t175 ^ 2 / 0.2e1, t158 * t171 + t199 * t175, t158 * t172 - t208 * t175, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t170, -t163 * t170, t170 ^ 2 / 0.2e1 (-t192 * t155 + t195 * t156) * t170 + t154 * t163 -(t195 * t155 + t192 * t156) * t170 + t154 * t164;];
T_reg  = t1;
