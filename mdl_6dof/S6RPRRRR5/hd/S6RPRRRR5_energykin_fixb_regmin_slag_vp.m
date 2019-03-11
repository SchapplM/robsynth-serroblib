% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:49
% EndTime: 2019-03-09 07:09:49
% DurationCPUTime: 0.14s
% Computational Cost: add. (513->55), mult. (1283->111), div. (0->0), fcn. (1016->10), ass. (0->47)
t209 = cos(qJ(3));
t208 = cos(qJ(5));
t207 = pkin(7) + qJ(2);
t188 = qJD(3) + qJD(4);
t189 = sin(pkin(11));
t190 = cos(pkin(11));
t194 = sin(qJ(3));
t179 = (t189 * t209 + t190 * t194) * qJD(1);
t203 = qJD(1) * t189;
t180 = t207 * t203;
t202 = qJD(1) * t190;
t181 = t207 * t202;
t199 = -t180 * t209 - t194 * t181;
t164 = qJD(3) * pkin(3) - t179 * pkin(8) + t199;
t178 = t194 * t203 - t202 * t209;
t204 = -t194 * t180 + t209 * t181;
t166 = -t178 * pkin(8) + t204;
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t205 = t193 * t164 + t196 * t166;
t155 = t188 * pkin(9) + t205;
t171 = t196 * t178 + t193 * t179;
t172 = -t193 * t178 + t196 * t179;
t182 = qJD(2) + (-pkin(2) * t190 - pkin(1)) * qJD(1);
t173 = t178 * pkin(3) + t182;
t158 = t171 * pkin(4) - t172 * pkin(9) + t173;
t192 = sin(qJ(5));
t206 = t208 * t155 + t192 * t158;
t201 = -t192 * t155 + t208 * t158;
t200 = t196 * t164 - t193 * t166;
t170 = qJD(5) + t171;
t154 = -t188 * pkin(4) - t200;
t197 = qJD(1) ^ 2;
t195 = cos(qJ(6));
t191 = sin(qJ(6));
t187 = t190 ^ 2;
t186 = t189 ^ 2;
t185 = -qJD(1) * pkin(1) + qJD(2);
t169 = qJD(6) + t170;
t168 = t172 * t208 + t192 * t188;
t167 = t192 * t172 - t188 * t208;
t160 = -t191 * t167 + t195 * t168;
t159 = t195 * t167 + t191 * t168;
t152 = t167 * pkin(5) + t154;
t151 = -t167 * pkin(10) + t206;
t150 = t170 * pkin(5) - t168 * pkin(10) + t201;
t1 = [t197 / 0.2e1, 0, 0, -t185 * t202, t185 * t203 (t186 + t187) * t197 * qJ(2), t185 ^ 2 / 0.2e1 + (t187 / 0.2e1 + t186 / 0.2e1) * qJ(2) ^ 2 * t197, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * qJD(3), -t178 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t199 + t182 * t178, -qJD(3) * t204 + t182 * t179, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t188, -t171 * t188, t188 ^ 2 / 0.2e1, t173 * t171 + t188 * t200, t173 * t172 - t188 * t205, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t170, -t167 * t170, t170 ^ 2 / 0.2e1, t154 * t167 + t170 * t201, t154 * t168 - t170 * t206, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t169, -t159 * t169, t169 ^ 2 / 0.2e1 (t195 * t150 - t191 * t151) * t169 + t152 * t159 -(t191 * t150 + t195 * t151) * t169 + t152 * t160;];
T_reg  = t1;
