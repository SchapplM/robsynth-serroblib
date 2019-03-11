% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:16
% EndTime: 2019-03-08 19:41:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (339->50), mult. (826->102), div. (0->0), fcn. (634->12), ass. (0->42)
t205 = cos(pkin(12));
t193 = sin(qJ(2));
t203 = qJD(1) * sin(pkin(6));
t181 = qJD(2) * qJ(3) + t193 * t203;
t189 = cos(pkin(11));
t202 = qJD(1) * cos(pkin(6));
t183 = t189 * t202;
t187 = sin(pkin(11));
t168 = t183 + (-pkin(8) * qJD(2) - t181) * t187;
t174 = t189 * t181 + t187 * t202;
t200 = qJD(2) * t189;
t169 = pkin(8) * t200 + t174;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t204 = t192 * t168 + t195 * t169;
t159 = qJD(4) * qJ(5) + t204;
t196 = cos(qJ(2));
t197 = -t196 * t203 + qJD(3);
t175 = (-pkin(3) * t189 - pkin(2)) * qJD(2) + t197;
t201 = qJD(2) * t187;
t177 = t192 * t201 - t195 * t200;
t178 = (t187 * t195 + t189 * t192) * qJD(2);
t164 = t177 * pkin(4) - t178 * qJ(5) + t175;
t186 = sin(pkin(12));
t155 = t205 * t159 + t186 * t164;
t199 = qJD(2) * t203;
t154 = -t186 * t159 + t205 * t164;
t198 = t195 * t168 - t192 * t169;
t158 = -qJD(4) * pkin(4) + qJD(5) - t198;
t194 = cos(qJ(6));
t191 = sin(qJ(6));
t180 = -qJD(2) * pkin(2) + t197;
t176 = qJD(6) + t177;
t173 = -t187 * t181 + t183;
t172 = t186 * qJD(4) + t205 * t178;
t171 = -t205 * qJD(4) + t186 * t178;
t161 = -t191 * t171 + t194 * t172;
t160 = t194 * t171 + t191 * t172;
t156 = t171 * pkin(5) + t158;
t153 = -t171 * pkin(9) + t155;
t152 = t177 * pkin(5) - t172 * pkin(9) + t154;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t196 * t199, -t193 * t199, -t180 * t200, t180 * t201 (-t173 * t187 + t174 * t189) * qJD(2), t174 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1, t178 ^ 2 / 0.2e1, -t178 * t177, t178 * qJD(4), -t177 * qJD(4), qJD(4) ^ 2 / 0.2e1, t198 * qJD(4) + t175 * t177, -t204 * qJD(4) + t175 * t178, t154 * t177 + t158 * t171, -t155 * t177 + t158 * t172, -t154 * t172 - t155 * t171, t155 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t176, -t160 * t176, t176 ^ 2 / 0.2e1 (t194 * t152 - t191 * t153) * t176 + t156 * t160 -(t191 * t152 + t194 * t153) * t176 + t156 * t161;];
T_reg  = t1;
