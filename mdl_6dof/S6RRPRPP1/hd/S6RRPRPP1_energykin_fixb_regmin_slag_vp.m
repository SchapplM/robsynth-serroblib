% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:51
% EndTime: 2019-03-09 09:47:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (560->48), mult. (1299->99), div. (0->0), fcn. (918->8), ass. (0->43)
t197 = qJD(1) ^ 2;
t208 = t197 / 0.2e1;
t207 = cos(qJ(4));
t206 = pkin(7) + qJ(3);
t196 = cos(qJ(2));
t205 = t196 * t197;
t191 = sin(pkin(9));
t193 = cos(pkin(9));
t195 = sin(qJ(2));
t182 = (t191 * t196 + t193 * t195) * qJD(1);
t194 = sin(qJ(4));
t178 = t194 * qJD(2) + t207 * t182;
t202 = qJD(1) * t196;
t203 = qJD(1) * t195;
t181 = -t191 * t203 + t193 * t202;
t180 = qJD(4) - t181;
t187 = qJD(3) + (-pkin(2) * t196 - pkin(1)) * qJD(1);
t170 = -t181 * pkin(3) - t182 * pkin(8) + t187;
t185 = qJD(2) * pkin(2) - t206 * t203;
t186 = t206 * t202;
t175 = t191 * t185 + t193 * t186;
t173 = qJD(2) * pkin(8) + t175;
t198 = t207 * t170 - t194 * t173;
t162 = t180 * pkin(4) - t178 * qJ(5) + t198;
t177 = -t207 * qJD(2) + t194 * t182;
t204 = t194 * t170 + t207 * t173;
t164 = -t177 * qJ(5) + t204;
t190 = sin(pkin(10));
t192 = cos(pkin(10));
t159 = t190 * t162 + t192 * t164;
t201 = qJD(1) * qJD(2);
t200 = t195 * t201;
t199 = t196 * t201;
t174 = t193 * t185 - t191 * t186;
t158 = t192 * t162 - t190 * t164;
t172 = -qJD(2) * pkin(3) - t174;
t165 = t177 * pkin(4) + qJD(5) + t172;
t167 = -t190 * t177 + t192 * t178;
t166 = t192 * t177 + t190 * t178;
t160 = t166 * pkin(5) - t167 * qJ(6) + t165;
t157 = t180 * qJ(6) + t159;
t156 = -t180 * pkin(5) + qJD(6) - t158;
t1 = [t208, 0, 0, t195 ^ 2 * t208, t195 * t205, t200, t199, qJD(2) ^ 2 / 0.2e1, pkin(1) * t205 - pkin(7) * t200, -t197 * pkin(1) * t195 - pkin(7) * t199, -t174 * t182 + t175 * t181, t175 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t178 ^ 2 / 0.2e1, -t178 * t177, t178 * t180, -t177 * t180, t180 ^ 2 / 0.2e1, t172 * t177 + t198 * t180, t172 * t178 - t204 * t180, -t158 * t167 - t159 * t166, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, -t156 * t180 + t160 * t166, t156 * t167 - t157 * t166, t157 * t180 - t160 * t167, t157 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1;];
T_reg  = t1;
