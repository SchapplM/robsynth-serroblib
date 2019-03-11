% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:06
% EndTime: 2019-03-09 09:32:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (293->55), mult. (722->109), div. (0->0), fcn. (488->8), ass. (0->42)
t187 = sin(qJ(2));
t205 = qJ(3) * t187;
t183 = sin(pkin(6));
t191 = qJD(1) ^ 2;
t204 = t183 ^ 2 * t191;
t199 = cos(pkin(6)) * qJD(1);
t181 = qJD(2) + t199;
t190 = cos(qJ(2));
t200 = qJD(1) * t183;
t195 = t190 * t200;
t198 = pkin(1) * t199;
t201 = pkin(8) * t195 + t187 * t198;
t165 = -t181 * qJ(3) - t201;
t162 = pkin(3) * t195 + qJD(4) - t165;
t156 = pkin(4) * t195 - t181 * pkin(9) + t162;
t194 = -pkin(1) + (-pkin(2) - qJ(4)) * t190;
t159 = ((pkin(9) - qJ(3)) * t187 + t194) * t200;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t203 = t186 * t156 + t189 * t159;
t196 = t187 * t200;
t202 = -pkin(8) * t196 + t190 * t198;
t197 = t190 * t204;
t168 = t186 * t181 - t189 * t196;
t164 = -t181 * pkin(2) + qJD(3) - t202;
t193 = t181 * qJ(4) - t164;
t192 = t189 * t156 - t186 * t159;
t155 = (-pkin(3) - pkin(4)) * t196 + t193;
t188 = cos(qJ(6));
t185 = sin(qJ(6));
t172 = qJD(5) + t195;
t169 = t189 * t181 + t186 * t196;
t167 = qJD(6) + t168;
t166 = (-pkin(2) * t190 - pkin(1) - t205) * t200;
t163 = (t194 - t205) * t200;
t161 = t188 * t169 + t185 * t172;
t160 = t185 * t169 - t188 * t172;
t157 = pkin(3) * t196 - t193;
t153 = t168 * pkin(5) - t169 * pkin(10) + t155;
t152 = t172 * pkin(10) + t203;
t151 = -t172 * pkin(5) - t192;
t1 = [t191 / 0.2e1, 0, 0, t187 ^ 2 * t204 / 0.2e1, t187 * t197, t181 * t196, t181 * t195, t181 ^ 2 / 0.2e1, pkin(1) * t197 + t202 * t181, -pkin(1) * t187 * t204 - t201 * t181 (t164 * t187 - t165 * t190) * t200, t164 * t181 + t166 * t195, -t165 * t181 - t166 * t196, t166 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 (t157 * t187 + t162 * t190) * t200, t162 * t181 - t163 * t196, -t157 * t181 - t163 * t195, t163 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t172, -t168 * t172, t172 ^ 2 / 0.2e1, t155 * t168 + t192 * t172, t155 * t169 - t203 * t172, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t167, -t160 * t167, t167 ^ 2 / 0.2e1 (-t185 * t152 + t188 * t153) * t167 + t151 * t160 -(t188 * t152 + t185 * t153) * t167 + t151 * t161;];
T_reg  = t1;
