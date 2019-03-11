% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:34
% EndTime: 2019-03-08 21:26:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (233->42), mult. (564->90), div. (0->0), fcn. (408->10), ass. (0->40)
t190 = qJD(2) ^ 2;
t202 = t190 / 0.2e1;
t201 = cos(qJ(5));
t187 = sin(qJ(2));
t198 = qJD(1) * sin(pkin(6));
t175 = qJD(2) * pkin(8) + t187 * t198;
t188 = cos(qJ(3));
t197 = qJD(1) * cos(pkin(6));
t179 = t188 * t197;
t186 = sin(qJ(3));
t166 = qJD(3) * pkin(3) + t179 + (-qJ(4) * qJD(2) - t175) * t186;
t195 = qJD(2) * t188;
t199 = t188 * t175 + t186 * t197;
t167 = qJ(4) * t195 + t199;
t181 = sin(pkin(11));
t183 = cos(pkin(11));
t159 = t181 * t166 + t183 * t167;
t157 = qJD(3) * pkin(9) + t159;
t189 = cos(qJ(2));
t193 = t189 * t198;
t170 = -t193 + qJD(4) + (-pkin(3) * t188 - pkin(2)) * qJD(2);
t196 = qJD(2) * t186;
t172 = -t181 * t196 + t183 * t195;
t173 = (t181 * t188 + t183 * t186) * qJD(2);
t162 = -t172 * pkin(4) - t173 * pkin(9) + t170;
t185 = sin(qJ(5));
t200 = t201 * t157 + t185 * t162;
t194 = qJD(2) * qJD(3);
t192 = qJD(2) * t198;
t191 = -t185 * t157 + t201 * t162;
t158 = t183 * t166 - t181 * t167;
t156 = -qJD(3) * pkin(4) - t158;
t176 = -qJD(2) * pkin(2) - t193;
t171 = qJD(5) - t172;
t169 = t185 * qJD(3) + t201 * t173;
t168 = -t201 * qJD(3) + t185 * t173;
t154 = t168 * pkin(5) + qJD(6) + t156;
t153 = -qJ(6) * t168 + t200;
t152 = t171 * pkin(5) - t169 * qJ(6) + t191;
t1 = [qJD(1) ^ 2 / 0.2e1, t202, t189 * t192, -t187 * t192, t186 ^ 2 * t202, t186 * t190 * t188, t186 * t194, t188 * t194, qJD(3) ^ 2 / 0.2e1 (-t186 * t175 + t179) * qJD(3) - t176 * t195, -t199 * qJD(3) + t176 * t196, -t158 * t173 + t159 * t172, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t171, -t168 * t171, t171 ^ 2 / 0.2e1, t156 * t168 + t191 * t171, t156 * t169 - t200 * t171, -t152 * t169 - t153 * t168, t153 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1;];
T_reg  = t1;
