% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:45
% EndTime: 2019-03-09 04:40:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (533->51), mult. (1277->101), div. (0->0), fcn. (924->8), ass. (0->41)
t200 = cos(qJ(4));
t199 = pkin(7) + qJ(2);
t185 = sin(pkin(9));
t187 = cos(pkin(9));
t189 = sin(qJ(3));
t190 = cos(qJ(3));
t175 = (t185 * t190 + t187 * t189) * qJD(1);
t188 = sin(qJ(4));
t169 = t188 * qJD(3) + t200 * t175;
t195 = qJD(1) * t187;
t196 = qJD(1) * t185;
t174 = t189 * t196 - t190 * t195;
t170 = qJD(4) + t174;
t178 = qJD(2) + (-pkin(2) * t187 - pkin(1)) * qJD(1);
t163 = t174 * pkin(3) - t175 * pkin(8) + t178;
t176 = t199 * t196;
t177 = t199 * t195;
t197 = -t189 * t176 + t190 * t177;
t166 = qJD(3) * pkin(8) + t197;
t194 = t200 * t163 - t188 * t166;
t155 = t170 * pkin(4) - t169 * qJ(5) + t194;
t168 = -t200 * qJD(3) + t188 * t175;
t198 = t188 * t163 + t200 * t166;
t157 = -t168 * qJ(5) + t198;
t184 = sin(pkin(10));
t186 = cos(pkin(10));
t152 = t184 * t155 + t186 * t157;
t193 = -t190 * t176 - t189 * t177;
t151 = t186 * t155 - t184 * t157;
t165 = -qJD(3) * pkin(3) - t193;
t158 = t168 * pkin(4) + qJD(5) + t165;
t191 = qJD(1) ^ 2;
t183 = t187 ^ 2;
t182 = t185 ^ 2;
t180 = -qJD(1) * pkin(1) + qJD(2);
t160 = -t184 * t168 + t186 * t169;
t159 = t186 * t168 + t184 * t169;
t153 = t159 * pkin(5) - t160 * qJ(6) + t158;
t150 = t170 * qJ(6) + t152;
t149 = -t170 * pkin(5) + qJD(6) - t151;
t1 = [t191 / 0.2e1, 0, 0, -t180 * t195, t180 * t196 (t182 + t183) * t191 * qJ(2), t180 ^ 2 / 0.2e1 + (t183 / 0.2e1 + t182 / 0.2e1) * qJ(2) ^ 2 * t191, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * qJD(3), -t174 * qJD(3), qJD(3) ^ 2 / 0.2e1, t193 * qJD(3) + t178 * t174, -t197 * qJD(3) + t178 * t175, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t170, -t168 * t170, t170 ^ 2 / 0.2e1, t165 * t168 + t194 * t170, t165 * t169 - t198 * t170, -t151 * t160 - t152 * t159, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, -t149 * t170 + t153 * t159, t149 * t160 - t150 * t159, t150 * t170 - t153 * t160, t150 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1;];
T_reg  = t1;
