% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:46:47
% EndTime: 2019-03-09 20:46:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (622->49), mult. (1303->102), div. (0->0), fcn. (927->8), ass. (0->44)
t202 = -pkin(8) - pkin(7);
t189 = qJD(1) ^ 2;
t201 = t189 / 0.2e1;
t200 = cos(qJ(4));
t188 = cos(qJ(2));
t199 = t188 * t189;
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t187 = cos(qJ(3));
t173 = (t185 * t188 + t186 * t187) * qJD(1);
t181 = qJD(2) + qJD(3);
t184 = sin(qJ(4));
t169 = t200 * t173 + t184 * t181;
t195 = qJD(1) * t188;
t196 = qJD(1) * t186;
t172 = t185 * t196 - t187 * t195;
t171 = qJD(4) + t172;
t178 = (-pkin(2) * t188 - pkin(1)) * qJD(1);
t163 = t172 * pkin(3) - t173 * pkin(9) + t178;
t176 = qJD(2) * pkin(2) + t202 * t196;
t177 = t202 * t195;
t197 = t185 * t176 - t187 * t177;
t166 = t181 * pkin(9) + t197;
t191 = t200 * t163 - t184 * t166;
t155 = t171 * pkin(4) - t169 * qJ(5) + t191;
t168 = t184 * t173 - t200 * t181;
t198 = t184 * t163 + t200 * t166;
t157 = -t168 * qJ(5) + t198;
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t152 = t182 * t155 + t183 * t157;
t194 = qJD(1) * qJD(2);
t193 = t186 * t194;
t192 = t188 * t194;
t190 = t187 * t176 + t185 * t177;
t151 = t183 * t155 - t182 * t157;
t165 = -t181 * pkin(3) - t190;
t158 = t168 * pkin(4) + qJD(5) + t165;
t160 = -t182 * t168 + t183 * t169;
t159 = t183 * t168 + t182 * t169;
t153 = t159 * pkin(5) - t160 * qJ(6) + t158;
t150 = t171 * qJ(6) + t152;
t149 = -t171 * pkin(5) + qJD(6) - t151;
t1 = [t201, 0, 0, t186 ^ 2 * t201, t186 * t199, t193, t192, qJD(2) ^ 2 / 0.2e1, pkin(1) * t199 - pkin(7) * t193, -t189 * pkin(1) * t186 - pkin(7) * t192, t173 ^ 2 / 0.2e1, -t173 * t172, t173 * t181, -t172 * t181, t181 ^ 2 / 0.2e1, t178 * t172 + t190 * t181, t178 * t173 - t197 * t181, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t171, -t168 * t171, t171 ^ 2 / 0.2e1, t165 * t168 + t191 * t171, t165 * t169 - t198 * t171, -t151 * t160 - t152 * t159, t152 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, -t149 * t171 + t153 * t159, t149 * t160 - t150 * t159, t150 * t171 - t153 * t160, t150 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1;];
T_reg  = t1;
