% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:17
% EndTime: 2019-03-09 16:37:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (583->49), mult. (1363->102), div. (0->0), fcn. (1000->8), ass. (0->44)
t202 = -pkin(8) - pkin(7);
t189 = qJD(1) ^ 2;
t201 = t189 / 0.2e1;
t200 = cos(qJ(3));
t188 = cos(qJ(2));
t199 = t188 * t189;
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t174 = (t185 * t188 + t200 * t186) * qJD(1);
t181 = qJD(2) + qJD(3);
t196 = qJD(1) * t186;
t176 = qJD(2) * pkin(2) + t202 * t196;
t195 = qJD(1) * t188;
t177 = t202 * t195;
t191 = t200 * t176 + t185 * t177;
t160 = t181 * pkin(3) - t174 * qJ(4) + t191;
t173 = t185 * t196 - t200 * t195;
t197 = t185 * t176 - t200 * t177;
t163 = -t173 * qJ(4) + t197;
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t156 = t182 * t160 + t183 * t163;
t154 = t181 * pkin(9) + t156;
t167 = -t183 * t173 - t182 * t174;
t168 = -t182 * t173 + t183 * t174;
t178 = (-pkin(2) * t188 - pkin(1)) * qJD(1);
t169 = t173 * pkin(3) + qJD(4) + t178;
t158 = -t167 * pkin(4) - t168 * pkin(9) + t169;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t198 = t187 * t154 + t184 * t158;
t194 = qJD(1) * qJD(2);
t193 = t186 * t194;
t192 = t188 * t194;
t155 = t183 * t160 - t182 * t163;
t190 = -t184 * t154 + t187 * t158;
t153 = -t181 * pkin(4) - t155;
t166 = qJD(5) - t167;
t165 = t187 * t168 + t184 * t181;
t164 = t184 * t168 - t187 * t181;
t151 = t164 * pkin(5) - t165 * qJ(6) + t153;
t150 = t166 * qJ(6) + t198;
t149 = -t166 * pkin(5) + qJD(6) - t190;
t1 = [t201, 0, 0, t186 ^ 2 * t201, t186 * t199, t193, t192, qJD(2) ^ 2 / 0.2e1, pkin(1) * t199 - pkin(7) * t193, -t189 * pkin(1) * t186 - pkin(7) * t192, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t181, -t173 * t181, t181 ^ 2 / 0.2e1, t178 * t173 + t191 * t181, t178 * t174 - t197 * t181, -t155 * t168 + t156 * t167, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t166, -t164 * t166, t166 ^ 2 / 0.2e1, t153 * t164 + t190 * t166, t153 * t165 - t198 * t166, -t149 * t166 + t151 * t164, t149 * t165 - t150 * t164, t150 * t166 - t151 * t165, t150 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1;];
T_reg  = t1;
