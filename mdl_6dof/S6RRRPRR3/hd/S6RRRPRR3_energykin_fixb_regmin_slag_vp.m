% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:22
% EndTime: 2019-03-09 18:13:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (366->51), mult. (805->105), div. (0->0), fcn. (562->8), ass. (0->44)
t202 = -pkin(8) - pkin(7);
t189 = qJD(1) ^ 2;
t201 = t189 / 0.2e1;
t188 = cos(qJ(2));
t200 = t188 * t189;
t183 = sin(qJ(3));
t184 = sin(qJ(2));
t187 = cos(qJ(3));
t167 = (t183 * t188 + t184 * t187) * qJD(1);
t179 = qJD(2) + qJD(3);
t197 = qJD(1) * t184;
t170 = qJD(2) * pkin(2) + t202 * t197;
t196 = qJD(1) * t188;
t171 = t202 * t196;
t192 = t187 * t170 + t183 * t171;
t191 = qJD(4) - t192;
t150 = -t167 * pkin(9) + (-pkin(3) - pkin(4)) * t179 + t191;
t198 = t183 * t170 - t187 * t171;
t159 = t179 * qJ(4) + t198;
t166 = t183 * t197 - t187 * t196;
t153 = t166 * pkin(9) + t159;
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t199 = t182 * t150 + t186 * t153;
t172 = -qJD(1) * pkin(1) - pkin(2) * t196;
t195 = qJD(1) * qJD(2);
t194 = t184 * t195;
t193 = t188 * t195;
t160 = -t186 * t166 + t182 * t167;
t156 = t166 * pkin(3) - t167 * qJ(4) + t172;
t190 = t186 * t150 - t182 * t153;
t151 = -t166 * pkin(4) - t156;
t185 = cos(qJ(6));
t181 = sin(qJ(6));
t177 = -qJD(5) + t179;
t161 = t182 * t166 + t186 * t167;
t158 = qJD(6) + t160;
t157 = -t179 * pkin(3) + t191;
t155 = t185 * t161 - t181 * t177;
t154 = t181 * t161 + t185 * t177;
t148 = -t177 * pkin(10) + t199;
t147 = t177 * pkin(5) - t190;
t146 = t160 * pkin(5) - t161 * pkin(10) + t151;
t1 = [t201, 0, 0, t184 ^ 2 * t201, t184 * t200, t194, t193, qJD(2) ^ 2 / 0.2e1, pkin(1) * t200 - pkin(7) * t194, -t189 * pkin(1) * t184 - pkin(7) * t193, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t179, -t166 * t179, t179 ^ 2 / 0.2e1, t172 * t166 + t192 * t179, t172 * t167 - t198 * t179, t156 * t166 - t157 * t179, t157 * t167 - t159 * t166, -t156 * t167 + t159 * t179, t159 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, -t161 * t177, t160 * t177, t177 ^ 2 / 0.2e1, t151 * t160 - t190 * t177, t151 * t161 + t199 * t177, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t158, -t154 * t158, t158 ^ 2 / 0.2e1 (t185 * t146 - t181 * t148) * t158 + t147 * t154 -(t181 * t146 + t185 * t148) * t158 + t147 * t155;];
T_reg  = t1;
