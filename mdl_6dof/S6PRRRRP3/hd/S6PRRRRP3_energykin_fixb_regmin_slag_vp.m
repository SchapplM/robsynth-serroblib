% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:06
% EndTime: 2019-03-09 00:11:06
% DurationCPUTime: 0.14s
% Computational Cost: add. (253->42), mult. (572->91), div. (0->0), fcn. (412->10), ass. (0->41)
t185 = qJD(2) ^ 2;
t200 = t185 / 0.2e1;
t199 = cos(qJ(5));
t179 = sin(qJ(4));
t182 = cos(qJ(4));
t180 = sin(qJ(3));
t193 = qJD(2) * t180;
t167 = t179 * qJD(3) + t182 * t193;
t183 = cos(qJ(3));
t192 = t183 * qJD(2);
t173 = -qJD(4) + t192;
t181 = sin(qJ(2));
t195 = qJD(1) * sin(pkin(6));
t168 = qJD(2) * pkin(8) + t181 * t195;
t194 = qJD(1) * cos(pkin(6));
t196 = t183 * t168 + t180 * t194;
t159 = qJD(3) * pkin(9) + t196;
t184 = cos(qJ(2));
t190 = t184 * t195;
t162 = -t190 + (-pkin(3) * t183 - pkin(9) * t180 - pkin(2)) * qJD(2);
t187 = -t179 * t159 + t182 * t162;
t150 = -t173 * pkin(4) - t167 * pkin(10) + t187;
t166 = -t182 * qJD(3) + t179 * t193;
t197 = t182 * t159 + t179 * t162;
t153 = -t166 * pkin(10) + t197;
t178 = sin(qJ(5));
t198 = t178 * t150 + t199 * t153;
t191 = qJD(2) * qJD(3);
t189 = qJD(2) * t195;
t188 = t199 * t150 - t178 * t153;
t186 = -t180 * t168 + t183 * t194;
t158 = -qJD(3) * pkin(3) - t186;
t154 = t166 * pkin(4) + t158;
t170 = -qJD(5) + t173;
t169 = -qJD(2) * pkin(2) - t190;
t156 = -t178 * t166 + t199 * t167;
t155 = t199 * t166 + t178 * t167;
t151 = t155 * pkin(5) + qJD(6) + t154;
t147 = -t155 * qJ(6) + t198;
t146 = -t170 * pkin(5) - t156 * qJ(6) + t188;
t1 = [qJD(1) ^ 2 / 0.2e1, t200, t184 * t189, -t181 * t189, t180 ^ 2 * t200, t180 * t185 * t183, t180 * t191, t183 * t191, qJD(3) ^ 2 / 0.2e1, t186 * qJD(3) - t169 * t192, -t196 * qJD(3) + t169 * t193, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t173, t166 * t173, t173 ^ 2 / 0.2e1, t158 * t166 - t187 * t173, t158 * t167 + t197 * t173, t156 ^ 2 / 0.2e1, -t156 * t155, -t156 * t170, t155 * t170, t170 ^ 2 / 0.2e1, t154 * t155 - t188 * t170, t154 * t156 + t198 * t170, -t146 * t156 - t147 * t155, t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1;];
T_reg  = t1;
