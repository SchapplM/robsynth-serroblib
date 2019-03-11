% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:22
% EndTime: 2019-03-09 02:48:22
% DurationCPUTime: 0.13s
% Computational Cost: add. (398->54), mult. (993->105), div. (0->0), fcn. (698->8), ass. (0->42)
t204 = -pkin(4) - pkin(5);
t203 = cos(qJ(3));
t202 = pkin(7) + qJ(2);
t191 = sin(qJ(3));
t189 = cos(pkin(9));
t198 = qJD(1) * t189;
t187 = sin(pkin(9));
t199 = qJD(1) * t187;
t176 = t191 * t199 - t203 * t198;
t177 = (t203 * t187 + t189 * t191) * qJD(1);
t180 = qJD(2) + (-pkin(2) * t189 - pkin(1)) * qJD(1);
t161 = t176 * pkin(3) - t177 * qJ(4) + t180;
t178 = t202 * t199;
t179 = t202 * t198;
t200 = -t191 * t178 + t203 * t179;
t165 = qJD(3) * qJ(4) + t200;
t186 = sin(pkin(10));
t188 = cos(pkin(10));
t157 = t186 * t161 + t188 * t165;
t201 = -t203 * t178 - t191 * t179;
t154 = t176 * qJ(5) + t157;
t156 = t188 * t161 - t186 * t165;
t197 = qJD(5) - t156;
t196 = qJD(3) * pkin(3) - qJD(4) + t201;
t168 = t186 * qJD(3) + t188 * t177;
t195 = t168 * qJ(5) + t196;
t193 = qJD(1) ^ 2;
t192 = cos(qJ(6));
t190 = sin(qJ(6));
t185 = t189 ^ 2;
t184 = t187 ^ 2;
t182 = -qJD(1) * pkin(1) + qJD(2);
t170 = -qJD(6) + t176;
t167 = -t188 * qJD(3) + t186 * t177;
t159 = t190 * t167 + t192 * t168;
t158 = -t192 * t167 + t190 * t168;
t155 = t167 * pkin(4) - t195;
t153 = -t176 * pkin(4) + t197;
t152 = t204 * t167 + t195;
t151 = t167 * pkin(8) + t154;
t150 = -t168 * pkin(8) + t204 * t176 + t197;
t1 = [t193 / 0.2e1, 0, 0, -t182 * t198, t182 * t199 (t184 + t185) * t193 * qJ(2), t182 ^ 2 / 0.2e1 + (t185 / 0.2e1 + t184 / 0.2e1) * qJ(2) ^ 2 * t193, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * qJD(3), -t176 * qJD(3), qJD(3) ^ 2 / 0.2e1, t201 * qJD(3) + t180 * t176, -t200 * qJD(3) + t180 * t177, t156 * t176 - t167 * t196, -t157 * t176 - t168 * t196, -t156 * t168 - t157 * t167, t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t196 ^ 2 / 0.2e1, -t153 * t176 + t155 * t167, t153 * t168 - t154 * t167, t154 * t176 - t155 * t168, t154 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1, t159 ^ 2 / 0.2e1, -t159 * t158, -t159 * t170, t158 * t170, t170 ^ 2 / 0.2e1 -(t192 * t150 - t190 * t151) * t170 + t152 * t158 (t190 * t150 + t192 * t151) * t170 + t152 * t159;];
T_reg  = t1;
