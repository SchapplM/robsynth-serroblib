% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:34
% EndTime: 2019-03-09 03:49:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (313->52), mult. (789->103), div. (0->0), fcn. (559->8), ass. (0->41)
t199 = pkin(7) + qJ(2);
t181 = sin(pkin(10));
t182 = cos(pkin(10));
t186 = sin(qJ(3));
t189 = cos(qJ(3));
t168 = (t181 * t189 + t182 * t186) * qJD(1);
t196 = qJD(1) * t181;
t169 = t199 * t196;
t195 = qJD(1) * t182;
t170 = t199 * t195;
t194 = -t189 * t169 - t186 * t170;
t193 = qJD(4) - t194;
t150 = -t168 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(3) + t193;
t197 = -t186 * t169 + t189 * t170;
t160 = qJD(3) * qJ(4) + t197;
t167 = t186 * t196 - t189 * t195;
t152 = t167 * pkin(8) + t160;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t198 = t185 * t150 + t188 * t152;
t175 = -qJD(1) * pkin(1) + qJD(2);
t171 = -pkin(2) * t195 + t175;
t157 = -t188 * t167 + t185 * t168;
t155 = t167 * pkin(3) - t168 * qJ(4) + t171;
t192 = t188 * t150 - t185 * t152;
t148 = -t167 * pkin(4) - t155;
t190 = qJD(1) ^ 2;
t187 = cos(qJ(6));
t184 = sin(qJ(6));
t178 = qJD(3) - qJD(5);
t177 = t182 ^ 2;
t176 = t181 ^ 2;
t159 = -qJD(3) * pkin(3) + t193;
t158 = t185 * t167 + t188 * t168;
t156 = qJD(6) + t157;
t154 = t187 * t158 - t184 * t178;
t153 = t184 * t158 + t187 * t178;
t147 = -t178 * pkin(9) + t198;
t146 = t178 * pkin(5) - t192;
t145 = t157 * pkin(5) - t158 * pkin(9) + t148;
t1 = [t190 / 0.2e1, 0, 0, -t175 * t195, t175 * t196 (t176 + t177) * t190 * qJ(2), t175 ^ 2 / 0.2e1 + (t177 / 0.2e1 + t176 / 0.2e1) * qJ(2) ^ 2 * t190, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * qJD(3), -t167 * qJD(3), qJD(3) ^ 2 / 0.2e1, t194 * qJD(3) + t171 * t167, -t197 * qJD(3) + t171 * t168, -t159 * qJD(3) + t155 * t167, t159 * t168 - t160 * t167, t160 * qJD(3) - t155 * t168, t160 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t158 ^ 2 / 0.2e1, -t158 * t157, -t158 * t178, t157 * t178, t178 ^ 2 / 0.2e1, t148 * t157 - t192 * t178, t148 * t158 + t198 * t178, t154 ^ 2 / 0.2e1, -t154 * t153, t154 * t156, -t153 * t156, t156 ^ 2 / 0.2e1 (t187 * t145 - t184 * t147) * t156 + t146 * t153 -(t184 * t145 + t187 * t147) * t156 + t146 * t154;];
T_reg  = t1;
