% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:40
% EndTime: 2019-03-10 01:15:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (613->50), mult. (1265->105), div. (0->0), fcn. (925->8), ass. (0->45)
t198 = -pkin(8) - pkin(7);
t183 = qJD(1) ^ 2;
t197 = t183 / 0.2e1;
t196 = cos(qJ(4));
t182 = cos(qJ(2));
t195 = t182 * t183;
t178 = sin(qJ(3));
t179 = sin(qJ(2));
t181 = cos(qJ(3));
t167 = (t178 * t182 + t179 * t181) * qJD(1);
t175 = qJD(2) + qJD(3);
t177 = sin(qJ(4));
t162 = t196 * t167 + t177 * t175;
t190 = qJD(1) * t182;
t191 = qJD(1) * t179;
t166 = t178 * t191 - t181 * t190;
t165 = qJD(4) + t166;
t172 = (-pkin(2) * t182 - pkin(1)) * qJD(1);
t156 = t166 * pkin(3) - t167 * pkin(9) + t172;
t170 = qJD(2) * pkin(2) + t198 * t191;
t171 = t198 * t190;
t192 = t178 * t170 - t181 * t171;
t159 = t175 * pkin(9) + t192;
t186 = t196 * t156 - t177 * t159;
t148 = t165 * pkin(4) - t162 * pkin(10) + t186;
t161 = t177 * t167 - t196 * t175;
t193 = t177 * t156 + t196 * t159;
t150 = -t161 * pkin(10) + t193;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t194 = t176 * t148 + t180 * t150;
t189 = qJD(1) * qJD(2);
t188 = t179 * t189;
t187 = t182 * t189;
t185 = t181 * t170 + t178 * t171;
t184 = t180 * t148 - t176 * t150;
t158 = -t175 * pkin(3) - t185;
t151 = t161 * pkin(4) + t158;
t163 = qJD(5) + t165;
t153 = -t176 * t161 + t180 * t162;
t152 = t180 * t161 + t176 * t162;
t146 = t152 * pkin(5) - t153 * qJ(6) + t151;
t145 = t163 * qJ(6) + t194;
t144 = -t163 * pkin(5) + qJD(6) - t184;
t1 = [t197, 0, 0, t179 ^ 2 * t197, t179 * t195, t188, t187, qJD(2) ^ 2 / 0.2e1, pkin(1) * t195 - pkin(7) * t188, -t183 * pkin(1) * t179 - pkin(7) * t187, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t175, -t166 * t175, t175 ^ 2 / 0.2e1, t172 * t166 + t185 * t175, t172 * t167 - t192 * t175, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t165, -t161 * t165, t165 ^ 2 / 0.2e1, t158 * t161 + t186 * t165, t158 * t162 - t193 * t165, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t163, -t152 * t163, t163 ^ 2 / 0.2e1, t151 * t152 + t184 * t163, t151 * t153 - t194 * t163, -t144 * t163 + t146 * t152, t144 * t153 - t145 * t152, t145 * t163 - t146 * t153, t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
