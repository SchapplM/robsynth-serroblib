% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:53
% EndTime: 2019-03-09 07:05:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (547->55), mult. (1435->111), div. (0->0), fcn. (1152->10), ass. (0->47)
t203 = cos(qJ(3));
t202 = cos(qJ(4));
t201 = pkin(7) + qJ(2);
t188 = sin(qJ(3));
t184 = cos(pkin(11));
t196 = qJD(1) * t184;
t183 = sin(pkin(11));
t197 = qJD(1) * t183;
t171 = t188 * t197 - t203 * t196;
t172 = (t203 * t183 + t184 * t188) * qJD(1);
t187 = sin(qJ(4));
t165 = -t187 * t171 + t202 * t172;
t182 = qJD(3) + qJD(4);
t173 = t201 * t197;
t174 = t201 * t196;
t194 = -t203 * t173 - t188 * t174;
t161 = qJD(3) * pkin(3) - t172 * pkin(8) + t194;
t198 = -t188 * t173 + t203 * t174;
t162 = -t171 * pkin(8) + t198;
t195 = t202 * t161 - t187 * t162;
t149 = t182 * pkin(4) - t165 * pkin(9) + t195;
t164 = t202 * t171 + t187 * t172;
t199 = t187 * t161 + t202 * t162;
t151 = -t164 * pkin(9) + t199;
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t200 = t186 * t149 + t190 * t151;
t155 = t190 * t164 + t186 * t165;
t193 = t190 * t149 - t186 * t151;
t175 = qJD(2) + (-pkin(2) * t184 - pkin(1)) * qJD(1);
t166 = t171 * pkin(3) + t175;
t157 = t164 * pkin(4) + t166;
t191 = qJD(1) ^ 2;
t189 = cos(qJ(6));
t185 = sin(qJ(6));
t181 = t184 ^ 2;
t180 = t183 ^ 2;
t179 = qJD(5) + t182;
t178 = -qJD(1) * pkin(1) + qJD(2);
t156 = -t186 * t164 + t190 * t165;
t154 = qJD(6) + t155;
t153 = t189 * t156 + t185 * t179;
t152 = t185 * t156 - t189 * t179;
t147 = t155 * pkin(5) - t156 * pkin(10) + t157;
t146 = t179 * pkin(10) + t200;
t145 = -t179 * pkin(5) - t193;
t1 = [t191 / 0.2e1, 0, 0, -t178 * t196, t178 * t197 (t180 + t181) * t191 * qJ(2), t178 ^ 2 / 0.2e1 + (t181 / 0.2e1 + t180 / 0.2e1) * qJ(2) ^ 2 * t191, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * qJD(3), -t171 * qJD(3), qJD(3) ^ 2 / 0.2e1, t194 * qJD(3) + t175 * t171, -t198 * qJD(3) + t175 * t172, t165 ^ 2 / 0.2e1, -t165 * t164, t182 * t165, -t182 * t164, t182 ^ 2 / 0.2e1, t166 * t164 + t195 * t182, t166 * t165 - t199 * t182, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t179, -t155 * t179, t179 ^ 2 / 0.2e1, t157 * t155 + t193 * t179, t157 * t156 - t200 * t179, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t154, -t152 * t154, t154 ^ 2 / 0.2e1 (-t185 * t146 + t189 * t147) * t154 + t145 * t152 -(t189 * t146 + t185 * t147) * t154 + t145 * t153;];
T_reg  = t1;
