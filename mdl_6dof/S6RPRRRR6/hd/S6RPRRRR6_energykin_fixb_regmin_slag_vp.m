% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:32
% EndTime: 2019-03-09 07:14:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (499->55), mult. (1227->111), div. (0->0), fcn. (960->10), ass. (0->47)
t214 = cos(qJ(5));
t213 = pkin(7) + qJ(2);
t194 = sin(pkin(11));
t195 = cos(pkin(11));
t199 = sin(qJ(3));
t202 = cos(qJ(3));
t185 = (t194 * t202 + t195 * t199) * qJD(1);
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t177 = t198 * qJD(3) + t201 * t185;
t208 = qJD(1) * t195;
t209 = qJD(1) * t194;
t184 = t199 * t209 - t202 * t208;
t180 = qJD(4) + t184;
t188 = qJD(2) + (-pkin(2) * t195 - pkin(1)) * qJD(1);
t171 = t184 * pkin(3) - t185 * pkin(8) + t188;
t186 = t213 * t209;
t187 = t213 * t208;
t210 = -t199 * t186 + t202 * t187;
t174 = qJD(3) * pkin(8) + t210;
t206 = t201 * t171 - t198 * t174;
t159 = t180 * pkin(4) - t177 * pkin(9) + t206;
t176 = -t201 * qJD(3) + t198 * t185;
t211 = t198 * t171 + t201 * t174;
t164 = -t176 * pkin(9) + t211;
t197 = sin(qJ(5));
t212 = t197 * t159 + t214 * t164;
t207 = t214 * t159 - t197 * t164;
t205 = -t202 * t186 - t199 * t187;
t173 = -qJD(3) * pkin(3) - t205;
t179 = qJD(5) + t180;
t165 = t176 * pkin(4) + t173;
t203 = qJD(1) ^ 2;
t200 = cos(qJ(6));
t196 = sin(qJ(6));
t193 = t195 ^ 2;
t192 = t194 ^ 2;
t190 = -qJD(1) * pkin(1) + qJD(2);
t178 = qJD(6) + t179;
t168 = -t197 * t176 + t177 * t214;
t167 = t176 * t214 + t197 * t177;
t162 = t167 * pkin(5) + t165;
t161 = -t196 * t167 + t200 * t168;
t160 = t200 * t167 + t196 * t168;
t156 = -t167 * pkin(10) + t212;
t155 = t179 * pkin(5) - t168 * pkin(10) + t207;
t1 = [t203 / 0.2e1, 0, 0, -t190 * t208, t190 * t209 (t192 + t193) * t203 * qJ(2), t190 ^ 2 / 0.2e1 + (t193 / 0.2e1 + t192 / 0.2e1) * qJ(2) ^ 2 * t203, t185 ^ 2 / 0.2e1, -t185 * t184, t185 * qJD(3), -t184 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t205 + t188 * t184, -qJD(3) * t210 + t188 * t185, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t180, -t176 * t180, t180 ^ 2 / 0.2e1, t173 * t176 + t180 * t206, t173 * t177 - t180 * t211, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t179, -t167 * t179, t179 ^ 2 / 0.2e1, t165 * t167 + t179 * t207, t165 * t168 - t179 * t212, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t178, -t160 * t178, t178 ^ 2 / 0.2e1 (t200 * t155 - t196 * t156) * t178 + t162 * t160 -(t196 * t155 + t200 * t156) * t178 + t162 * t161;];
T_reg  = t1;
