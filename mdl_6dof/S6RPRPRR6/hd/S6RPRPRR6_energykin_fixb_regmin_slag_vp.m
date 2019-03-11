% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:52
% EndTime: 2019-03-09 03:52:52
% DurationCPUTime: 0.16s
% Computational Cost: add. (525->56), mult. (1313->112), div. (0->0), fcn. (1004->10), ass. (0->46)
t211 = cos(qJ(5));
t210 = pkin(7) + qJ(2);
t199 = sin(qJ(3));
t201 = cos(qJ(3));
t196 = cos(pkin(10));
t206 = qJD(1) * t196;
t194 = sin(pkin(10));
t207 = qJD(1) * t194;
t183 = t199 * t207 - t201 * t206;
t184 = (t194 * t201 + t196 * t199) * qJD(1);
t187 = qJD(2) + (-pkin(2) * t196 - pkin(1)) * qJD(1);
t171 = t183 * pkin(3) - t184 * qJ(4) + t187;
t185 = t210 * t207;
t186 = t210 * t206;
t208 = -t185 * t199 + t186 * t201;
t174 = qJD(3) * qJ(4) + t208;
t193 = sin(pkin(11));
t195 = cos(pkin(11));
t163 = t171 * t195 - t174 * t193;
t177 = qJD(3) * t193 + t184 * t195;
t157 = pkin(4) * t183 - pkin(8) * t177 + t163;
t164 = t171 * t193 + t174 * t195;
t176 = -qJD(3) * t195 + t184 * t193;
t162 = -pkin(8) * t176 + t164;
t198 = sin(qJ(5));
t209 = t157 * t198 + t162 * t211;
t205 = t157 * t211 - t162 * t198;
t204 = -t185 * t201 - t186 * t199;
t179 = qJD(5) + t183;
t173 = -qJD(3) * pkin(3) + qJD(4) - t204;
t165 = pkin(4) * t176 + t173;
t202 = qJD(1) ^ 2;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t192 = t196 ^ 2;
t191 = t194 ^ 2;
t189 = -qJD(1) * pkin(1) + qJD(2);
t178 = qJD(6) + t179;
t168 = -t176 * t198 + t177 * t211;
t167 = t176 * t211 + t177 * t198;
t160 = -t167 * t197 + t168 * t200;
t159 = t167 * t200 + t168 * t197;
t158 = pkin(5) * t167 + t165;
t154 = -pkin(9) * t167 + t209;
t153 = pkin(5) * t179 - pkin(9) * t168 + t205;
t1 = [t202 / 0.2e1, 0, 0, -t189 * t206, t189 * t207 (t191 + t192) * t202 * qJ(2), t189 ^ 2 / 0.2e1 + (t192 / 0.2e1 + t191 / 0.2e1) * qJ(2) ^ 2 * t202, t184 ^ 2 / 0.2e1, -t184 * t183, t184 * qJD(3), -t183 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t204 + t183 * t187, -qJD(3) * t208 + t184 * t187, t163 * t183 + t173 * t176, -t164 * t183 + t173 * t177, -t163 * t177 - t164 * t176, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t179, -t167 * t179, t179 ^ 2 / 0.2e1, t165 * t167 + t179 * t205, t165 * t168 - t179 * t209, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t178, -t159 * t178, t178 ^ 2 / 0.2e1 (t153 * t200 - t154 * t197) * t178 + t158 * t159 -(t153 * t197 + t154 * t200) * t178 + t158 * t160;];
T_reg  = t1;
