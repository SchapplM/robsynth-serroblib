% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:42
% EndTime: 2019-03-09 05:17:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (513->54), mult. (1262->108), div. (0->0), fcn. (962->10), ass. (0->46)
t211 = cos(qJ(4));
t210 = pkin(7) + qJ(2);
t194 = sin(pkin(10));
t196 = cos(pkin(10));
t199 = sin(qJ(3));
t201 = cos(qJ(3));
t184 = (t194 * t201 + t196 * t199) * qJD(1);
t198 = sin(qJ(4));
t177 = t198 * qJD(3) + t184 * t211;
t206 = qJD(1) * t196;
t207 = qJD(1) * t194;
t183 = t199 * t207 - t201 * t206;
t179 = qJD(4) + t183;
t187 = qJD(2) + (-pkin(2) * t196 - pkin(1)) * qJD(1);
t172 = t183 * pkin(3) - t184 * pkin(8) + t187;
t185 = t210 * t207;
t186 = t210 * t206;
t208 = -t199 * t185 + t201 * t186;
t175 = qJD(3) * pkin(8) + t208;
t205 = t172 * t211 - t198 * t175;
t160 = t179 * pkin(4) - t177 * qJ(5) + t205;
t176 = -qJD(3) * t211 + t198 * t184;
t209 = t198 * t172 + t175 * t211;
t165 = -t176 * qJ(5) + t209;
t193 = sin(pkin(11));
t195 = cos(pkin(11));
t157 = t193 * t160 + t195 * t165;
t156 = t195 * t160 - t193 * t165;
t204 = -t201 * t185 - t199 * t186;
t174 = -qJD(3) * pkin(3) - t204;
t166 = t176 * pkin(4) + qJD(5) + t174;
t202 = qJD(1) ^ 2;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t192 = t196 ^ 2;
t191 = t194 ^ 2;
t189 = -qJD(1) * pkin(1) + qJD(2);
t178 = qJD(6) + t179;
t169 = -t193 * t176 + t195 * t177;
t168 = -t195 * t176 - t193 * t177;
t163 = t197 * t168 + t200 * t169;
t162 = -t200 * t168 + t197 * t169;
t161 = -t168 * pkin(5) + t166;
t155 = t168 * pkin(9) + t157;
t154 = t179 * pkin(5) - t169 * pkin(9) + t156;
t1 = [t202 / 0.2e1, 0, 0, -t189 * t206, t189 * t207 (t191 + t192) * t202 * qJ(2), t189 ^ 2 / 0.2e1 + (t192 / 0.2e1 + t191 / 0.2e1) * qJ(2) ^ 2 * t202, t184 ^ 2 / 0.2e1, -t184 * t183, t184 * qJD(3), -t183 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t204 + t187 * t183, -qJD(3) * t208 + t187 * t184, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t179, -t176 * t179, t179 ^ 2 / 0.2e1, t174 * t176 + t179 * t205, t174 * t177 - t179 * t209, -t156 * t169 + t157 * t168, t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t178, -t162 * t178, t178 ^ 2 / 0.2e1 (t200 * t154 - t197 * t155) * t178 + t161 * t162 -(t197 * t154 + t200 * t155) * t178 + t161 * t163;];
T_reg  = t1;
