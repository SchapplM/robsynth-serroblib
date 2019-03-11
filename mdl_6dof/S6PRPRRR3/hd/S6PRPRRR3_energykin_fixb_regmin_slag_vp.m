% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:11
% EndTime: 2019-03-08 20:34:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (299->49), mult. (754->101), div. (0->0), fcn. (600->12), ass. (0->43)
t206 = cos(qJ(4));
t185 = sin(pkin(12));
t187 = cos(pkin(12));
t191 = sin(qJ(4));
t176 = (t206 * t185 + t187 * t191) * qJD(2);
t192 = sin(qJ(2));
t203 = qJD(1) * sin(pkin(6));
t179 = qJD(2) * qJ(3) + t192 * t203;
t202 = qJD(1) * cos(pkin(6));
t181 = t187 * t202;
t169 = t181 + (-pkin(8) * qJD(2) - t179) * t185;
t172 = t187 * t179 + t185 * t202;
t200 = qJD(2) * t187;
t170 = pkin(8) * t200 + t172;
t198 = t206 * t169 - t191 * t170;
t158 = qJD(4) * pkin(4) - t176 * pkin(9) + t198;
t201 = qJD(2) * t185;
t175 = t191 * t201 - t206 * t200;
t204 = t191 * t169 + t206 * t170;
t159 = -t175 * pkin(9) + t204;
t190 = sin(qJ(5));
t194 = cos(qJ(5));
t205 = t190 * t158 + t194 * t159;
t199 = qJD(2) * t203;
t163 = t194 * t175 + t190 * t176;
t195 = cos(qJ(2));
t197 = -t195 * t203 + qJD(3);
t196 = t194 * t158 - t190 * t159;
t174 = (-pkin(3) * t187 - pkin(2)) * qJD(2) + t197;
t165 = t175 * pkin(4) + t174;
t193 = cos(qJ(6));
t189 = sin(qJ(6));
t184 = qJD(4) + qJD(5);
t178 = -qJD(2) * pkin(2) + t197;
t171 = -t185 * t179 + t181;
t164 = -t190 * t175 + t194 * t176;
t162 = qJD(6) + t163;
t161 = t193 * t164 + t189 * t184;
t160 = t189 * t164 - t193 * t184;
t155 = t163 * pkin(5) - t164 * pkin(10) + t165;
t154 = t184 * pkin(10) + t205;
t153 = -t184 * pkin(5) - t196;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t195 * t199, -t192 * t199, -t178 * t200, t178 * t201 (-t171 * t185 + t172 * t187) * qJD(2), t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * qJD(4), -t175 * qJD(4), qJD(4) ^ 2 / 0.2e1, t198 * qJD(4) + t174 * t175, -t204 * qJD(4) + t174 * t176, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t184, -t163 * t184, t184 ^ 2 / 0.2e1, t165 * t163 + t196 * t184, t165 * t164 - t205 * t184, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t162, -t160 * t162, t162 ^ 2 / 0.2e1 (-t189 * t154 + t193 * t155) * t162 + t153 * t160 -(t193 * t154 + t189 * t155) * t162 + t153 * t161;];
T_reg  = t1;
