% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:28
% EndTime: 2019-03-09 13:14:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (574->52), mult. (1457->109), div. (0->0), fcn. (1146->10), ass. (0->47)
t210 = qJD(1) * (pkin(7) + qJ(3));
t198 = qJD(1) ^ 2;
t209 = t198 / 0.2e1;
t208 = cos(qJ(4));
t197 = cos(qJ(2));
t206 = t197 * t198;
t189 = sin(pkin(11));
t190 = cos(pkin(11));
t194 = sin(qJ(2));
t180 = (-t189 * t194 + t190 * t197) * qJD(1);
t181 = (t189 * t197 + t190 * t194) * qJD(1);
t193 = sin(qJ(4));
t173 = t193 * t180 + t208 * t181;
t188 = qJD(2) + qJD(4);
t183 = qJD(2) * pkin(2) - t194 * t210;
t184 = t197 * t210;
t174 = t190 * t183 - t189 * t184;
t169 = qJD(2) * pkin(3) - t181 * pkin(8) + t174;
t175 = t189 * t183 + t190 * t184;
t170 = t180 * pkin(8) + t175;
t200 = t208 * t169 - t193 * t170;
t157 = t188 * pkin(4) - t173 * pkin(9) + t200;
t172 = -t208 * t180 + t193 * t181;
t204 = t193 * t169 + t208 * t170;
t159 = -t172 * pkin(9) + t204;
t192 = sin(qJ(5));
t196 = cos(qJ(5));
t205 = t192 * t157 + t196 * t159;
t203 = qJD(1) * qJD(2);
t202 = t194 * t203;
t201 = t197 * t203;
t163 = t196 * t172 + t192 * t173;
t199 = t196 * t157 - t192 * t159;
t185 = qJD(3) + (-pkin(2) * t197 - pkin(1)) * qJD(1);
t176 = -t180 * pkin(3) + t185;
t165 = t172 * pkin(4) + t176;
t195 = cos(qJ(6));
t191 = sin(qJ(6));
t187 = qJD(5) + t188;
t164 = -t192 * t172 + t196 * t173;
t162 = qJD(6) + t163;
t161 = t195 * t164 + t191 * t187;
t160 = t191 * t164 - t195 * t187;
t155 = t163 * pkin(5) - t164 * pkin(10) + t165;
t154 = t187 * pkin(10) + t205;
t153 = -t187 * pkin(5) - t199;
t1 = [t209, 0, 0, t194 ^ 2 * t209, t194 * t206, t202, t201, qJD(2) ^ 2 / 0.2e1, pkin(1) * t206 - pkin(7) * t202, -t198 * pkin(1) * t194 - pkin(7) * t201, -t174 * t181 + t175 * t180, t175 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1, t173 ^ 2 / 0.2e1, -t173 * t172, t173 * t188, -t172 * t188, t188 ^ 2 / 0.2e1, t176 * t172 + t200 * t188, t176 * t173 - t204 * t188, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t187, -t163 * t187, t187 ^ 2 / 0.2e1, t165 * t163 + t199 * t187, t165 * t164 - t205 * t187, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t162, -t160 * t162, t162 ^ 2 / 0.2e1 (-t191 * t154 + t195 * t155) * t162 + t153 * t160 -(t195 * t154 + t191 * t155) * t162 + t153 * t161;];
T_reg  = t1;
