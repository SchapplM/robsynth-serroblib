% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:24
% EndTime: 2019-03-09 13:53:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (329->52), mult. (715->109), div. (0->0), fcn. (479->8), ass. (0->43)
t196 = qJD(1) ^ 2;
t208 = t196 / 0.2e1;
t195 = cos(qJ(2));
t207 = t195 * t196;
t190 = sin(qJ(4));
t191 = sin(qJ(2));
t194 = cos(qJ(4));
t171 = (-t190 * t195 + t191 * t194) * qJD(1);
t184 = qJD(2) - qJD(4);
t204 = qJD(1) * t191;
t202 = pkin(7) * t204 + qJD(3);
t167 = -pkin(8) * t204 + (-pkin(2) - pkin(3)) * qJD(2) + t202;
t203 = qJD(1) * t195;
t175 = pkin(7) * t203 + qJD(2) * qJ(3);
t172 = -pkin(8) * t203 + t175;
t198 = t194 * t167 - t190 * t172;
t155 = -t184 * pkin(4) - t171 * pkin(9) + t198;
t170 = (t190 * t191 + t194 * t195) * qJD(1);
t205 = t190 * t167 + t194 * t172;
t157 = -t170 * pkin(9) + t205;
t189 = sin(qJ(5));
t193 = cos(qJ(5));
t206 = t189 * t155 + t193 * t157;
t201 = qJD(1) * qJD(2);
t173 = -qJD(1) * pkin(1) - pkin(2) * t203 - qJ(3) * t204;
t200 = t191 * t201;
t199 = t195 * t201;
t161 = t193 * t170 + t189 * t171;
t166 = pkin(3) * t203 - t173;
t197 = t193 * t155 - t189 * t157;
t163 = t170 * pkin(4) + t166;
t192 = cos(qJ(6));
t188 = sin(qJ(6));
t182 = -qJD(5) + t184;
t174 = -qJD(2) * pkin(2) + t202;
t162 = -t189 * t170 + t193 * t171;
t160 = qJD(6) + t161;
t159 = t192 * t162 - t188 * t182;
t158 = t188 * t162 + t192 * t182;
t153 = t161 * pkin(5) - t162 * pkin(10) + t163;
t152 = -t182 * pkin(10) + t206;
t151 = t182 * pkin(5) - t197;
t1 = [t208, 0, 0, t191 ^ 2 * t208, t191 * t207, t200, t199, qJD(2) ^ 2 / 0.2e1, pkin(1) * t207 - pkin(7) * t200, -t196 * pkin(1) * t191 - pkin(7) * t199, -t174 * qJD(2) - t173 * t203 (t174 * t191 + t175 * t195) * qJD(1), t175 * qJD(2) - t173 * t204, t175 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t184, t170 * t184, t184 ^ 2 / 0.2e1, t166 * t170 - t198 * t184, t166 * t171 + t205 * t184, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t182, t161 * t182, t182 ^ 2 / 0.2e1, t163 * t161 - t197 * t182, t163 * t162 + t206 * t182, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t160, -t158 * t160, t160 ^ 2 / 0.2e1 (-t188 * t152 + t192 * t153) * t160 + t151 * t158 -(t192 * t152 + t188 * t153) * t160 + t151 * t159;];
T_reg  = t1;
