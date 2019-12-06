% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:05
% EndTime: 2019-12-05 17:26:05
% DurationCPUTime: 0.15s
% Computational Cost: add. (198->38), mult. (509->89), div. (0->0), fcn. (405->12), ass. (0->41)
t197 = cos(qJ(2));
t209 = qJD(1) * sin(pkin(5));
t178 = qJD(2) * pkin(2) + t197 * t209;
t186 = sin(pkin(6));
t188 = cos(pkin(6));
t208 = qJD(1) * cos(pkin(5));
t199 = t178 * t188 + t186 * t208;
t198 = qJD(2) ^ 2;
t211 = t186 ^ 2 * t198;
t184 = t188 * qJD(2) + qJD(3);
t193 = sin(qJ(2));
t207 = qJD(2) * t186;
t176 = pkin(8) * t207 + t193 * t209;
t192 = sin(qJ(3));
t196 = cos(qJ(3));
t206 = t196 * t176 + t199 * t192;
t164 = t184 * pkin(9) + t206;
t183 = t188 * t208;
t166 = t183 + (-t178 + (-pkin(3) * t196 - pkin(9) * t192) * qJD(2)) * t186;
t191 = sin(qJ(4));
t195 = cos(qJ(4));
t210 = t195 * t164 + t191 * t166;
t204 = t192 * t207;
t203 = (-t186 * t178 + t183) * t207;
t202 = t196 * t207;
t201 = qJD(2) * t209;
t200 = -t191 * t164 + t195 * t166;
t171 = -t195 * t184 + t191 * t204;
t174 = t192 * t176;
t163 = -t184 * pkin(3) - t199 * t196 + t174;
t194 = cos(qJ(5));
t190 = sin(qJ(5));
t181 = -qJD(4) + t202;
t172 = t191 * t184 + t195 * t204;
t170 = qJD(5) + t171;
t168 = t194 * t172 - t190 * t181;
t167 = t190 * t172 + t194 * t181;
t161 = t171 * pkin(4) - t172 * pkin(10) + t163;
t160 = -t181 * pkin(10) + t210;
t159 = t181 * pkin(4) - t200;
t1 = [qJD(1) ^ 2 / 0.2e1, t198 / 0.2e1, t197 * t201, -t193 * t201, t192 ^ 2 * t211 / 0.2e1, t192 * t196 * t211, t184 * t204, t184 * t202, t184 ^ 2 / 0.2e1, -t174 * t184 + (t199 * t184 - t203) * t196, -t206 * t184 + t192 * t203, t172 ^ 2 / 0.2e1, -t172 * t171, -t172 * t181, t171 * t181, t181 ^ 2 / 0.2e1, t163 * t171 - t200 * t181, t163 * t172 + t210 * t181, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t170, -t167 * t170, t170 ^ 2 / 0.2e1, (-t190 * t160 + t194 * t161) * t170 + t159 * t167, -(t194 * t160 + t190 * t161) * t170 + t159 * t168;];
T_reg = t1;
