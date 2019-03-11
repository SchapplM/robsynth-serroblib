% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:25
% EndTime: 2019-03-08 21:16:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (250->46), mult. (575->96), div. (0->0), fcn. (394->10), ass. (0->39)
t194 = qJD(2) ^ 2;
t206 = t194 / 0.2e1;
t190 = sin(qJ(2));
t203 = qJD(1) * sin(pkin(6));
t177 = qJD(2) * pkin(8) + t190 * t203;
t189 = sin(qJ(3));
t192 = cos(qJ(3));
t202 = qJD(1) * cos(pkin(6));
t204 = t192 * t177 + t189 * t202;
t170 = qJD(3) * qJ(4) + t204;
t193 = cos(qJ(2));
t198 = t193 * t203;
t171 = -t198 + (-pkin(3) * t192 - qJ(4) * t189 - pkin(2)) * qJD(2);
t184 = sin(pkin(11));
t186 = cos(pkin(11));
t163 = t186 * t170 + t184 * t171;
t205 = -t189 * t177 + t192 * t202;
t201 = qJD(2) * t189;
t200 = t192 * qJD(2);
t199 = qJD(2) * qJD(3);
t197 = qJD(2) * t203;
t162 = -t184 * t170 + t186 * t171;
t196 = qJD(3) * pkin(3) - qJD(4) + t205;
t160 = -qJ(5) * t200 + t163;
t159 = pkin(4) * t200 + qJD(5) - t162;
t176 = t184 * qJD(3) + t186 * t201;
t195 = t176 * qJ(5) + t196;
t191 = cos(qJ(6));
t188 = sin(qJ(6));
t181 = qJD(6) + t200;
t178 = -qJD(2) * pkin(2) - t198;
t175 = -t186 * qJD(3) + t184 * t201;
t165 = t188 * t175 + t191 * t176;
t164 = -t191 * t175 + t188 * t176;
t161 = t175 * pkin(4) - t195;
t158 = (-pkin(4) - pkin(5)) * t175 + t195;
t157 = t175 * pkin(9) + t160;
t156 = pkin(5) * t200 - t176 * pkin(9) + t159;
t1 = [qJD(1) ^ 2 / 0.2e1, t206, t193 * t197, -t190 * t197, t189 ^ 2 * t206, t189 * t194 * t192, t189 * t199, t192 * t199, qJD(3) ^ 2 / 0.2e1, t205 * qJD(3) - t178 * t200, -t204 * qJD(3) + t178 * t201, -t162 * t200 - t175 * t196, t163 * t200 - t176 * t196, -t162 * t176 - t163 * t175, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t196 ^ 2 / 0.2e1, t159 * t200 + t161 * t175, t159 * t176 - t160 * t175, -t160 * t200 - t161 * t176, t160 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t181, -t164 * t181, t181 ^ 2 / 0.2e1 (t156 * t191 - t157 * t188) * t181 + t158 * t164 -(t156 * t188 + t157 * t191) * t181 + t158 * t165;];
T_reg  = t1;
