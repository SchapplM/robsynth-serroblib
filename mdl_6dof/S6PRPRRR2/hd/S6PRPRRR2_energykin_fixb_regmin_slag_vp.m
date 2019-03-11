% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR2
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:54
% EndTime: 2019-03-08 20:29:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (187->40), mult. (433->90), div. (0->0), fcn. (315->12), ass. (0->41)
t194 = qJD(2) ^ 2;
t206 = t194 / 0.2e1;
t205 = cos(qJ(5));
t193 = cos(qJ(2));
t202 = qJD(1) * sin(pkin(6));
t174 = qJD(2) * pkin(2) + t193 * t202;
t184 = sin(pkin(12));
t186 = cos(pkin(12));
t190 = sin(qJ(2));
t198 = t190 * t202;
t169 = t184 * t174 + t186 * t198;
t167 = qJD(2) * pkin(8) + t169;
t180 = cos(pkin(6)) * qJD(1) + qJD(3);
t189 = sin(qJ(4));
t192 = cos(qJ(4));
t203 = t192 * t167 + t189 * t180;
t158 = qJD(4) * pkin(9) + t203;
t168 = t186 * t174 - t184 * t198;
t161 = (-pkin(4) * t192 - pkin(9) * t189 - pkin(3)) * qJD(2) - t168;
t188 = sin(qJ(5));
t204 = t205 * t158 + t188 * t161;
t201 = qJD(2) * t189;
t200 = t192 * qJD(2);
t199 = qJD(2) * qJD(4);
t197 = qJD(2) * t202;
t196 = -t188 * t158 + t205 * t161;
t195 = -t189 * t167 + t192 * t180;
t181 = -qJD(5) + t200;
t157 = -qJD(4) * pkin(4) - t195;
t191 = cos(qJ(6));
t187 = sin(qJ(6));
t178 = -qJD(6) + t181;
t173 = t188 * qJD(4) + t205 * t201;
t172 = -t205 * qJD(4) + t188 * t201;
t166 = -qJD(2) * pkin(3) - t168;
t163 = -t187 * t172 + t191 * t173;
t162 = t191 * t172 + t187 * t173;
t155 = t172 * pkin(5) + t157;
t154 = -t172 * pkin(10) + t204;
t153 = -t181 * pkin(5) - t173 * pkin(10) + t196;
t1 = [qJD(1) ^ 2 / 0.2e1, t206, t193 * t197, -t190 * t197, t169 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1, t189 ^ 2 * t206, t189 * t194 * t192, t189 * t199, t192 * t199, qJD(4) ^ 2 / 0.2e1, t195 * qJD(4) - t166 * t200, -t203 * qJD(4) + t166 * t201, t173 ^ 2 / 0.2e1, -t173 * t172, -t173 * t181, t172 * t181, t181 ^ 2 / 0.2e1, t157 * t172 - t196 * t181, t157 * t173 + t204 * t181, t163 ^ 2 / 0.2e1, -t163 * t162, -t163 * t178, t162 * t178, t178 ^ 2 / 0.2e1 -(t191 * t153 - t187 * t154) * t178 + t155 * t162 (t187 * t153 + t191 * t154) * t178 + t155 * t163;];
T_reg  = t1;
