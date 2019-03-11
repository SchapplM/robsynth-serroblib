% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:05
% EndTime: 2019-03-08 20:25:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (191->41), mult. (442->92), div. (0->0), fcn. (324->12), ass. (0->41)
t192 = qJD(2) ^ 2;
t202 = t192 / 0.2e1;
t191 = cos(qJ(2));
t199 = qJD(1) * sin(pkin(6));
t172 = qJD(2) * pkin(2) + t191 * t199;
t181 = sin(pkin(12));
t183 = cos(pkin(12));
t187 = sin(qJ(2));
t195 = t187 * t199;
t165 = t181 * t172 + t183 * t195;
t163 = qJD(2) * pkin(8) + t165;
t178 = cos(pkin(6)) * qJD(1) + qJD(3);
t190 = cos(qJ(4));
t176 = t190 * t178;
t186 = sin(qJ(4));
t158 = qJD(4) * pkin(4) + t176 + (-pkin(9) * qJD(2) - t163) * t186;
t197 = qJD(2) * t190;
t200 = t190 * t163 + t186 * t178;
t159 = pkin(9) * t197 + t200;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t201 = t185 * t158 + t189 * t159;
t198 = qJD(2) * t186;
t196 = qJD(2) * qJD(4);
t194 = qJD(2) * t199;
t164 = t183 * t172 - t181 * t195;
t169 = t185 * t198 - t189 * t197;
t193 = t189 * t158 - t185 * t159;
t160 = (-pkin(4) * t190 - pkin(3)) * qJD(2) - t164;
t188 = cos(qJ(6));
t184 = sin(qJ(6));
t180 = qJD(4) + qJD(5);
t170 = (t185 * t190 + t186 * t189) * qJD(2);
t168 = qJD(6) + t169;
t167 = t188 * t170 + t184 * t180;
t166 = t184 * t170 - t188 * t180;
t162 = -qJD(2) * pkin(3) - t164;
t155 = t169 * pkin(5) - t170 * pkin(10) + t160;
t154 = t180 * pkin(10) + t201;
t153 = -t180 * pkin(5) - t193;
t1 = [qJD(1) ^ 2 / 0.2e1, t202, t191 * t194, -t187 * t194, t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1, t186 ^ 2 * t202, t186 * t192 * t190, t186 * t196, t190 * t196, qJD(4) ^ 2 / 0.2e1 (-t186 * t163 + t176) * qJD(4) - t162 * t197, -t200 * qJD(4) + t162 * t198, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t180, -t169 * t180, t180 ^ 2 / 0.2e1, t160 * t169 + t193 * t180, t160 * t170 - t201 * t180, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t168, -t166 * t168, t168 ^ 2 / 0.2e1 (-t184 * t154 + t188 * t155) * t168 + t153 * t166 -(t188 * t154 + t184 * t155) * t168 + t153 * t167;];
T_reg  = t1;
