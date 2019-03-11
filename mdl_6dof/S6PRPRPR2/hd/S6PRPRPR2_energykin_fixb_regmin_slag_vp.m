% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:53
% EndTime: 2019-03-08 19:32:53
% DurationCPUTime: 0.09s
% Computational Cost: add. (234->41), mult. (538->91), div. (0->0), fcn. (388->12), ass. (0->40)
t193 = qJD(2) ^ 2;
t203 = t193 / 0.2e1;
t202 = cos(pkin(12));
t192 = cos(qJ(2));
t200 = qJD(1) * sin(pkin(6));
t175 = qJD(2) * pkin(2) + t192 * t200;
t184 = sin(pkin(11));
t186 = cos(pkin(11));
t189 = sin(qJ(2));
t196 = t189 * t200;
t170 = t184 * t175 + t186 * t196;
t168 = qJD(2) * pkin(8) + t170;
t179 = cos(pkin(6)) * qJD(1) + qJD(3);
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t201 = t191 * t168 + t188 * t179;
t159 = qJD(4) * qJ(5) + t201;
t169 = t186 * t175 - t184 * t196;
t162 = (-pkin(4) * t191 - qJ(5) * t188 - pkin(3)) * qJD(2) - t169;
t183 = sin(pkin(12));
t155 = t202 * t159 + t183 * t162;
t199 = qJD(2) * t188;
t198 = t191 * qJD(2);
t197 = qJD(2) * qJD(4);
t195 = qJD(2) * t200;
t154 = -t183 * t159 + t202 * t162;
t194 = -t188 * t168 + t191 * t179;
t158 = -qJD(4) * pkin(4) + qJD(5) - t194;
t190 = cos(qJ(6));
t187 = sin(qJ(6));
t180 = -qJD(6) + t198;
t174 = t183 * qJD(4) + t202 * t199;
t173 = -t202 * qJD(4) + t183 * t199;
t167 = -qJD(2) * pkin(3) - t169;
t164 = -t187 * t173 + t190 * t174;
t163 = t190 * t173 + t187 * t174;
t156 = t173 * pkin(5) + t158;
t153 = -t173 * pkin(9) + t155;
t152 = -pkin(5) * t198 - t174 * pkin(9) + t154;
t1 = [qJD(1) ^ 2 / 0.2e1, t203, t192 * t195, -t189 * t195, t170 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1, t188 ^ 2 * t203, t188 * t193 * t191, t188 * t197, t191 * t197, qJD(4) ^ 2 / 0.2e1, t194 * qJD(4) - t167 * t198, -t201 * qJD(4) + t167 * t199, -t154 * t198 + t158 * t173, t155 * t198 + t158 * t174, -t154 * t174 - t155 * t173, t155 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, t164 ^ 2 / 0.2e1, -t164 * t163, -t164 * t180, t163 * t180, t180 ^ 2 / 0.2e1 -(t190 * t152 - t187 * t153) * t180 + t156 * t163 (t187 * t152 + t190 * t153) * t180 + t156 * t164;];
T_reg  = t1;
