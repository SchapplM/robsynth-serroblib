% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:10
% EndTime: 2019-03-08 22:34:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (191->46), mult. (429->97), div. (0->0), fcn. (284->10), ass. (0->42)
t205 = -pkin(3) - pkin(9);
t190 = qJD(2) ^ 2;
t204 = t190 / 0.2e1;
t185 = sin(qJ(2));
t201 = qJD(1) * sin(pkin(6));
t173 = qJD(2) * pkin(8) + t185 * t201;
t184 = sin(qJ(3));
t188 = cos(qJ(3));
t200 = qJD(1) * cos(pkin(6));
t192 = -t184 * t173 + t188 * t200;
t191 = qJD(4) - t192;
t198 = t184 * qJD(2);
t159 = pkin(4) * t198 + t205 * qJD(3) + t191;
t194 = -qJ(4) * t184 - pkin(2);
t189 = cos(qJ(2));
t196 = t189 * t201;
t164 = -t196 + (t205 * t188 + t194) * qJD(2);
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t203 = t183 * t159 + t187 * t164;
t202 = t188 * t173 + t184 * t200;
t199 = qJD(2) * t188;
t197 = qJD(2) * qJD(3);
t166 = -qJD(3) * qJ(4) - t202;
t195 = qJD(2) * t201;
t193 = t187 * t159 - t183 * t164;
t177 = qJD(5) + t198;
t162 = pkin(4) * t199 - t166;
t186 = cos(qJ(6));
t182 = sin(qJ(6));
t175 = qJD(6) + t177;
t174 = -qJD(2) * pkin(2) - t196;
t172 = t187 * qJD(3) - t183 * t199;
t171 = t183 * qJD(3) + t187 * t199;
t167 = -t196 + (-pkin(3) * t188 + t194) * qJD(2);
t165 = -qJD(3) * pkin(3) + t191;
t161 = -t182 * t171 + t186 * t172;
t160 = t186 * t171 + t182 * t172;
t156 = t171 * pkin(5) + t162;
t155 = -t171 * pkin(10) + t203;
t154 = t177 * pkin(5) - t172 * pkin(10) + t193;
t1 = [qJD(1) ^ 2 / 0.2e1, t204, t189 * t195, -t185 * t195, t184 ^ 2 * t204, t184 * t190 * t188, t184 * t197, t188 * t197, qJD(3) ^ 2 / 0.2e1, t192 * qJD(3) - t174 * t199, -t202 * qJD(3) + t174 * t198 (t165 * t184 - t166 * t188) * qJD(2), t165 * qJD(3) + t167 * t199, -t166 * qJD(3) - t167 * t198, t167 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t177, -t171 * t177, t177 ^ 2 / 0.2e1, t162 * t171 + t193 * t177, t162 * t172 - t203 * t177, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t175, -t160 * t175, t175 ^ 2 / 0.2e1 (t186 * t154 - t182 * t155) * t175 + t156 * t160 -(t182 * t154 + t186 * t155) * t175 + t156 * t161;];
T_reg  = t1;
