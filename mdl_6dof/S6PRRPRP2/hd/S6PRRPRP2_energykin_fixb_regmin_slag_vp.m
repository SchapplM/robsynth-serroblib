% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:14
% EndTime: 2019-03-08 21:32:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (312->44), mult. (729->94), div. (0->0), fcn. (532->10), ass. (0->40)
t198 = qJD(2) ^ 2;
t209 = t198 / 0.2e1;
t194 = sin(qJ(2));
t206 = qJD(1) * sin(pkin(6));
t182 = qJD(2) * pkin(8) + t194 * t206;
t196 = cos(qJ(3));
t205 = qJD(1) * cos(pkin(6));
t186 = t196 * t205;
t193 = sin(qJ(3));
t173 = qJD(3) * pkin(3) + t186 + (-qJ(4) * qJD(2) - t182) * t193;
t203 = qJD(2) * t196;
t207 = t196 * t182 + t193 * t205;
t174 = qJ(4) * t203 + t207;
t188 = sin(pkin(11));
t190 = cos(pkin(11));
t167 = t188 * t173 + t190 * t174;
t165 = qJD(3) * pkin(9) + t167;
t197 = cos(qJ(2));
t201 = t197 * t206;
t177 = -t201 + qJD(4) + (-pkin(3) * t196 - pkin(2)) * qJD(2);
t204 = qJD(2) * t193;
t179 = -t188 * t204 + t190 * t203;
t180 = (t188 * t196 + t190 * t193) * qJD(2);
t169 = -t179 * pkin(4) - t180 * pkin(9) + t177;
t192 = sin(qJ(5));
t195 = cos(qJ(5));
t208 = t195 * t165 + t192 * t169;
t202 = qJD(2) * qJD(3);
t200 = qJD(2) * t206;
t166 = t190 * t173 - t188 * t174;
t199 = -t192 * t165 + t195 * t169;
t164 = -qJD(3) * pkin(4) - t166;
t183 = -qJD(2) * pkin(2) - t201;
t178 = qJD(5) - t179;
t176 = t192 * qJD(3) + t195 * t180;
t175 = -t195 * qJD(3) + t192 * t180;
t162 = t175 * pkin(5) - t176 * qJ(6) + t164;
t161 = t178 * qJ(6) + t208;
t160 = -t178 * pkin(5) + qJD(6) - t199;
t1 = [qJD(1) ^ 2 / 0.2e1, t209, t197 * t200, -t194 * t200, t193 ^ 2 * t209, t193 * t198 * t196, t193 * t202, t196 * t202, qJD(3) ^ 2 / 0.2e1 (-t193 * t182 + t186) * qJD(3) - t183 * t203, -t207 * qJD(3) + t183 * t204, -t166 * t180 + t167 * t179, t167 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * t178, -t175 * t178, t178 ^ 2 / 0.2e1, t164 * t175 + t199 * t178, t164 * t176 - t208 * t178, -t160 * t178 + t162 * t175, t160 * t176 - t161 * t175, t161 * t178 - t162 * t176, t161 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1;];
T_reg  = t1;
