% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP3
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:20
% EndTime: 2019-03-08 21:38:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (348->45), mult. (781->96), div. (0->0), fcn. (557->10), ass. (0->40)
t195 = qJD(2) ^ 2;
t207 = t195 / 0.2e1;
t191 = sin(qJ(2));
t204 = qJD(1) * sin(pkin(6));
t179 = qJD(2) * pkin(8) + t191 * t204;
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t203 = qJD(1) * cos(pkin(6));
t205 = t193 * t179 + t190 * t203;
t172 = qJD(3) * qJ(4) + t205;
t194 = cos(qJ(2));
t199 = t194 * t204;
t173 = -t199 + (-pkin(3) * t193 - qJ(4) * t190 - pkin(2)) * qJD(2);
t185 = sin(pkin(11));
t187 = cos(pkin(11));
t163 = -t185 * t172 + t187 * t173;
t202 = qJD(2) * t190;
t178 = t185 * qJD(3) + t187 * t202;
t201 = t193 * qJD(2);
t160 = -pkin(4) * t201 - t178 * pkin(9) + t163;
t164 = t187 * t172 + t185 * t173;
t177 = -t187 * qJD(3) + t185 * t202;
t162 = -t177 * pkin(9) + t164;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t206 = t189 * t160 + t192 * t162;
t200 = qJD(2) * qJD(3);
t198 = qJD(2) * t204;
t197 = -t190 * t179 + t193 * t203;
t196 = t192 * t160 - t189 * t162;
t169 = -qJD(3) * pkin(3) + qJD(4) - t197;
t165 = t177 * pkin(4) + t169;
t182 = -qJD(5) + t201;
t180 = -qJD(2) * pkin(2) - t199;
t167 = -t189 * t177 + t192 * t178;
t166 = t192 * t177 + t189 * t178;
t158 = t166 * pkin(5) - t167 * qJ(6) + t165;
t157 = -t182 * qJ(6) + t206;
t156 = t182 * pkin(5) + qJD(6) - t196;
t1 = [qJD(1) ^ 2 / 0.2e1, t207, t194 * t198, -t191 * t198, t190 ^ 2 * t207, t190 * t195 * t193, t190 * t200, t193 * t200, qJD(3) ^ 2 / 0.2e1, qJD(3) * t197 - t180 * t201, -qJD(3) * t205 + t180 * t202, -t163 * t201 + t169 * t177, t164 * t201 + t169 * t178, -t163 * t178 - t164 * t177, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t182, t166 * t182, t182 ^ 2 / 0.2e1, t165 * t166 - t182 * t196, t165 * t167 + t182 * t206, t156 * t182 + t158 * t166, t156 * t167 - t157 * t166, -t157 * t182 - t158 * t167, t157 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1;];
T_reg  = t1;
