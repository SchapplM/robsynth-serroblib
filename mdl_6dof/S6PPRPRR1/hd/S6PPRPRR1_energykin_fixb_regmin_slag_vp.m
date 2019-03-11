% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:37
% EndTime: 2019-03-08 18:43:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (171->35), mult. (452->86), div. (0->0), fcn. (377->14), ass. (0->40)
t195 = qJD(3) ^ 2;
t206 = t195 / 0.2e1;
t177 = cos(pkin(6)) * qJD(1) + qJD(2);
t184 = sin(pkin(7));
t205 = t177 * t184;
t183 = sin(pkin(12));
t188 = cos(pkin(7));
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t187 = cos(pkin(12));
t185 = sin(pkin(6));
t203 = qJD(1) * t185;
t199 = t187 * t203;
t197 = -t191 * t183 * t203 + (t188 * t199 + t205) * t194;
t169 = qJD(3) * pkin(3) + t197;
t170 = t191 * t205 + (t187 * t188 * t191 + t183 * t194) * t203;
t182 = sin(pkin(13));
t186 = cos(pkin(13));
t165 = t182 * t169 + t186 * t170;
t163 = qJD(3) * pkin(9) + t165;
t172 = t188 * t177 - t184 * t199 + qJD(4);
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t204 = t193 * t163 + t190 * t172;
t202 = qJD(3) * t190;
t201 = t193 * qJD(3);
t200 = qJD(3) * qJD(5);
t164 = t186 * t169 - t182 * t170;
t198 = -t190 * t163 + t193 * t172;
t196 = qJD(1) ^ 2;
t192 = cos(qJ(6));
t189 = sin(qJ(6));
t178 = -qJD(6) + t201;
t176 = t189 * qJD(5) + t192 * t202;
t175 = -t192 * qJD(5) + t189 * t202;
t162 = -qJD(3) * pkin(4) - t164;
t160 = (-pkin(5) * t193 - pkin(10) * t190 - pkin(4)) * qJD(3) - t164;
t159 = qJD(5) * pkin(10) + t204;
t158 = -qJD(5) * pkin(5) - t198;
t1 = [t196 / 0.2e1, t177 ^ 2 / 0.2e1 + (t183 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1) * t196 * t185 ^ 2, t206, t197 * qJD(3), -t170 * qJD(3), t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t190 ^ 2 * t206, t190 * t195 * t193, t190 * t200, t193 * t200, qJD(5) ^ 2 / 0.2e1, t198 * qJD(5) - t162 * t201, -t204 * qJD(5) + t162 * t202, t176 ^ 2 / 0.2e1, -t176 * t175, -t176 * t178, t175 * t178, t178 ^ 2 / 0.2e1 -(-t189 * t159 + t192 * t160) * t178 + t158 * t175 (t192 * t159 + t189 * t160) * t178 + t158 * t176;];
T_reg  = t1;
