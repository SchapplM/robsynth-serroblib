% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR4
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:27
% EndTime: 2019-03-08 22:14:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (191->44), mult. (428->98), div. (0->0), fcn. (285->10), ass. (0->40)
t194 = qJD(2) ^ 2;
t207 = t194 / 0.2e1;
t189 = sin(qJ(2));
t204 = qJD(1) * sin(pkin(6));
t172 = qJD(2) * pkin(8) + t189 * t204;
t188 = sin(qJ(3));
t192 = cos(qJ(3));
t203 = qJD(1) * cos(pkin(6));
t197 = -t188 * t172 + t192 * t203;
t195 = qJD(4) - t197;
t202 = qJD(2) * t188;
t158 = -pkin(9) * t202 + (-pkin(3) - pkin(4)) * qJD(3) + t195;
t205 = t192 * t172 + t188 * t203;
t163 = qJD(3) * qJ(4) + t205;
t201 = qJD(2) * t192;
t160 = -pkin(9) * t201 + t163;
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t206 = t187 * t158 + t191 * t160;
t193 = cos(qJ(2));
t173 = -qJD(2) * pkin(2) - t193 * t204;
t200 = t173 * qJD(2);
t199 = qJD(2) * qJD(3);
t198 = qJD(2) * t204;
t166 = -pkin(3) * t201 - qJ(4) * t202 + t173;
t161 = pkin(4) * t201 - t166;
t196 = t191 * t158 - t187 * t160;
t168 = (t187 * t188 + t191 * t192) * qJD(2);
t190 = cos(qJ(6));
t186 = sin(qJ(6));
t180 = qJD(3) - qJD(5);
t169 = (-t187 * t192 + t188 * t191) * qJD(2);
t167 = qJD(6) + t168;
t165 = t190 * t169 - t186 * t180;
t164 = t186 * t169 + t190 * t180;
t162 = -qJD(3) * pkin(3) + t195;
t156 = t168 * pkin(5) - t169 * pkin(10) + t161;
t155 = -t180 * pkin(10) + t206;
t154 = t180 * pkin(5) - t196;
t1 = [qJD(1) ^ 2 / 0.2e1, t207, t193 * t198, -t189 * t198, t188 ^ 2 * t207, t188 * t194 * t192, t188 * t199, t192 * t199, qJD(3) ^ 2 / 0.2e1, t197 * qJD(3) - t192 * t200, -t205 * qJD(3) + t188 * t200, -t162 * qJD(3) - t166 * t201 (t162 * t188 + t163 * t192) * qJD(2), t163 * qJD(3) - t166 * t202, t163 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t180, t168 * t180, t180 ^ 2 / 0.2e1, t161 * t168 - t196 * t180, t161 * t169 + t206 * t180, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t167, -t164 * t167, t167 ^ 2 / 0.2e1 (-t186 * t155 + t190 * t156) * t167 + t154 * t164 -(t190 * t155 + t186 * t156) * t167 + t154 * t165;];
T_reg  = t1;
