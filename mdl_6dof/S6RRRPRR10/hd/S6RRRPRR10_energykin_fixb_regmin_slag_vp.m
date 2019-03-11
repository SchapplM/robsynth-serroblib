% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:55
% EndTime: 2019-03-09 19:21:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (370->52), mult. (773->106), div. (0->0), fcn. (530->8), ass. (0->43)
t197 = qJD(1) ^ 2;
t210 = t197 / 0.2e1;
t209 = cos(qJ(5));
t196 = cos(qJ(2));
t208 = t196 * t197;
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t193 = sin(qJ(2));
t205 = qJD(1) * t193;
t177 = t192 * qJD(2) + t195 * t205;
t204 = t196 * qJD(1);
t185 = -qJD(3) + t204;
t173 = (-pkin(2) * t196 - pkin(8) * t193 - pkin(1)) * qJD(1);
t182 = pkin(7) * t204 + qJD(2) * pkin(8);
t199 = t195 * t173 - t192 * t182;
t198 = qJD(4) - t199;
t161 = -t177 * pkin(9) + (pkin(3) + pkin(4)) * t185 + t198;
t206 = t192 * t173 + t195 * t182;
t166 = -t185 * qJ(4) + t206;
t176 = -t195 * qJD(2) + t192 * t205;
t163 = t176 * pkin(9) + t166;
t191 = sin(qJ(5));
t207 = t191 * t161 + t209 * t163;
t181 = -qJD(2) * pkin(2) + pkin(7) * t205;
t203 = qJD(1) * qJD(2);
t202 = t193 * t203;
t201 = t196 * t203;
t200 = t209 * t161 - t191 * t163;
t167 = t176 * pkin(3) - t177 * qJ(4) + t181;
t184 = qJD(5) + t185;
t164 = -t176 * pkin(4) - t167;
t194 = cos(qJ(6));
t190 = sin(qJ(6));
t180 = qJD(6) + t184;
t170 = t191 * t176 + t209 * t177;
t169 = -t209 * t176 + t191 * t177;
t165 = t185 * pkin(3) + t198;
t158 = -t190 * t169 + t194 * t170;
t157 = t194 * t169 + t190 * t170;
t156 = t169 * pkin(5) + t164;
t155 = -t169 * pkin(10) + t207;
t154 = t184 * pkin(5) - t170 * pkin(10) + t200;
t1 = [t210, 0, 0, t193 ^ 2 * t210, t193 * t208, t202, t201, qJD(2) ^ 2 / 0.2e1, pkin(1) * t208 - pkin(7) * t202, -t197 * pkin(1) * t193 - pkin(7) * t201, t177 ^ 2 / 0.2e1, -t177 * t176, -t177 * t185, t176 * t185, t185 ^ 2 / 0.2e1, t181 * t176 - t199 * t185, t181 * t177 + t206 * t185, t165 * t185 + t167 * t176, t165 * t177 - t166 * t176, -t166 * t185 - t167 * t177, t166 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t184, -t169 * t184, t184 ^ 2 / 0.2e1, t164 * t169 + t200 * t184, t164 * t170 - t207 * t184, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t180, -t157 * t180, t180 ^ 2 / 0.2e1 (t194 * t154 - t190 * t155) * t180 + t156 * t157 -(t190 * t154 + t194 * t155) * t180 + t156 * t158;];
T_reg  = t1;
