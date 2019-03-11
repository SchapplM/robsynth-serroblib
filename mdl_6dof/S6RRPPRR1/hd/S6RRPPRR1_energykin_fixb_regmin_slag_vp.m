% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:03
% EndTime: 2019-03-09 08:48:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (330->50), mult. (801->102), div. (0->0), fcn. (553->8), ass. (0->43)
t199 = qJD(1) ^ 2;
t210 = t199 / 0.2e1;
t209 = pkin(7) + qJ(3);
t198 = cos(qJ(2));
t208 = t198 * t199;
t190 = sin(pkin(10));
t191 = cos(pkin(10));
t195 = sin(qJ(2));
t178 = (t190 * t198 + t191 * t195) * qJD(1);
t206 = qJD(1) * t195;
t181 = qJD(2) * pkin(2) - t209 * t206;
t205 = qJD(1) * t198;
t182 = t209 * t205;
t171 = t191 * t181 - t190 * t182;
t201 = qJD(4) - t171;
t160 = -t178 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(2) + t201;
t172 = t190 * t181 + t191 * t182;
t170 = qJD(2) * qJ(4) + t172;
t177 = t190 * t206 - t191 * t205;
t162 = t177 * pkin(8) + t170;
t194 = sin(qJ(5));
t197 = cos(qJ(5));
t207 = t194 * t160 + t197 * t162;
t204 = qJD(1) * qJD(2);
t183 = -qJD(1) * pkin(1) - pkin(2) * t205 + qJD(3);
t203 = t195 * t204;
t202 = t198 * t204;
t168 = -t197 * t177 + t194 * t178;
t165 = t177 * pkin(3) - t178 * qJ(4) + t183;
t200 = t197 * t160 - t194 * t162;
t158 = -t177 * pkin(4) - t165;
t196 = cos(qJ(6));
t193 = sin(qJ(6));
t187 = qJD(2) - qJD(5);
t169 = t194 * t177 + t197 * t178;
t167 = -qJD(2) * pkin(3) + t201;
t166 = qJD(6) + t168;
t164 = t196 * t169 - t193 * t187;
t163 = t193 * t169 + t196 * t187;
t157 = -t187 * pkin(9) + t207;
t156 = t187 * pkin(5) - t200;
t155 = t168 * pkin(5) - t169 * pkin(9) + t158;
t1 = [t210, 0, 0, t195 ^ 2 * t210, t195 * t208, t203, t202, qJD(2) ^ 2 / 0.2e1, pkin(1) * t208 - pkin(7) * t203, -t199 * pkin(1) * t195 - pkin(7) * t202, -t171 * t178 - t172 * t177, t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1, -t167 * qJD(2) + t165 * t177, t167 * t178 - t170 * t177, t170 * qJD(2) - t165 * t178, t170 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t187, t168 * t187, t187 ^ 2 / 0.2e1, t158 * t168 - t200 * t187, t158 * t169 + t207 * t187, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t166, -t163 * t166, t166 ^ 2 / 0.2e1 (t196 * t155 - t193 * t157) * t166 + t156 * t163 -(t193 * t155 + t196 * t157) * t166 + t156 * t164;];
T_reg  = t1;
