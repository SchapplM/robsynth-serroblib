% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:17:45
% EndTime: 2019-03-09 18:17:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (631->54), mult. (1336->113), div. (0->0), fcn. (1007->10), ass. (0->49)
t213 = -pkin(8) - pkin(7);
t200 = qJD(1) ^ 2;
t212 = t200 / 0.2e1;
t211 = cos(qJ(5));
t199 = cos(qJ(2));
t210 = t199 * t200;
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t206 = qJD(1) * t199;
t196 = sin(qJ(2));
t207 = qJD(1) * t196;
t181 = t195 * t207 - t198 * t206;
t182 = (t195 * t199 + t196 * t198) * qJD(1);
t187 = (-pkin(2) * t199 - pkin(1)) * qJD(1);
t171 = t181 * pkin(3) - t182 * qJ(4) + t187;
t190 = qJD(2) + qJD(3);
t185 = qJD(2) * pkin(2) + t213 * t207;
t186 = t213 * t206;
t208 = t195 * t185 - t198 * t186;
t174 = t190 * qJ(4) + t208;
t191 = sin(pkin(11));
t192 = cos(pkin(11));
t163 = t192 * t171 - t191 * t174;
t177 = t192 * t182 + t191 * t190;
t157 = t181 * pkin(4) - t177 * pkin(9) + t163;
t164 = t191 * t171 + t192 * t174;
t176 = t191 * t182 - t192 * t190;
t162 = -t176 * pkin(9) + t164;
t194 = sin(qJ(5));
t209 = t194 * t157 + t211 * t162;
t205 = qJD(1) * qJD(2);
t204 = t196 * t205;
t203 = t199 * t205;
t202 = t211 * t157 - t194 * t162;
t201 = t198 * t185 + t195 * t186;
t180 = qJD(5) + t181;
t173 = -t190 * pkin(3) + qJD(4) - t201;
t165 = t176 * pkin(4) + t173;
t197 = cos(qJ(6));
t193 = sin(qJ(6));
t178 = qJD(6) + t180;
t168 = -t194 * t176 + t211 * t177;
t167 = t211 * t176 + t194 * t177;
t160 = -t193 * t167 + t197 * t168;
t159 = t197 * t167 + t193 * t168;
t158 = t167 * pkin(5) + t165;
t154 = -t167 * pkin(10) + t209;
t153 = t180 * pkin(5) - t168 * pkin(10) + t202;
t1 = [t212, 0, 0, t196 ^ 2 * t212, t196 * t210, t204, t203, qJD(2) ^ 2 / 0.2e1, pkin(1) * t210 - pkin(7) * t204, -t200 * pkin(1) * t196 - pkin(7) * t203, t182 ^ 2 / 0.2e1, -t182 * t181, t182 * t190, -t181 * t190, t190 ^ 2 / 0.2e1, t187 * t181 + t201 * t190, t187 * t182 - t208 * t190, t163 * t181 + t173 * t176, -t164 * t181 + t173 * t177, -t163 * t177 - t164 * t176, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t180, -t167 * t180, t180 ^ 2 / 0.2e1, t165 * t167 + t202 * t180, t165 * t168 - t209 * t180, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t178, -t159 * t178, t178 ^ 2 / 0.2e1 (t197 * t153 - t193 * t154) * t178 + t158 * t159 -(t193 * t153 + t197 * t154) * t178 + t158 * t160;];
T_reg  = t1;
