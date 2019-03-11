% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:37
% EndTime: 2019-03-10 03:42:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (602->53), mult. (1245->112), div. (0->0), fcn. (963->10), ass. (0->50)
t215 = -pkin(8) - pkin(7);
t200 = qJD(1) ^ 2;
t214 = t200 / 0.2e1;
t213 = cos(qJ(5));
t199 = cos(qJ(2));
t212 = t199 * t200;
t194 = sin(qJ(3));
t195 = sin(qJ(2));
t198 = cos(qJ(3));
t182 = (t194 * t199 + t195 * t198) * qJD(1);
t190 = qJD(2) + qJD(3);
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t176 = t197 * t182 + t193 * t190;
t207 = qJD(1) * t199;
t208 = qJD(1) * t195;
t181 = t194 * t208 - t198 * t207;
t180 = qJD(4) + t181;
t187 = (-pkin(2) * t199 - pkin(1)) * qJD(1);
t170 = t181 * pkin(3) - t182 * pkin(9) + t187;
t185 = qJD(2) * pkin(2) + t215 * t208;
t186 = t215 * t207;
t209 = t194 * t185 - t198 * t186;
t173 = t190 * pkin(9) + t209;
t202 = t197 * t170 - t193 * t173;
t158 = t180 * pkin(4) - t176 * pkin(10) + t202;
t175 = t193 * t182 - t197 * t190;
t210 = t193 * t170 + t197 * t173;
t163 = -t175 * pkin(10) + t210;
t192 = sin(qJ(5));
t211 = t192 * t158 + t213 * t163;
t206 = qJD(1) * qJD(2);
t205 = t195 * t206;
t204 = t199 * t206;
t203 = t213 * t158 - t192 * t163;
t201 = t198 * t185 + t194 * t186;
t172 = -t190 * pkin(3) - t201;
t178 = qJD(5) + t180;
t164 = t175 * pkin(4) + t172;
t196 = cos(qJ(6));
t191 = sin(qJ(6));
t177 = qJD(6) + t178;
t167 = -t192 * t175 + t213 * t176;
t166 = t213 * t175 + t192 * t176;
t161 = -t191 * t166 + t196 * t167;
t160 = t196 * t166 + t191 * t167;
t159 = t166 * pkin(5) + t164;
t155 = -t166 * pkin(11) + t211;
t154 = t178 * pkin(5) - t167 * pkin(11) + t203;
t1 = [t214, 0, 0, t195 ^ 2 * t214, t195 * t212, t205, t204, qJD(2) ^ 2 / 0.2e1, pkin(1) * t212 - pkin(7) * t205, -t200 * pkin(1) * t195 - pkin(7) * t204, t182 ^ 2 / 0.2e1, -t182 * t181, t182 * t190, -t181 * t190, t190 ^ 2 / 0.2e1, t187 * t181 + t201 * t190, t187 * t182 - t209 * t190, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * t180, -t175 * t180, t180 ^ 2 / 0.2e1, t172 * t175 + t202 * t180, t172 * t176 - t210 * t180, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t178, -t166 * t178, t178 ^ 2 / 0.2e1, t164 * t166 + t203 * t178, t164 * t167 - t211 * t178, t161 ^ 2 / 0.2e1, -t161 * t160, t161 * t177, -t160 * t177, t177 ^ 2 / 0.2e1 (t196 * t154 - t191 * t155) * t177 + t159 * t160 -(t191 * t154 + t196 * t155) * t177 + t159 * t161;];
T_reg  = t1;
