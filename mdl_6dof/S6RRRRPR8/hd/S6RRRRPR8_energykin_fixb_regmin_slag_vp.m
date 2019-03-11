% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:47
% EndTime: 2019-03-09 22:44:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (454->52), mult. (955->106), div. (0->0), fcn. (673->8), ass. (0->45)
t213 = pkin(4) + pkin(5);
t197 = qJD(1) ^ 2;
t212 = t197 / 0.2e1;
t211 = cos(qJ(3));
t196 = cos(qJ(2));
t210 = t196 * t197;
t192 = sin(qJ(3));
t193 = sin(qJ(2));
t207 = qJD(1) * t193;
t178 = t192 * qJD(2) + t211 * t207;
t206 = t196 * qJD(1);
t186 = -qJD(3) + t206;
t176 = (-pkin(2) * t196 - pkin(8) * t193 - pkin(1)) * qJD(1);
t183 = pkin(7) * t206 + qJD(2) * pkin(8);
t201 = t211 * t176 - t192 * t183;
t165 = -t186 * pkin(3) - t178 * pkin(9) + t201;
t177 = -t211 * qJD(2) + t192 * t207;
t208 = t192 * t176 + t211 * t183;
t168 = -t177 * pkin(9) + t208;
t191 = sin(qJ(4));
t195 = cos(qJ(4));
t209 = t191 * t165 + t195 * t168;
t205 = qJD(1) * qJD(2);
t184 = -qJD(4) + t186;
t160 = -t184 * qJ(5) + t209;
t204 = t193 * t205;
t203 = t196 * t205;
t182 = -qJD(2) * pkin(2) + pkin(7) * t207;
t202 = t195 * t165 - t191 * t168;
t200 = qJD(5) - t202;
t199 = -t177 * pkin(3) - t182;
t171 = -t191 * t177 + t195 * t178;
t198 = t171 * qJ(5) + t199;
t194 = cos(qJ(6));
t190 = sin(qJ(6));
t181 = qJD(6) + t184;
t170 = t195 * t177 + t191 * t178;
t163 = t190 * t170 + t194 * t171;
t162 = -t194 * t170 + t190 * t171;
t161 = t170 * pkin(4) - t198;
t159 = t184 * pkin(4) + t200;
t158 = -t213 * t170 + t198;
t157 = t170 * pkin(10) + t160;
t156 = -t171 * pkin(10) + t213 * t184 + t200;
t1 = [t212, 0, 0, t193 ^ 2 * t212, t193 * t210, t204, t203, qJD(2) ^ 2 / 0.2e1, pkin(1) * t210 - pkin(7) * t204, -t197 * pkin(1) * t193 - pkin(7) * t203, t178 ^ 2 / 0.2e1, -t178 * t177, -t178 * t186, t177 * t186, t186 ^ 2 / 0.2e1, t182 * t177 - t201 * t186, t182 * t178 + t208 * t186, t171 ^ 2 / 0.2e1, -t171 * t170, -t171 * t184, t170 * t184, t184 ^ 2 / 0.2e1, -t170 * t199 - t202 * t184, -t171 * t199 + t209 * t184, t159 * t184 + t161 * t170, t159 * t171 - t160 * t170, -t160 * t184 - t161 * t171, t160 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t181, -t162 * t181, t181 ^ 2 / 0.2e1 (t194 * t156 - t190 * t157) * t181 + t158 * t162 -(t190 * t156 + t194 * t157) * t181 + t158 * t163;];
T_reg  = t1;
