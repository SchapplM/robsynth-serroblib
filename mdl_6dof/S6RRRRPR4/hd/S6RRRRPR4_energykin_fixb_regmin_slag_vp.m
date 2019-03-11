% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:19
% EndTime: 2019-03-09 22:09:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (611->52), mult. (1283->109), div. (0->0), fcn. (965->10), ass. (0->49)
t212 = -pkin(8) - pkin(7);
t199 = qJD(1) ^ 2;
t211 = t199 / 0.2e1;
t210 = cos(qJ(4));
t198 = cos(qJ(2));
t209 = t198 * t199;
t194 = sin(qJ(3));
t195 = sin(qJ(2));
t197 = cos(qJ(3));
t181 = (t194 * t198 + t195 * t197) * qJD(1);
t189 = qJD(2) + qJD(3);
t193 = sin(qJ(4));
t176 = t210 * t181 + t193 * t189;
t205 = qJD(1) * t198;
t206 = qJD(1) * t195;
t180 = t194 * t206 - t197 * t205;
t179 = qJD(4) + t180;
t186 = (-pkin(2) * t198 - pkin(1)) * qJD(1);
t171 = t180 * pkin(3) - t181 * pkin(9) + t186;
t184 = qJD(2) * pkin(2) + t212 * t206;
t185 = t212 * t205;
t207 = t194 * t184 - t197 * t185;
t174 = t189 * pkin(9) + t207;
t201 = t210 * t171 - t193 * t174;
t159 = t179 * pkin(4) - t176 * qJ(5) + t201;
t175 = t193 * t181 - t210 * t189;
t208 = t193 * t171 + t210 * t174;
t164 = -t175 * qJ(5) + t208;
t190 = sin(pkin(11));
t191 = cos(pkin(11));
t156 = t190 * t159 + t191 * t164;
t204 = qJD(1) * qJD(2);
t203 = t195 * t204;
t202 = t198 * t204;
t155 = t191 * t159 - t190 * t164;
t200 = t197 * t184 + t194 * t185;
t173 = -t189 * pkin(3) - t200;
t165 = t175 * pkin(4) + qJD(5) + t173;
t196 = cos(qJ(6));
t192 = sin(qJ(6));
t177 = qJD(6) + t179;
t168 = -t190 * t175 + t191 * t176;
t167 = -t191 * t175 - t190 * t176;
t162 = t192 * t167 + t196 * t168;
t161 = -t196 * t167 + t192 * t168;
t160 = -t167 * pkin(5) + t165;
t154 = t167 * pkin(10) + t156;
t153 = t179 * pkin(5) - t168 * pkin(10) + t155;
t1 = [t211, 0, 0, t195 ^ 2 * t211, t195 * t209, t203, t202, qJD(2) ^ 2 / 0.2e1, pkin(1) * t209 - pkin(7) * t203, -t199 * pkin(1) * t195 - pkin(7) * t202, t181 ^ 2 / 0.2e1, -t181 * t180, t181 * t189, -t180 * t189, t189 ^ 2 / 0.2e1, t186 * t180 + t200 * t189, t186 * t181 - t207 * t189, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * t179, -t175 * t179, t179 ^ 2 / 0.2e1, t173 * t175 + t201 * t179, t173 * t176 - t208 * t179, -t155 * t168 + t156 * t167, t156 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t177, -t161 * t177, t177 ^ 2 / 0.2e1 (t196 * t153 - t192 * t154) * t177 + t160 * t161 -(t192 * t153 + t196 * t154) * t177 + t160 * t162;];
T_reg  = t1;
