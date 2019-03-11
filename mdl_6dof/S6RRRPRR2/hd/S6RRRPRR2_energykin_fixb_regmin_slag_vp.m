% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR2
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
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:38
% EndTime: 2019-03-09 18:08:39
% DurationCPUTime: 0.16s
% Computational Cost: add. (556->52), mult. (1314->109), div. (0->0), fcn. (1008->10), ass. (0->49)
t216 = -pkin(8) - pkin(7);
t202 = qJD(1) ^ 2;
t215 = t202 / 0.2e1;
t214 = cos(qJ(3));
t213 = cos(qJ(5));
t201 = cos(qJ(2));
t212 = t201 * t202;
t198 = sin(qJ(3));
t199 = sin(qJ(2));
t186 = (t198 * t201 + t199 * t214) * qJD(1);
t193 = qJD(2) + qJD(3);
t209 = qJD(1) * t199;
t188 = qJD(2) * pkin(2) + t209 * t216;
t208 = qJD(1) * t201;
t189 = t216 * t208;
t203 = t188 * t214 + t189 * t198;
t170 = pkin(3) * t193 - qJ(4) * t186 + t203;
t185 = t198 * t209 - t208 * t214;
t210 = t188 * t198 - t189 * t214;
t174 = -qJ(4) * t185 + t210;
t194 = sin(pkin(11));
t195 = cos(pkin(11));
t163 = t170 * t194 + t174 * t195;
t161 = pkin(9) * t193 + t163;
t179 = -t185 * t195 - t186 * t194;
t180 = -t185 * t194 + t186 * t195;
t190 = (-pkin(2) * t201 - pkin(1)) * qJD(1);
t181 = t185 * pkin(3) + qJD(4) + t190;
t166 = -t179 * pkin(4) - t180 * pkin(9) + t181;
t197 = sin(qJ(5));
t211 = t161 * t213 + t166 * t197;
t207 = qJD(1) * qJD(2);
t206 = t199 * t207;
t205 = t201 * t207;
t204 = -t161 * t197 + t166 * t213;
t162 = t170 * t195 - t174 * t194;
t178 = qJD(5) - t179;
t160 = -pkin(4) * t193 - t162;
t200 = cos(qJ(6));
t196 = sin(qJ(6));
t177 = qJD(6) + t178;
t176 = t180 * t213 + t193 * t197;
t175 = t180 * t197 - t193 * t213;
t168 = -t175 * t196 + t176 * t200;
t167 = t175 * t200 + t176 * t196;
t158 = pkin(5) * t175 + t160;
t157 = -pkin(10) * t175 + t211;
t156 = pkin(5) * t178 - pkin(10) * t176 + t204;
t1 = [t215, 0, 0, t199 ^ 2 * t215, t199 * t212, t206, t205, qJD(2) ^ 2 / 0.2e1, pkin(1) * t212 - pkin(7) * t206, -pkin(1) * t199 * t202 - pkin(7) * t205, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t193, -t185 * t193, t193 ^ 2 / 0.2e1, t185 * t190 + t193 * t203, t186 * t190 - t193 * t210, -t162 * t180 + t163 * t179, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * t178, -t175 * t178, t178 ^ 2 / 0.2e1, t160 * t175 + t204 * t178, t160 * t176 - t178 * t211, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t177, -t167 * t177, t177 ^ 2 / 0.2e1 (t156 * t200 - t157 * t196) * t177 + t158 * t167 -(t156 * t196 + t157 * t200) * t177 + t158 * t168;];
T_reg  = t1;
