% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:08
% EndTime: 2019-03-09 13:19:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (536->52), mult. (1301->109), div. (0->0), fcn. (1010->10), ass. (0->47)
t215 = qJD(1) * (pkin(7) + qJ(3));
t203 = qJD(1) ^ 2;
t214 = t203 / 0.2e1;
t213 = cos(qJ(5));
t202 = cos(qJ(2));
t211 = t202 * t203;
t193 = qJD(2) + qJD(4);
t199 = sin(qJ(2));
t189 = qJD(2) * pkin(2) - t199 * t215;
t190 = t202 * t215;
t194 = sin(pkin(11));
t195 = cos(pkin(11));
t180 = t195 * t189 - t194 * t190;
t187 = (t194 * t202 + t195 * t199) * qJD(1);
t171 = qJD(2) * pkin(3) - t187 * pkin(8) + t180;
t181 = t194 * t189 + t195 * t190;
t186 = (-t194 * t199 + t195 * t202) * qJD(1);
t173 = t186 * pkin(8) + t181;
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t209 = t198 * t171 + t201 * t173;
t162 = t193 * pkin(9) + t209;
t178 = -t201 * t186 + t198 * t187;
t179 = t198 * t186 + t201 * t187;
t191 = qJD(3) + (-pkin(2) * t202 - pkin(1)) * qJD(1);
t182 = -t186 * pkin(3) + t191;
t165 = t178 * pkin(4) - t179 * pkin(9) + t182;
t197 = sin(qJ(5));
t210 = t213 * t162 + t197 * t165;
t208 = qJD(1) * qJD(2);
t207 = t199 * t208;
t206 = t202 * t208;
t205 = -t197 * t162 + t213 * t165;
t204 = t201 * t171 - t198 * t173;
t177 = qJD(5) + t178;
t161 = -t193 * pkin(4) - t204;
t200 = cos(qJ(6));
t196 = sin(qJ(6));
t176 = qJD(6) + t177;
t175 = t213 * t179 + t197 * t193;
t174 = t197 * t179 - t213 * t193;
t167 = -t196 * t174 + t200 * t175;
t166 = t200 * t174 + t196 * t175;
t159 = t174 * pkin(5) + t161;
t158 = -t174 * pkin(10) + t210;
t157 = t177 * pkin(5) - t175 * pkin(10) + t205;
t1 = [t214, 0, 0, t199 ^ 2 * t214, t199 * t211, t207, t206, qJD(2) ^ 2 / 0.2e1, pkin(1) * t211 - pkin(7) * t207, -t203 * pkin(1) * t199 - pkin(7) * t206, -t180 * t187 + t181 * t186, t181 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t193, -t178 * t193, t193 ^ 2 / 0.2e1, t182 * t178 + t204 * t193, t182 * t179 - t209 * t193, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t177, -t174 * t177, t177 ^ 2 / 0.2e1, t161 * t174 + t205 * t177, t161 * t175 - t210 * t177, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t176, -t166 * t176, t176 ^ 2 / 0.2e1 (t200 * t157 - t196 * t158) * t176 + t159 * t166 -(t196 * t157 + t200 * t158) * t176 + t159 * t167;];
T_reg  = t1;
