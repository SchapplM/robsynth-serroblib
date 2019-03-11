% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:21
% EndTime: 2019-03-08 23:03:21
% DurationCPUTime: 0.13s
% Computational Cost: add. (334->47), mult. (797->101), div. (0->0), fcn. (610->12), ass. (0->45)
t202 = qJD(2) ^ 2;
t214 = t202 / 0.2e1;
t213 = cos(qJ(4));
t196 = sin(qJ(4));
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t182 = (t196 * t200 + t213 * t197) * qJD(2);
t190 = qJD(3) + qJD(4);
t198 = sin(qJ(2));
t210 = qJD(1) * sin(pkin(6));
t184 = qJD(2) * pkin(8) + t198 * t210;
t209 = qJD(1) * cos(pkin(6));
t187 = t200 * t209;
t177 = qJD(3) * pkin(3) + t187 + (-pkin(9) * qJD(2) - t184) * t197;
t207 = qJD(2) * t200;
t211 = t200 * t184 + t197 * t209;
t178 = pkin(9) * t207 + t211;
t203 = t213 * t177 - t196 * t178;
t165 = t190 * pkin(4) - t182 * qJ(5) + t203;
t208 = qJD(2) * t197;
t181 = t196 * t208 - t213 * t207;
t212 = t196 * t177 + t213 * t178;
t167 = -t181 * qJ(5) + t212;
t191 = sin(pkin(12));
t193 = cos(pkin(12));
t162 = t191 * t165 + t193 * t167;
t206 = qJD(2) * qJD(3);
t201 = cos(qJ(2));
t205 = t201 * t210;
t204 = qJD(2) * t210;
t171 = -t193 * t181 - t191 * t182;
t161 = t193 * t165 - t191 * t167;
t180 = -t205 + (-pkin(3) * t200 - pkin(2)) * qJD(2);
t173 = t181 * pkin(4) + qJD(5) + t180;
t199 = cos(qJ(6));
t195 = sin(qJ(6));
t185 = -qJD(2) * pkin(2) - t205;
t172 = -t191 * t181 + t193 * t182;
t170 = qJD(6) - t171;
t169 = t199 * t172 + t195 * t190;
t168 = t195 * t172 - t199 * t190;
t163 = -t171 * pkin(5) - t172 * pkin(10) + t173;
t160 = t190 * pkin(10) + t162;
t159 = -t190 * pkin(5) - t161;
t1 = [qJD(1) ^ 2 / 0.2e1, t214, t201 * t204, -t198 * t204, t197 ^ 2 * t214, t197 * t202 * t200, t197 * t206, t200 * t206, qJD(3) ^ 2 / 0.2e1 (-t197 * t184 + t187) * qJD(3) - t185 * t207, -t211 * qJD(3) + t185 * t208, t182 ^ 2 / 0.2e1, -t182 * t181, t182 * t190, -t181 * t190, t190 ^ 2 / 0.2e1, t180 * t181 + t203 * t190, t180 * t182 - t212 * t190, -t161 * t172 + t162 * t171, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t170, -t168 * t170, t170 ^ 2 / 0.2e1 (-t195 * t160 + t199 * t163) * t170 + t159 * t168 -(t199 * t160 + t195 * t163) * t170 + t159 * t169;];
T_reg  = t1;
