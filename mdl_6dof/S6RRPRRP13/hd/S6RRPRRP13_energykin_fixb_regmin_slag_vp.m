% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:32
% EndTime: 2019-03-09 13:01:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (367->52), mult. (872->107), div. (0->0), fcn. (620->8), ass. (0->45)
t214 = -pkin(2) - pkin(9);
t213 = cos(qJ(5));
t191 = sin(pkin(6));
t198 = qJD(1) ^ 2;
t212 = t191 ^ 2 * t198;
t195 = sin(qJ(2));
t208 = qJD(1) * t191;
t203 = t195 * t208;
t183 = qJD(4) + t203;
t185 = pkin(8) * t203;
t192 = cos(pkin(6));
t207 = t192 * qJD(1);
t189 = qJD(2) + t207;
t197 = cos(qJ(2));
t169 = qJD(3) + t185 + t214 * t189 + (-pkin(1) * t192 * t197 + pkin(3) * t191 * t195) * qJD(1);
t202 = -qJ(3) * t195 - pkin(1);
t175 = (t214 * t197 + t202) * t208;
t194 = sin(qJ(4));
t196 = cos(qJ(4));
t210 = t194 * t169 + t196 * t175;
t164 = t183 * pkin(10) + t210;
t204 = t197 * t208;
t206 = pkin(1) * t207;
t209 = pkin(8) * t204 + t195 * t206;
t177 = -t189 * qJ(3) - t209;
t174 = pkin(3) * t204 - t177;
t180 = t194 * t189 + t196 * t204;
t181 = t196 * t189 - t194 * t204;
t167 = t180 * pkin(4) - t181 * pkin(10) + t174;
t193 = sin(qJ(5));
t211 = t213 * t164 + t193 * t167;
t205 = t197 * t212;
t201 = -t193 * t164 + t213 * t167;
t200 = t196 * t169 - t194 * t175;
t199 = t197 * t206 - t185;
t163 = -t183 * pkin(4) - t200;
t179 = qJD(5) + t180;
t178 = (-pkin(2) * t197 + t202) * t208;
t176 = -t189 * pkin(2) + qJD(3) - t199;
t171 = t213 * t181 + t193 * t183;
t170 = t193 * t181 - t213 * t183;
t161 = t170 * pkin(5) + qJD(6) + t163;
t160 = -t170 * qJ(6) + t211;
t159 = t179 * pkin(5) - t171 * qJ(6) + t201;
t1 = [t198 / 0.2e1, 0, 0, t195 ^ 2 * t212 / 0.2e1, t195 * t205, t189 * t203, t189 * t204, t189 ^ 2 / 0.2e1, pkin(1) * t205 + t199 * t189, -pkin(1) * t195 * t212 - t209 * t189 (t176 * t195 - t177 * t197) * t208, t176 * t189 + t178 * t204, -t177 * t189 - t178 * t203, t178 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1, t181 ^ 2 / 0.2e1, -t181 * t180, t181 * t183, -t180 * t183, t183 ^ 2 / 0.2e1, t174 * t180 + t200 * t183, t174 * t181 - t210 * t183, t171 ^ 2 / 0.2e1, -t171 * t170, t171 * t179, -t170 * t179, t179 ^ 2 / 0.2e1, t163 * t170 + t201 * t179, t163 * t171 - t211 * t179, -t159 * t171 - t160 * t170, t160 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1;];
T_reg  = t1;
