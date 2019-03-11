% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:48
% EndTime: 2019-03-09 14:32:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (341->52), mult. (717->108), div. (0->0), fcn. (473->8), ass. (0->45)
t215 = -pkin(2) - pkin(8);
t200 = qJD(1) ^ 2;
t214 = t200 / 0.2e1;
t213 = cos(qJ(5));
t199 = cos(qJ(2));
t212 = t199 * t200;
t195 = sin(qJ(4));
t198 = cos(qJ(4));
t209 = qJD(1) * t199;
t183 = t198 * qJD(2) - t195 * t209;
t196 = sin(qJ(2));
t208 = t196 * qJD(1);
t188 = qJD(4) + t208;
t203 = -qJ(3) * t196 - pkin(1);
t177 = (t215 * t199 + t203) * qJD(1);
t207 = pkin(7) * t208 + qJD(3);
t178 = pkin(3) * t208 + t215 * qJD(2) + t207;
t201 = -t195 * t177 + t198 * t178;
t166 = t188 * pkin(4) - t183 * pkin(9) + t201;
t182 = t195 * qJD(2) + t198 * t209;
t210 = t198 * t177 + t195 * t178;
t169 = -t182 * pkin(9) + t210;
t194 = sin(qJ(5));
t211 = t194 * t166 + t213 * t169;
t186 = -pkin(7) * t209 - qJD(2) * qJ(3);
t206 = qJD(1) * qJD(2);
t180 = pkin(3) * t209 - t186;
t205 = t196 * t206;
t204 = t199 * t206;
t202 = t213 * t166 - t194 * t169;
t187 = qJD(5) + t188;
t173 = t182 * pkin(4) + t180;
t197 = cos(qJ(6));
t193 = sin(qJ(6));
t185 = qJD(6) + t187;
t184 = -qJD(2) * pkin(2) + t207;
t181 = (-pkin(2) * t199 + t203) * qJD(1);
t172 = -t194 * t182 + t213 * t183;
t171 = t213 * t182 + t194 * t183;
t168 = t171 * pkin(5) + t173;
t165 = -t193 * t171 + t197 * t172;
t164 = t197 * t171 + t193 * t172;
t161 = -t171 * pkin(10) + t211;
t160 = t187 * pkin(5) - t172 * pkin(10) + t202;
t1 = [t214, 0, 0, t196 ^ 2 * t214, t196 * t212, t205, t204, qJD(2) ^ 2 / 0.2e1, pkin(1) * t212 - pkin(7) * t205, -t200 * pkin(1) * t196 - pkin(7) * t204 (t184 * t196 - t186 * t199) * qJD(1), t184 * qJD(2) + t181 * t209, -t186 * qJD(2) - t181 * t208, t181 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t184 ^ 2 / 0.2e1, t183 ^ 2 / 0.2e1, -t183 * t182, t183 * t188, -t182 * t188, t188 ^ 2 / 0.2e1, t180 * t182 + t201 * t188, t180 * t183 - t210 * t188, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t187, -t171 * t187, t187 ^ 2 / 0.2e1, t173 * t171 + t202 * t187, t173 * t172 - t211 * t187, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t185, -t164 * t185, t185 ^ 2 / 0.2e1 (t197 * t160 - t193 * t161) * t185 + t168 * t164 -(t193 * t160 + t197 * t161) * t185 + t168 * t165;];
T_reg  = t1;
