% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP14
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:50
% EndTime: 2019-03-09 13:09:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (480->54), mult. (1104->111), div. (0->0), fcn. (791->8), ass. (0->45)
t218 = -pkin(2) - pkin(9);
t195 = sin(pkin(6));
t203 = qJD(1) ^ 2;
t217 = t195 ^ 2 * t203;
t199 = sin(qJ(2));
t213 = qJD(1) * t195;
t208 = t199 * t213;
t187 = qJD(4) + t208;
t189 = pkin(8) * t208;
t196 = cos(pkin(6));
t212 = t196 * qJD(1);
t193 = qJD(2) + t212;
t202 = cos(qJ(2));
t173 = qJD(3) + t189 + t218 * t193 + (-pkin(1) * t196 * t202 + pkin(3) * t195 * t199) * qJD(1);
t207 = -qJ(3) * t199 - pkin(1);
t179 = (t218 * t202 + t207) * t213;
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t215 = t198 * t173 + t201 * t179;
t169 = t187 * pkin(10) + t215;
t209 = t202 * t213;
t211 = pkin(1) * t212;
t214 = pkin(8) * t209 + t199 * t211;
t181 = -t193 * qJ(3) - t214;
t178 = pkin(3) * t209 - t181;
t184 = t198 * t193 + t201 * t209;
t185 = t201 * t193 - t198 * t209;
t171 = t184 * pkin(4) - t185 * pkin(10) + t178;
t197 = sin(qJ(5));
t200 = cos(qJ(5));
t216 = t200 * t169 + t197 * t171;
t210 = t202 * t217;
t206 = t201 * t173 - t198 * t179;
t205 = -t197 * t169 + t200 * t171;
t204 = t202 * t211 - t189;
t168 = -t187 * pkin(4) - t206;
t183 = qJD(5) + t184;
t182 = (-pkin(2) * t202 + t207) * t213;
t180 = -t193 * pkin(2) + qJD(3) - t204;
t175 = t200 * t185 + t197 * t187;
t174 = t197 * t185 - t200 * t187;
t166 = t174 * pkin(5) - t175 * qJ(6) + t168;
t165 = t183 * qJ(6) + t216;
t164 = -t183 * pkin(5) + qJD(6) - t205;
t1 = [t203 / 0.2e1, 0, 0, t199 ^ 2 * t217 / 0.2e1, t199 * t210, t193 * t208, t193 * t209, t193 ^ 2 / 0.2e1, pkin(1) * t210 + t204 * t193, -pkin(1) * t199 * t217 - t214 * t193 (t180 * t199 - t181 * t202) * t213, t180 * t193 + t182 * t209, -t181 * t193 - t182 * t208, t182 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1, t185 ^ 2 / 0.2e1, -t185 * t184, t185 * t187, -t184 * t187, t187 ^ 2 / 0.2e1, t178 * t184 + t206 * t187, t178 * t185 - t215 * t187, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t183, -t174 * t183, t183 ^ 2 / 0.2e1, t168 * t174 + t205 * t183, t168 * t175 - t216 * t183, -t164 * t183 + t166 * t174, t164 * t175 - t165 * t174, t165 * t183 - t166 * t175, t165 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1;];
T_reg  = t1;
