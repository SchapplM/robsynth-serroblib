% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:02
% EndTime: 2019-03-09 17:49:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (426->50), mult. (996->101), div. (0->0), fcn. (734->8), ass. (0->44)
t218 = pkin(3) + pkin(10);
t217 = cos(qJ(5));
t194 = sin(pkin(6));
t201 = qJD(1) ^ 2;
t216 = t194 ^ 2 * t201;
t211 = cos(pkin(6)) * qJD(1);
t192 = qJD(2) + t211;
t197 = sin(qJ(3));
t199 = cos(qJ(3));
t198 = sin(qJ(2));
t212 = qJD(1) * t194;
t208 = t198 * t212;
t184 = t197 * t192 + t199 * t208;
t200 = cos(qJ(2));
t207 = t200 * t212;
t187 = -qJD(3) + t207;
t210 = pkin(1) * t211;
t213 = pkin(8) * t207 + t198 * t210;
t179 = t192 * pkin(9) + t213;
t181 = (-pkin(2) * t200 - pkin(9) * t198 - pkin(1)) * t212;
t205 = -t197 * t179 + t199 * t181;
t204 = qJD(4) - t205;
t166 = t184 * pkin(4) + t218 * t187 + t204;
t183 = -t199 * t192 + t197 * t208;
t203 = -pkin(8) * t208 + t200 * t210;
t178 = -t192 * pkin(2) - t203;
t202 = -t184 * qJ(4) + t178;
t168 = t218 * t183 + t202;
t196 = sin(qJ(5));
t215 = t196 * t166 + t217 * t168;
t214 = t199 * t179 + t197 * t181;
t209 = t200 * t216;
t172 = t187 * qJ(4) - t214;
t206 = t217 * t166 - t196 * t168;
t169 = -t183 * pkin(4) - t172;
t182 = qJD(5) + t184;
t174 = t196 * t183 - t217 * t187;
t173 = -t217 * t183 - t196 * t187;
t171 = t187 * pkin(3) + t204;
t170 = t183 * pkin(3) + t202;
t163 = t173 * pkin(5) + qJD(6) + t169;
t162 = -t173 * qJ(6) + t215;
t161 = t182 * pkin(5) - t174 * qJ(6) + t206;
t1 = [t201 / 0.2e1, 0, 0, t198 ^ 2 * t216 / 0.2e1, t198 * t209, t192 * t208, t192 * t207, t192 ^ 2 / 0.2e1, pkin(1) * t209 + t203 * t192, -pkin(1) * t198 * t216 - t213 * t192, t184 ^ 2 / 0.2e1, -t184 * t183, -t184 * t187, t183 * t187, t187 ^ 2 / 0.2e1, t178 * t183 - t205 * t187, t178 * t184 + t214 * t187, t171 * t184 + t172 * t183, -t170 * t183 - t171 * t187, -t170 * t184 + t172 * t187, t170 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t182, -t173 * t182, t182 ^ 2 / 0.2e1, t169 * t173 + t206 * t182, t169 * t174 - t215 * t182, -t161 * t174 - t162 * t173, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1;];
T_reg  = t1;
