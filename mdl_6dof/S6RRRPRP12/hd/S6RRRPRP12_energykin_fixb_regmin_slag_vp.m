% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP12
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
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:16
% EndTime: 2019-03-09 17:59:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (541->52), mult. (1232->105), div. (0->0), fcn. (915->8), ass. (0->44)
t219 = pkin(3) + pkin(10);
t195 = sin(pkin(6));
t203 = qJD(1) ^ 2;
t218 = t195 ^ 2 * t203;
t213 = cos(pkin(6)) * qJD(1);
t193 = qJD(2) + t213;
t198 = sin(qJ(3));
t201 = cos(qJ(3));
t199 = sin(qJ(2));
t214 = qJD(1) * t195;
t210 = t199 * t214;
t185 = t198 * t193 + t201 * t210;
t202 = cos(qJ(2));
t209 = t202 * t214;
t188 = -qJD(3) + t209;
t212 = pkin(1) * t213;
t215 = pkin(8) * t209 + t199 * t212;
t180 = t193 * pkin(9) + t215;
t182 = (-pkin(2) * t202 - pkin(9) * t199 - pkin(1)) * t214;
t208 = -t198 * t180 + t201 * t182;
t207 = qJD(4) - t208;
t167 = t185 * pkin(4) + t219 * t188 + t207;
t184 = -t201 * t193 + t198 * t210;
t205 = -pkin(8) * t210 + t202 * t212;
t179 = -t193 * pkin(2) - t205;
t204 = -t185 * qJ(4) + t179;
t169 = t219 * t184 + t204;
t197 = sin(qJ(5));
t200 = cos(qJ(5));
t217 = t197 * t167 + t200 * t169;
t216 = t201 * t180 + t198 * t182;
t211 = t202 * t218;
t173 = t188 * qJ(4) - t216;
t206 = t200 * t167 - t197 * t169;
t170 = -t184 * pkin(4) - t173;
t183 = qJD(5) + t185;
t175 = t197 * t184 - t200 * t188;
t174 = -t200 * t184 - t197 * t188;
t172 = t188 * pkin(3) + t207;
t171 = t184 * pkin(3) + t204;
t165 = t174 * pkin(5) - t175 * qJ(6) + t170;
t164 = t183 * qJ(6) + t217;
t163 = -t183 * pkin(5) + qJD(6) - t206;
t1 = [t203 / 0.2e1, 0, 0, t199 ^ 2 * t218 / 0.2e1, t199 * t211, t193 * t210, t193 * t209, t193 ^ 2 / 0.2e1, pkin(1) * t211 + t205 * t193, -pkin(1) * t199 * t218 - t215 * t193, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t188, t184 * t188, t188 ^ 2 / 0.2e1, t179 * t184 - t208 * t188, t179 * t185 + t216 * t188, t172 * t185 + t173 * t184, -t171 * t184 - t172 * t188, -t171 * t185 + t173 * t188, t171 ^ 2 / 0.2e1 + t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t183, -t174 * t183, t183 ^ 2 / 0.2e1, t170 * t174 + t206 * t183, t170 * t175 - t217 * t183, -t163 * t183 + t165 * t174, t163 * t175 - t164 * t174, t164 * t183 - t165 * t175, t164 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1;];
T_reg  = t1;
