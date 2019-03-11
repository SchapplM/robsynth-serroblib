% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:28
% EndTime: 2019-03-09 21:38:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (637->52), mult. (1443->105), div. (0->0), fcn. (1095->8), ass. (0->44)
t220 = -pkin(4) - pkin(5);
t219 = cos(qJ(4));
t196 = sin(pkin(6));
t203 = qJD(1) ^ 2;
t218 = t196 ^ 2 * t203;
t212 = cos(pkin(6)) * qJD(1);
t194 = qJD(2) + t212;
t202 = cos(qJ(2));
t200 = sin(qJ(2));
t213 = qJD(1) * t196;
t209 = t200 * t213;
t211 = pkin(1) * t212;
t206 = -pkin(8) * t209 + t202 * t211;
t180 = -t194 * pkin(2) - t206;
t199 = sin(qJ(3));
t201 = cos(qJ(3));
t185 = -t201 * t194 + t199 * t209;
t186 = t199 * t194 + t201 * t209;
t169 = t185 * pkin(3) - t186 * pkin(10) + t180;
t208 = t202 * t213;
t189 = -qJD(3) + t208;
t214 = pkin(8) * t208 + t200 * t211;
t181 = t194 * pkin(9) + t214;
t183 = (-pkin(2) * t202 - pkin(9) * t200 - pkin(1)) * t213;
t215 = t201 * t181 + t199 * t183;
t173 = -pkin(10) * t189 + t215;
t198 = sin(qJ(4));
t217 = t198 * t169 + t219 * t173;
t216 = -t199 * t181 + t201 * t183;
t210 = t202 * t218;
t184 = qJD(4) + t185;
t166 = t184 * qJ(5) + t217;
t172 = pkin(3) * t189 - t216;
t207 = t219 * t169 - t198 * t173;
t205 = qJD(5) - t207;
t175 = t219 * t186 - t198 * t189;
t204 = qJ(5) * t175 - t172;
t174 = t198 * t186 + t219 * t189;
t167 = pkin(4) * t174 - t204;
t165 = -t184 * pkin(4) + t205;
t164 = t220 * t174 + qJD(6) + t204;
t163 = qJ(6) * t174 + t166;
t162 = -t175 * qJ(6) + t220 * t184 + t205;
t1 = [t203 / 0.2e1, 0, 0, t200 ^ 2 * t218 / 0.2e1, t200 * t210, t194 * t209, t194 * t208, t194 ^ 2 / 0.2e1, pkin(1) * t210 + t206 * t194, -pkin(1) * t200 * t218 - t214 * t194, t186 ^ 2 / 0.2e1, -t186 * t185, -t186 * t189, t185 * t189, t189 ^ 2 / 0.2e1, t180 * t185 - t216 * t189, t180 * t186 + t215 * t189, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t184, -t174 * t184, t184 ^ 2 / 0.2e1, t172 * t174 + t207 * t184, t172 * t175 - t217 * t184, -t165 * t184 + t167 * t174, t165 * t175 - t166 * t174, t166 * t184 - t167 * t175, t166 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, -t162 * t184 - t164 * t174, t163 * t184 + t164 * t175, -t162 * t175 + t163 * t174, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1;];
T_reg  = t1;
