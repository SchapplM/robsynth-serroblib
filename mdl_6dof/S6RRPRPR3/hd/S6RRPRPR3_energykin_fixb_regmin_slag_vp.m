% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:06
% EndTime: 2019-03-09 10:19:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (535->51), mult. (1279->106), div. (0->0), fcn. (956->10), ass. (0->48)
t207 = qJD(1) ^ 2;
t218 = t207 / 0.2e1;
t217 = cos(qJ(4));
t216 = pkin(7) + qJ(3);
t206 = cos(qJ(2));
t215 = t206 * t207;
t199 = sin(pkin(10));
t201 = cos(pkin(10));
t204 = sin(qJ(2));
t190 = (t199 * t206 + t201 * t204) * qJD(1);
t203 = sin(qJ(4));
t185 = t203 * qJD(2) + t217 * t190;
t212 = qJD(1) * t206;
t213 = qJD(1) * t204;
t189 = -t199 * t213 + t201 * t212;
t188 = qJD(4) - t189;
t195 = qJD(3) + (-pkin(2) * t206 - pkin(1)) * qJD(1);
t178 = -t189 * pkin(3) - t190 * pkin(8) + t195;
t193 = qJD(2) * pkin(2) - t216 * t213;
t194 = t216 * t212;
t183 = t199 * t193 + t201 * t194;
t181 = qJD(2) * pkin(8) + t183;
t208 = t217 * t178 - t203 * t181;
t166 = t188 * pkin(4) - t185 * qJ(5) + t208;
t184 = -t217 * qJD(2) + t203 * t190;
t214 = t203 * t178 + t217 * t181;
t171 = -t184 * qJ(5) + t214;
t198 = sin(pkin(11));
t200 = cos(pkin(11));
t163 = t198 * t166 + t200 * t171;
t211 = qJD(1) * qJD(2);
t210 = t204 * t211;
t209 = t206 * t211;
t162 = t200 * t166 - t198 * t171;
t182 = t201 * t193 - t199 * t194;
t180 = -qJD(2) * pkin(3) - t182;
t172 = t184 * pkin(4) + qJD(5) + t180;
t205 = cos(qJ(6));
t202 = sin(qJ(6));
t186 = qJD(6) + t188;
t175 = -t198 * t184 + t200 * t185;
t174 = -t200 * t184 - t198 * t185;
t169 = t202 * t174 + t205 * t175;
t168 = -t205 * t174 + t202 * t175;
t167 = -t174 * pkin(5) + t172;
t161 = t174 * pkin(9) + t163;
t160 = t188 * pkin(5) - t175 * pkin(9) + t162;
t1 = [t218, 0, 0, t204 ^ 2 * t218, t204 * t215, t210, t209, qJD(2) ^ 2 / 0.2e1, pkin(1) * t215 - pkin(7) * t210, -t207 * pkin(1) * t204 - pkin(7) * t209, -t182 * t190 + t183 * t189, t183 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t195 ^ 2 / 0.2e1, t185 ^ 2 / 0.2e1, -t185 * t184, t185 * t188, -t184 * t188, t188 ^ 2 / 0.2e1, t180 * t184 + t208 * t188, t180 * t185 - t214 * t188, -t162 * t175 + t163 * t174, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t186, -t168 * t186, t186 ^ 2 / 0.2e1 (t205 * t160 - t202 * t161) * t186 + t167 * t168 -(t202 * t160 + t205 * t161) * t186 + t167 * t169;];
T_reg  = t1;
