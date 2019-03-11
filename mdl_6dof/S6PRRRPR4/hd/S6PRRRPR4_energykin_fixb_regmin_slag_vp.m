% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR4
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:27
% EndTime: 2019-03-08 23:20:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (329->46), mult. (743->99), div. (0->0), fcn. (556->12), ass. (0->45)
t205 = qJD(2) ^ 2;
t218 = t205 / 0.2e1;
t217 = cos(qJ(4));
t199 = sin(qJ(4));
t200 = sin(qJ(3));
t212 = qJD(2) * t200;
t185 = t199 * qJD(3) + t217 * t212;
t203 = cos(qJ(3));
t211 = t203 * qJD(2);
t191 = -qJD(4) + t211;
t201 = sin(qJ(2));
t214 = qJD(1) * sin(pkin(6));
t186 = qJD(2) * pkin(8) + t201 * t214;
t213 = qJD(1) * cos(pkin(6));
t215 = t203 * t186 + t200 * t213;
t178 = qJD(3) * pkin(9) + t215;
t204 = cos(qJ(2));
t209 = t204 * t214;
t181 = -t209 + (-pkin(3) * t203 - pkin(9) * t200 - pkin(2)) * qJD(2);
t207 = -t199 * t178 + t217 * t181;
t166 = -t191 * pkin(4) - t185 * qJ(5) + t207;
t184 = -t217 * qJD(3) + t199 * t212;
t216 = t217 * t178 + t199 * t181;
t171 = -t184 * qJ(5) + t216;
t194 = sin(pkin(12));
t196 = cos(pkin(12));
t163 = t194 * t166 + t196 * t171;
t210 = qJD(2) * qJD(3);
t208 = qJD(2) * t214;
t162 = t196 * t166 - t194 * t171;
t206 = -t200 * t186 + t203 * t213;
t177 = -qJD(3) * pkin(3) - t206;
t172 = t184 * pkin(4) + qJD(5) + t177;
t202 = cos(qJ(6));
t198 = sin(qJ(6));
t188 = -qJD(6) + t191;
t187 = -qJD(2) * pkin(2) - t209;
t175 = -t194 * t184 + t196 * t185;
t174 = -t196 * t184 - t194 * t185;
t169 = t198 * t174 + t202 * t175;
t168 = -t202 * t174 + t198 * t175;
t167 = -t174 * pkin(5) + t172;
t161 = t174 * pkin(10) + t163;
t160 = -t191 * pkin(5) - t175 * pkin(10) + t162;
t1 = [qJD(1) ^ 2 / 0.2e1, t218, t204 * t208, -t201 * t208, t200 ^ 2 * t218, t200 * t205 * t203, t200 * t210, t203 * t210, qJD(3) ^ 2 / 0.2e1, t206 * qJD(3) - t187 * t211, -t215 * qJD(3) + t187 * t212, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t191, t184 * t191, t191 ^ 2 / 0.2e1, t177 * t184 - t207 * t191, t177 * t185 + t216 * t191, -t162 * t175 + t163 * t174, t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t188, t168 * t188, t188 ^ 2 / 0.2e1 -(t202 * t160 - t198 * t161) * t188 + t167 * t168 (t198 * t160 + t202 * t161) * t188 + t167 * t169;];
T_reg  = t1;
