% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:35
% EndTime: 2019-03-09 08:52:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (549->53), mult. (1332->110), div. (0->0), fcn. (998->10), ass. (0->48)
t205 = qJD(1) ^ 2;
t216 = t205 / 0.2e1;
t215 = cos(qJ(5));
t214 = pkin(7) + qJ(3);
t204 = cos(qJ(2));
t213 = t204 * t205;
t197 = sin(pkin(10));
t199 = cos(pkin(10));
t210 = qJD(1) * t204;
t202 = sin(qJ(2));
t211 = qJD(1) * t202;
t187 = t197 * t211 - t199 * t210;
t188 = (t197 * t204 + t199 * t202) * qJD(1);
t193 = qJD(3) + (-pkin(2) * t204 - pkin(1)) * qJD(1);
t175 = t187 * pkin(3) - t188 * qJ(4) + t193;
t191 = qJD(2) * pkin(2) - t214 * t211;
t192 = t214 * t210;
t180 = t197 * t191 + t199 * t192;
t178 = qJD(2) * qJ(4) + t180;
t196 = sin(pkin(11));
t198 = cos(pkin(11));
t167 = t198 * t175 - t196 * t178;
t183 = t196 * qJD(2) + t198 * t188;
t161 = t187 * pkin(4) - t183 * pkin(8) + t167;
t168 = t196 * t175 + t198 * t178;
t182 = -t198 * qJD(2) + t196 * t188;
t166 = -t182 * pkin(8) + t168;
t201 = sin(qJ(5));
t212 = t201 * t161 + t215 * t166;
t209 = qJD(1) * qJD(2);
t208 = t202 * t209;
t207 = t204 * t209;
t206 = t215 * t161 - t201 * t166;
t179 = t199 * t191 - t197 * t192;
t186 = qJD(5) + t187;
t177 = -qJD(2) * pkin(3) + qJD(4) - t179;
t169 = t182 * pkin(4) + t177;
t203 = cos(qJ(6));
t200 = sin(qJ(6));
t184 = qJD(6) + t186;
t172 = -t201 * t182 + t215 * t183;
t171 = t215 * t182 + t201 * t183;
t164 = -t200 * t171 + t203 * t172;
t163 = t203 * t171 + t200 * t172;
t162 = t171 * pkin(5) + t169;
t158 = -t171 * pkin(9) + t212;
t157 = t186 * pkin(5) - t172 * pkin(9) + t206;
t1 = [t216, 0, 0, t202 ^ 2 * t216, t202 * t213, t208, t207, qJD(2) ^ 2 / 0.2e1, pkin(1) * t213 - pkin(7) * t208, -t205 * pkin(1) * t202 - pkin(7) * t207, -t179 * t188 - t180 * t187, t180 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t167 * t187 + t177 * t182, -t168 * t187 + t177 * t183, -t167 * t183 - t168 * t182, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1, t172 ^ 2 / 0.2e1, -t172 * t171, t172 * t186, -t171 * t186, t186 ^ 2 / 0.2e1, t169 * t171 + t206 * t186, t169 * t172 - t212 * t186, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t184, -t163 * t184, t184 ^ 2 / 0.2e1 (t203 * t157 - t200 * t158) * t184 + t162 * t163 -(t200 * t157 + t203 * t158) * t184 + t162 * t164;];
T_reg  = t1;
