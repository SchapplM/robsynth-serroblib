% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:29
% EndTime: 2019-03-08 19:01:29
% DurationCPUTime: 0.16s
% Computational Cost: add. (227->42), mult. (568->97), div. (0->0), fcn. (468->14), ass. (0->46)
t196 = cos(pkin(6)) * qJD(1) + qJD(2);
t201 = sin(pkin(7));
t204 = cos(pkin(7));
t203 = cos(pkin(13));
t202 = sin(pkin(6));
t225 = qJD(1) * t202;
t219 = t203 * t225;
t231 = t196 * t201 + t204 * t219;
t208 = sin(qJ(3));
t212 = cos(qJ(3));
t200 = sin(pkin(13));
t220 = t200 * t225;
t230 = -t208 * t220 + t231 * t212;
t213 = qJD(3) ^ 2;
t229 = t213 / 0.2e1;
t221 = t231 * t208 + t212 * t220;
t182 = qJD(3) * pkin(9) + t221;
t187 = t204 * t196 - t201 * t219;
t211 = cos(qJ(4));
t186 = t211 * t187;
t207 = sin(qJ(4));
t177 = qJD(4) * pkin(4) + t186 + (-pkin(10) * qJD(3) - t182) * t207;
t223 = qJD(3) * t211;
t226 = t211 * t182 + t207 * t187;
t178 = pkin(10) * t223 + t226;
t206 = sin(qJ(5));
t210 = cos(qJ(5));
t227 = t206 * t177 + t210 * t178;
t224 = qJD(3) * t207;
t222 = qJD(3) * qJD(4);
t189 = t206 * t224 - t210 * t223;
t217 = t210 * t177 - t206 * t178;
t179 = (-pkin(4) * t211 - pkin(3)) * qJD(3) - t230;
t214 = qJD(1) ^ 2;
t209 = cos(qJ(6));
t205 = sin(qJ(6));
t199 = qJD(4) + qJD(5);
t190 = (t206 * t211 + t207 * t210) * qJD(3);
t188 = qJD(6) + t189;
t184 = t209 * t190 + t205 * t199;
t183 = t205 * t190 - t209 * t199;
t181 = -qJD(3) * pkin(3) - t230;
t175 = t189 * pkin(5) - t190 * pkin(11) + t179;
t173 = t199 * pkin(11) + t227;
t172 = -t199 * pkin(5) - t217;
t1 = [t214 / 0.2e1, t196 ^ 2 / 0.2e1 + (t200 ^ 2 / 0.2e1 + t203 ^ 2 / 0.2e1) * t214 * t202 ^ 2, t229, t230 * qJD(3), -t221 * qJD(3), t207 ^ 2 * t229, t207 * t213 * t211, t207 * t222, t211 * t222, qJD(4) ^ 2 / 0.2e1, -t181 * t223 + (-t207 * t182 + t186) * qJD(4), -t226 * qJD(4) + t181 * t224, t190 ^ 2 / 0.2e1, -t190 * t189, t190 * t199, -t189 * t199, t199 ^ 2 / 0.2e1, t179 * t189 + t217 * t199, t179 * t190 - t227 * t199, t184 ^ 2 / 0.2e1, -t184 * t183, t184 * t188, -t183 * t188, t188 ^ 2 / 0.2e1 (-t205 * t173 + t209 * t175) * t188 + t172 * t183 -(t209 * t173 + t205 * t175) * t188 + t172 * t184;];
T_reg  = t1;
