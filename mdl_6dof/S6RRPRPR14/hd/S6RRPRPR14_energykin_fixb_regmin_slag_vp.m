% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:55
% EndTime: 2019-03-09 11:36:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (363->55), mult. (860->111), div. (0->0), fcn. (599->8), ass. (0->46)
t222 = -pkin(2) - pkin(9);
t221 = pkin(4) + pkin(10);
t198 = sin(pkin(6));
t206 = qJD(1) ^ 2;
t220 = t198 ^ 2 * t206;
t202 = sin(qJ(2));
t217 = qJD(1) * t198;
t212 = t202 * t217;
t192 = pkin(8) * t212;
t199 = cos(pkin(6));
t216 = t199 * qJD(1);
t196 = qJD(2) + t216;
t205 = cos(qJ(2));
t175 = qJD(3) + t192 + t222 * t196 + (-pkin(1) * t199 * t205 + pkin(3) * t198 * t202) * qJD(1);
t211 = -qJ(3) * t202 - pkin(1);
t181 = (t222 * t205 + t211) * t217;
t201 = sin(qJ(4));
t204 = cos(qJ(4));
t219 = t201 * t175 + t204 * t181;
t213 = t205 * t217;
t215 = pkin(1) * t216;
t218 = pkin(8) * t213 + t202 * t215;
t214 = t205 * t220;
t183 = -t196 * qJ(3) - t218;
t210 = t204 * t175 - t201 * t181;
t188 = t204 * t196 - t201 * t213;
t180 = pkin(3) * t213 - t183;
t190 = qJD(4) + t212;
t172 = -t190 * qJ(5) - t219;
t209 = qJD(5) - t210;
t208 = t205 * t215 - t192;
t207 = -t188 * qJ(5) + t180;
t203 = cos(qJ(6));
t200 = sin(qJ(6));
t187 = t201 * t196 + t204 * t213;
t186 = qJD(6) + t188;
t185 = (-pkin(2) * t205 + t211) * t217;
t182 = -t196 * pkin(2) + qJD(3) - t208;
t177 = t200 * t187 + t203 * t190;
t176 = -t203 * t187 + t200 * t190;
t173 = t187 * pkin(4) + t207;
t171 = -t190 * pkin(4) + t209;
t170 = t221 * t187 + t207;
t169 = -t187 * pkin(5) - t172;
t168 = t188 * pkin(5) - t221 * t190 + t209;
t1 = [t206 / 0.2e1, 0, 0, t202 ^ 2 * t220 / 0.2e1, t202 * t214, t196 * t212, t196 * t213, t196 ^ 2 / 0.2e1, pkin(1) * t214 + t208 * t196, -pkin(1) * t202 * t220 - t218 * t196 (t182 * t202 - t183 * t205) * t217, t182 * t196 + t185 * t213, -t183 * t196 - t185 * t212, t185 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t188 ^ 2 / 0.2e1, -t188 * t187, t188 * t190, -t187 * t190, t190 ^ 2 / 0.2e1, t180 * t187 + t210 * t190, t180 * t188 - t219 * t190, t171 * t188 + t172 * t187, t171 * t190 - t173 * t187, -t172 * t190 - t173 * t188, t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t186, -t176 * t186, t186 ^ 2 / 0.2e1 (t203 * t168 - t200 * t170) * t186 + t169 * t176 -(t200 * t168 + t203 * t170) * t186 + t169 * t177;];
T_reg  = t1;
