% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:18
% EndTime: 2019-03-09 13:25:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (518->52), mult. (1241->109), div. (0->0), fcn. (954->10), ass. (0->49)
t210 = qJD(1) ^ 2;
t223 = t210 / 0.2e1;
t222 = cos(qJ(5));
t221 = pkin(7) + qJ(3);
t209 = cos(qJ(2));
t220 = t209 * t210;
t201 = sin(pkin(11));
t202 = cos(pkin(11));
t206 = sin(qJ(2));
t193 = (t201 * t209 + t202 * t206) * qJD(1);
t205 = sin(qJ(4));
t208 = cos(qJ(4));
t187 = t205 * qJD(2) + t208 * t193;
t216 = qJD(1) * t209;
t217 = qJD(1) * t206;
t192 = -t201 * t217 + t202 * t216;
t191 = qJD(4) - t192;
t198 = qJD(3) + (-pkin(2) * t209 - pkin(1)) * qJD(1);
t179 = -t192 * pkin(3) - t193 * pkin(8) + t198;
t196 = qJD(2) * pkin(2) - t217 * t221;
t197 = t221 * t216;
t184 = t201 * t196 + t202 * t197;
t182 = qJD(2) * pkin(8) + t184;
t211 = t208 * t179 - t205 * t182;
t167 = t191 * pkin(4) - t187 * pkin(9) + t211;
t186 = -t208 * qJD(2) + t205 * t193;
t218 = t205 * t179 + t208 * t182;
t172 = -t186 * pkin(9) + t218;
t204 = sin(qJ(5));
t219 = t204 * t167 + t222 * t172;
t215 = qJD(1) * qJD(2);
t214 = t206 * t215;
t213 = t209 * t215;
t212 = t222 * t167 - t204 * t172;
t183 = t202 * t196 - t201 * t197;
t181 = -qJD(2) * pkin(3) - t183;
t189 = qJD(5) + t191;
t173 = t186 * pkin(4) + t181;
t207 = cos(qJ(6));
t203 = sin(qJ(6));
t188 = qJD(6) + t189;
t176 = -t204 * t186 + t187 * t222;
t175 = t186 * t222 + t204 * t187;
t171 = -t203 * t175 + t207 * t176;
t170 = t207 * t175 + t203 * t176;
t168 = t175 * pkin(5) + t173;
t164 = -t175 * pkin(10) + t219;
t163 = t189 * pkin(5) - t176 * pkin(10) + t212;
t1 = [t223, 0, 0, t206 ^ 2 * t223, t206 * t220, t214, t213, qJD(2) ^ 2 / 0.2e1, pkin(1) * t220 - pkin(7) * t214, -t210 * pkin(1) * t206 - pkin(7) * t213, -t183 * t193 + t184 * t192, t184 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1 + t198 ^ 2 / 0.2e1, t187 ^ 2 / 0.2e1, -t187 * t186, t187 * t191, -t186 * t191, t191 ^ 2 / 0.2e1, t181 * t186 + t191 * t211, t181 * t187 - t191 * t218, t176 ^ 2 / 0.2e1, -t176 * t175, t176 * t189, -t175 * t189, t189 ^ 2 / 0.2e1, t173 * t175 + t189 * t212, t173 * t176 - t189 * t219, t171 ^ 2 / 0.2e1, -t171 * t170, t171 * t188, -t170 * t188, t188 ^ 2 / 0.2e1 (t207 * t163 - t203 * t164) * t188 + t168 * t170 -(t203 * t163 + t207 * t164) * t188 + t168 * t171;];
T_reg  = t1;
