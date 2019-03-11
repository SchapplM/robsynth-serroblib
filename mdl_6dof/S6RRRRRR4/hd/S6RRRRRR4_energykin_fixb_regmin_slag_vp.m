% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:51:24
% EndTime: 2019-03-10 03:51:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (626->54), mult. (1366->113), div. (0->0), fcn. (1060->10), ass. (0->49)
t207 = qJD(1) ^ 2;
t222 = t207 / 0.2e1;
t221 = cos(qJ(3));
t220 = cos(qJ(5));
t206 = cos(qJ(2));
t219 = t206 * t207;
t202 = sin(qJ(3));
t203 = sin(qJ(2));
t215 = qJD(1) * t203;
t184 = -t221 * qJD(2) + t202 * t215;
t185 = t202 * qJD(2) + t221 * t215;
t201 = sin(qJ(4));
t205 = cos(qJ(4));
t178 = -t201 * t184 + t205 * t185;
t214 = t206 * qJD(1);
t195 = -qJD(3) + t214;
t193 = -qJD(4) + t195;
t183 = (-pkin(2) * t206 - pkin(8) * t203 - pkin(1)) * qJD(1);
t192 = pkin(7) * t214 + qJD(2) * pkin(8);
t208 = t221 * t183 - t202 * t192;
t173 = -t195 * pkin(3) - t185 * pkin(9) + t208;
t216 = t202 * t183 + t221 * t192;
t175 = -t184 * pkin(9) + t216;
t209 = t205 * t173 - t201 * t175;
t163 = -t193 * pkin(4) - t178 * pkin(10) + t209;
t177 = t205 * t184 + t201 * t185;
t217 = t201 * t173 + t205 * t175;
t165 = -t177 * pkin(10) + t217;
t200 = sin(qJ(5));
t218 = t200 * t163 + t220 * t165;
t213 = qJD(1) * qJD(2);
t212 = t203 * t213;
t211 = t206 * t213;
t191 = -qJD(2) * pkin(2) + pkin(7) * t215;
t210 = t220 * t163 - t200 * t165;
t179 = t184 * pkin(3) + t191;
t189 = -qJD(5) + t193;
t170 = t177 * pkin(4) + t179;
t204 = cos(qJ(6));
t199 = sin(qJ(6));
t187 = -qJD(6) + t189;
t169 = -t200 * t177 + t220 * t178;
t168 = t220 * t177 + t200 * t178;
t166 = t168 * pkin(5) + t170;
t162 = -t199 * t168 + t204 * t169;
t161 = t204 * t168 + t199 * t169;
t158 = -t168 * pkin(11) + t218;
t157 = -t189 * pkin(5) - t169 * pkin(11) + t210;
t1 = [t222, 0, 0, t203 ^ 2 * t222, t203 * t219, t212, t211, qJD(2) ^ 2 / 0.2e1, pkin(1) * t219 - pkin(7) * t212, -t207 * pkin(1) * t203 - pkin(7) * t211, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t195, t184 * t195, t195 ^ 2 / 0.2e1, t191 * t184 - t208 * t195, t191 * t185 + t216 * t195, t178 ^ 2 / 0.2e1, -t178 * t177, -t178 * t193, t177 * t193, t193 ^ 2 / 0.2e1, t179 * t177 - t209 * t193, t179 * t178 + t217 * t193, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t189, t168 * t189, t189 ^ 2 / 0.2e1, t170 * t168 - t210 * t189, t170 * t169 + t218 * t189, t162 ^ 2 / 0.2e1, -t162 * t161, -t162 * t187, t161 * t187, t187 ^ 2 / 0.2e1 -(t204 * t157 - t199 * t158) * t187 + t166 * t161 (t199 * t157 + t204 * t158) * t187 + t166 * t162;];
T_reg  = t1;
