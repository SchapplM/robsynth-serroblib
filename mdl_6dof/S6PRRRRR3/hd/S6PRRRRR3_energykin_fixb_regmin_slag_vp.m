% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:31
% EndTime: 2019-03-09 00:51:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (320->47), mult. (714->102), div. (0->0), fcn. (547->12), ass. (0->46)
t207 = qJD(2) ^ 2;
t222 = t207 / 0.2e1;
t221 = cos(qJ(5));
t200 = sin(qJ(4));
t204 = cos(qJ(4));
t201 = sin(qJ(3));
t215 = qJD(2) * t201;
t185 = t200 * qJD(3) + t204 * t215;
t205 = cos(qJ(3));
t214 = t205 * qJD(2);
t193 = -qJD(4) + t214;
t202 = sin(qJ(2));
t217 = qJD(1) * sin(pkin(6));
t186 = qJD(2) * pkin(8) + t202 * t217;
t216 = qJD(1) * cos(pkin(6));
t218 = t205 * t186 + t201 * t216;
t177 = qJD(3) * pkin(9) + t218;
t206 = cos(qJ(2));
t212 = t206 * t217;
t180 = -t212 + (-pkin(3) * t205 - pkin(9) * t201 - pkin(2)) * qJD(2);
t209 = -t200 * t177 + t204 * t180;
t165 = -t193 * pkin(4) - t185 * pkin(10) + t209;
t184 = -t204 * qJD(3) + t200 * t215;
t219 = t204 * t177 + t200 * t180;
t170 = -t184 * pkin(10) + t219;
t199 = sin(qJ(5));
t220 = t199 * t165 + t221 * t170;
t213 = qJD(2) * qJD(3);
t211 = qJD(2) * t217;
t210 = t221 * t165 - t199 * t170;
t208 = -t201 * t186 + t205 * t216;
t190 = -qJD(5) + t193;
t176 = -qJD(3) * pkin(3) - t208;
t171 = t184 * pkin(4) + t176;
t203 = cos(qJ(6));
t198 = sin(qJ(6));
t188 = -qJD(6) + t190;
t187 = -qJD(2) * pkin(2) - t212;
t174 = -t199 * t184 + t221 * t185;
t173 = t221 * t184 + t199 * t185;
t169 = -t198 * t173 + t203 * t174;
t168 = t203 * t173 + t198 * t174;
t166 = t173 * pkin(5) + t171;
t162 = -t173 * pkin(11) + t220;
t161 = -t190 * pkin(5) - t174 * pkin(11) + t210;
t1 = [qJD(1) ^ 2 / 0.2e1, t222, t206 * t211, -t202 * t211, t201 ^ 2 * t222, t201 * t207 * t205, t201 * t213, t205 * t213, qJD(3) ^ 2 / 0.2e1, t208 * qJD(3) - t187 * t214, -t218 * qJD(3) + t187 * t215, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t193, t184 * t193, t193 ^ 2 / 0.2e1, t176 * t184 - t209 * t193, t176 * t185 + t219 * t193, t174 ^ 2 / 0.2e1, -t174 * t173, -t174 * t190, t173 * t190, t190 ^ 2 / 0.2e1, t171 * t173 - t210 * t190, t171 * t174 + t220 * t190, t169 ^ 2 / 0.2e1, -t169 * t168, -t169 * t188, t168 * t188, t188 ^ 2 / 0.2e1 -(t203 * t161 - t198 * t162) * t188 + t166 * t168 (t198 * t161 + t203 * t162) * t188 + t166 * t169;];
T_reg  = t1;
