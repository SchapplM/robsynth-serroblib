% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:33
% EndTime: 2019-03-09 11:29:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (547->58), mult. (1269->119), div. (0->0), fcn. (928->10), ass. (0->49)
t230 = -pkin(2) - pkin(9);
t229 = cos(pkin(11));
t208 = sin(pkin(6));
t216 = qJD(1) ^ 2;
t228 = t208 ^ 2 * t216;
t212 = sin(qJ(2));
t225 = qJD(1) * t208;
t220 = t212 * t225;
t199 = qJD(4) + t220;
t201 = pkin(8) * t220;
t209 = cos(pkin(6));
t224 = t209 * qJD(1);
t205 = qJD(2) + t224;
t215 = cos(qJ(2));
t184 = qJD(3) + t201 + t230 * t205 + (-pkin(1) * t209 * t215 + pkin(3) * t208 * t212) * qJD(1);
t219 = -qJ(3) * t212 - pkin(1);
t191 = (t230 * t215 + t219) * t225;
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t227 = t211 * t184 + t214 * t191;
t177 = t199 * qJ(5) + t227;
t221 = t215 * t225;
t223 = pkin(1) * t224;
t226 = pkin(8) * t221 + t212 * t223;
t193 = -t205 * qJ(3) - t226;
t190 = pkin(3) * t221 - t193;
t196 = t211 * t205 + t214 * t221;
t197 = t214 * t205 - t211 * t221;
t182 = t196 * pkin(4) - t197 * qJ(5) + t190;
t207 = sin(pkin(11));
t173 = t229 * t177 + t207 * t182;
t222 = t215 * t228;
t172 = -t207 * t177 + t229 * t182;
t218 = t214 * t184 - t211 * t191;
t217 = t215 * t223 - t201;
t176 = -t199 * pkin(4) + qJD(5) - t218;
t213 = cos(qJ(6));
t210 = sin(qJ(6));
t195 = qJD(6) + t196;
t194 = (-pkin(2) * t215 + t219) * t225;
t192 = -t205 * pkin(2) + qJD(3) - t217;
t187 = t229 * t197 + t207 * t199;
t186 = t207 * t197 - t229 * t199;
t179 = -t210 * t186 + t213 * t187;
t178 = t213 * t186 + t210 * t187;
t174 = t186 * pkin(5) + t176;
t171 = -t186 * pkin(10) + t173;
t170 = t196 * pkin(5) - t187 * pkin(10) + t172;
t1 = [t216 / 0.2e1, 0, 0, t212 ^ 2 * t228 / 0.2e1, t212 * t222, t205 * t220, t205 * t221, t205 ^ 2 / 0.2e1, pkin(1) * t222 + t217 * t205, -pkin(1) * t212 * t228 - t226 * t205 (t192 * t212 - t193 * t215) * t225, t192 * t205 + t194 * t221, -t193 * t205 - t194 * t220, t194 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1, t197 ^ 2 / 0.2e1, -t197 * t196, t197 * t199, -t196 * t199, t199 ^ 2 / 0.2e1, t190 * t196 + t218 * t199, t190 * t197 - t227 * t199, t172 * t196 + t176 * t186, -t173 * t196 + t176 * t187, -t172 * t187 - t173 * t186, t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t195, -t178 * t195, t195 ^ 2 / 0.2e1 (t213 * t170 - t210 * t171) * t195 + t174 * t178 -(t210 * t170 + t213 * t171) * t195 + t174 * t179;];
T_reg  = t1;
