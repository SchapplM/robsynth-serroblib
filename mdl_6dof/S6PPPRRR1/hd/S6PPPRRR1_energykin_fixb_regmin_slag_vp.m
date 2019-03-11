% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:34
% EndTime: 2019-03-08 18:40:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (231->36), mult. (608->89), div. (0->0), fcn. (545->16), ass. (0->43)
t198 = sin(pkin(14));
t199 = sin(pkin(13));
t203 = cos(pkin(14));
t202 = sin(pkin(6));
t222 = qJD(1) * t202;
t204 = cos(pkin(13));
t206 = cos(pkin(7));
t224 = t204 * t206;
t193 = cos(pkin(6)) * qJD(1) + qJD(2);
t201 = sin(pkin(7));
t225 = t193 * t201;
t187 = t203 * t225 + (-t198 * t199 + t203 * t224) * t222;
t190 = -t201 * t204 * t222 + t206 * t193 + qJD(3);
t200 = sin(pkin(8));
t205 = cos(pkin(8));
t230 = t187 * t205 + t190 * t200;
t188 = t198 * t225 + (t198 * t224 + t199 * t203) * t222;
t209 = sin(qJ(4));
t212 = cos(qJ(4));
t229 = -t209 * t188 + t230 * t212;
t213 = qJD(4) ^ 2;
t228 = t213 / 0.2e1;
t218 = t212 * t188 + t230 * t209;
t181 = qJD(4) * pkin(10) + t218;
t183 = -t200 * t187 + t205 * t190;
t208 = sin(qJ(5));
t211 = cos(qJ(5));
t223 = t211 * t181 + t208 * t183;
t221 = qJD(4) * t208;
t220 = t211 * qJD(4);
t219 = qJD(4) * qJD(5);
t217 = -t208 * t181 + t211 * t183;
t214 = qJD(1) ^ 2;
t210 = cos(qJ(6));
t207 = sin(qJ(6));
t194 = -qJD(6) + t220;
t192 = t207 * qJD(5) + t210 * t221;
t191 = -t210 * qJD(5) + t207 * t221;
t180 = -qJD(4) * pkin(4) - t229;
t178 = (-pkin(5) * t211 - pkin(11) * t208 - pkin(4)) * qJD(4) - t229;
t177 = qJD(5) * pkin(11) + t223;
t176 = -qJD(5) * pkin(5) - t217;
t1 = [t214 / 0.2e1, t193 ^ 2 / 0.2e1 + (t199 ^ 2 / 0.2e1 + t204 ^ 2 / 0.2e1) * t214 * t202 ^ 2, t188 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1 + t190 ^ 2 / 0.2e1, t228, t229 * qJD(4), -t218 * qJD(4), t208 ^ 2 * t228, t208 * t213 * t211, t208 * t219, t211 * t219, qJD(5) ^ 2 / 0.2e1, t217 * qJD(5) - t180 * t220, -t223 * qJD(5) + t180 * t221, t192 ^ 2 / 0.2e1, -t191 * t192, -t194 * t192, t191 * t194, t194 ^ 2 / 0.2e1 -(-t207 * t177 + t210 * t178) * t194 + t176 * t191 (t210 * t177 + t207 * t178) * t194 + t176 * t192;];
T_reg  = t1;
