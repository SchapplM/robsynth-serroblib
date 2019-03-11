% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRR2
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:23
% EndTime: 2019-03-08 19:05:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (221->41), mult. (551->95), div. (0->0), fcn. (451->14), ass. (0->46)
t199 = cos(pkin(6)) * qJD(1) + qJD(2);
t205 = sin(pkin(7));
t208 = cos(pkin(7));
t207 = cos(pkin(13));
t206 = sin(pkin(6));
t229 = qJD(1) * t206;
t223 = t207 * t229;
t236 = t199 * t205 + t208 * t223;
t212 = sin(qJ(3));
t215 = cos(qJ(3));
t204 = sin(pkin(13));
t224 = t204 * t229;
t235 = -t212 * t224 + t236 * t215;
t216 = qJD(3) ^ 2;
t234 = t216 / 0.2e1;
t233 = cos(qJ(5));
t225 = t236 * t212 + t215 * t224;
t185 = qJD(3) * pkin(9) + t225;
t189 = t208 * t199 - t205 * t223;
t211 = sin(qJ(4));
t214 = cos(qJ(4));
t230 = t214 * t185 + t211 * t189;
t178 = qJD(4) * pkin(10) + t230;
t181 = (-pkin(4) * t214 - pkin(10) * t211 - pkin(3)) * qJD(3) - t235;
t210 = sin(qJ(5));
t231 = t233 * t178 + t210 * t181;
t228 = qJD(3) * t211;
t227 = t214 * qJD(3);
t226 = qJD(3) * qJD(4);
t222 = -t210 * t178 + t233 * t181;
t221 = -t211 * t185 + t214 * t189;
t200 = -qJD(5) + t227;
t177 = -qJD(4) * pkin(4) - t221;
t217 = qJD(1) ^ 2;
t213 = cos(qJ(6));
t209 = sin(qJ(6));
t197 = -qJD(6) + t200;
t194 = t210 * qJD(4) + t233 * t228;
t193 = -t233 * qJD(4) + t210 * t228;
t187 = -t209 * t193 + t213 * t194;
t186 = t213 * t193 + t209 * t194;
t184 = -qJD(3) * pkin(3) - t235;
t175 = t193 * pkin(5) + t177;
t174 = -t193 * pkin(11) + t231;
t173 = -t200 * pkin(5) - t194 * pkin(11) + t222;
t1 = [t217 / 0.2e1, t199 ^ 2 / 0.2e1 + (t204 ^ 2 / 0.2e1 + t207 ^ 2 / 0.2e1) * t217 * t206 ^ 2, t234, t235 * qJD(3), -t225 * qJD(3), t211 ^ 2 * t234, t214 * t216 * t211, t211 * t226, t214 * t226, qJD(4) ^ 2 / 0.2e1, t221 * qJD(4) - t184 * t227, -t230 * qJD(4) + t184 * t228, t194 ^ 2 / 0.2e1, -t194 * t193, -t194 * t200, t193 * t200, t200 ^ 2 / 0.2e1, t177 * t193 - t222 * t200, t177 * t194 + t231 * t200, t187 ^ 2 / 0.2e1, -t187 * t186, -t187 * t197, t186 * t197, t197 ^ 2 / 0.2e1 -(t213 * t173 - t209 * t174) * t197 + t175 * t186 (t209 * t173 + t213 * t174) * t197 + t175 * t187;];
T_reg  = t1;
