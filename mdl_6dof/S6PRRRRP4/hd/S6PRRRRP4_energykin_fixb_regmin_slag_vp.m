% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:46
% EndTime: 2019-03-09 00:17:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (335->44), mult. (721->95), div. (0->0), fcn. (520->10), ass. (0->41)
t195 = qJD(2) ^ 2;
t210 = t195 / 0.2e1;
t209 = cos(qJ(4));
t189 = sin(qJ(4));
t190 = sin(qJ(3));
t203 = qJD(2) * t190;
t177 = t189 * qJD(3) + t209 * t203;
t193 = cos(qJ(3));
t202 = t193 * qJD(2);
t183 = -qJD(4) + t202;
t191 = sin(qJ(2));
t205 = qJD(1) * sin(pkin(6));
t178 = qJD(2) * pkin(8) + t191 * t205;
t204 = qJD(1) * cos(pkin(6));
t206 = t193 * t178 + t190 * t204;
t169 = qJD(3) * pkin(9) + t206;
t194 = cos(qJ(2));
t200 = t194 * t205;
t172 = -t200 + (-pkin(3) * t193 - pkin(9) * t190 - pkin(2)) * qJD(2);
t198 = -t189 * t169 + t209 * t172;
t161 = -t183 * pkin(4) - t177 * pkin(10) + t198;
t176 = -t209 * qJD(3) + t189 * t203;
t207 = t209 * t169 + t189 * t172;
t163 = -t176 * pkin(10) + t207;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t208 = t188 * t161 + t192 * t163;
t201 = qJD(2) * qJD(3);
t199 = qJD(2) * t205;
t197 = -t190 * t178 + t193 * t204;
t196 = t192 * t161 - t188 * t163;
t168 = -qJD(3) * pkin(3) - t197;
t164 = t176 * pkin(4) + t168;
t180 = -qJD(5) + t183;
t179 = -qJD(2) * pkin(2) - t200;
t166 = -t188 * t176 + t192 * t177;
t165 = t192 * t176 + t188 * t177;
t159 = t165 * pkin(5) - t166 * qJ(6) + t164;
t158 = -t180 * qJ(6) + t208;
t157 = t180 * pkin(5) + qJD(6) - t196;
t1 = [qJD(1) ^ 2 / 0.2e1, t210, t194 * t199, -t191 * t199, t190 ^ 2 * t210, t190 * t195 * t193, t190 * t201, t193 * t201, qJD(3) ^ 2 / 0.2e1, t197 * qJD(3) - t179 * t202, -t206 * qJD(3) + t179 * t203, t177 ^ 2 / 0.2e1, -t177 * t176, -t177 * t183, t176 * t183, t183 ^ 2 / 0.2e1, t168 * t176 - t198 * t183, t168 * t177 + t207 * t183, t166 ^ 2 / 0.2e1, -t166 * t165, -t166 * t180, t165 * t180, t180 ^ 2 / 0.2e1, t164 * t165 - t196 * t180, t164 * t166 + t208 * t180, t157 * t180 + t159 * t165, t157 * t166 - t158 * t165, -t158 * t180 - t159 * t166, t158 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1;];
T_reg  = t1;
