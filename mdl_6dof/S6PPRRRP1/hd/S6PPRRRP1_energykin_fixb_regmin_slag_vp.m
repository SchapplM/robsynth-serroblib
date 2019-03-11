% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:24
% EndTime: 2019-03-08 18:54:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (182->36), mult. (466->84), div. (0->0), fcn. (372->12), ass. (0->41)
t181 = cos(pkin(6)) * qJD(1) + qJD(2);
t187 = sin(pkin(7));
t190 = cos(pkin(7));
t189 = cos(pkin(12));
t188 = sin(pkin(6));
t209 = qJD(1) * t188;
t203 = t189 * t209;
t216 = t181 * t187 + t190 * t203;
t193 = sin(qJ(3));
t195 = cos(qJ(3));
t186 = sin(pkin(12));
t204 = t186 * t209;
t215 = -t193 * t204 + t216 * t195;
t196 = qJD(3) ^ 2;
t214 = t196 / 0.2e1;
t213 = cos(qJ(5));
t205 = t216 * t193 + t195 * t204;
t172 = qJD(3) * pkin(9) + t205;
t174 = t190 * t181 - t187 * t203;
t192 = sin(qJ(4));
t194 = cos(qJ(4));
t210 = t194 * t172 + t192 * t174;
t165 = qJD(4) * pkin(10) + t210;
t168 = (-pkin(4) * t194 - pkin(10) * t192 - pkin(3)) * qJD(3) - t215;
t191 = sin(qJ(5));
t211 = t213 * t165 + t191 * t168;
t208 = qJD(3) * t192;
t207 = t194 * qJD(3);
t206 = qJD(3) * qJD(4);
t202 = -t191 * t165 + t213 * t168;
t201 = -t192 * t172 + t194 * t174;
t164 = -qJD(4) * pkin(4) - t201;
t197 = qJD(1) ^ 2;
t182 = -qJD(5) + t207;
t178 = t191 * qJD(4) + t213 * t208;
t177 = -t213 * qJD(4) + t191 * t208;
t171 = -qJD(3) * pkin(3) - t215;
t162 = t177 * pkin(5) + qJD(6) + t164;
t161 = -t177 * qJ(6) + t211;
t160 = -t182 * pkin(5) - t178 * qJ(6) + t202;
t1 = [t197 / 0.2e1, t181 ^ 2 / 0.2e1 + (t186 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1) * t197 * t188 ^ 2, t214, t215 * qJD(3), -t205 * qJD(3), t192 ^ 2 * t214, t192 * t196 * t194, t192 * t206, t194 * t206, qJD(4) ^ 2 / 0.2e1, t201 * qJD(4) - t171 * t207, -t210 * qJD(4) + t171 * t208, t178 ^ 2 / 0.2e1, -t178 * t177, -t178 * t182, t177 * t182, t182 ^ 2 / 0.2e1, t164 * t177 - t202 * t182, t164 * t178 + t211 * t182, -t160 * t178 - t161 * t177, t161 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1;];
T_reg  = t1;
