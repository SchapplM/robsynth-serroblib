% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRP2
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
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:08
% EndTime: 2019-03-08 18:58:08
% DurationCPUTime: 0.15s
% Computational Cost: add. (246->38), mult. (600->88), div. (0->0), fcn. (483->12), ass. (0->41)
t188 = cos(pkin(6)) * qJD(1) + qJD(2);
t194 = sin(pkin(7));
t197 = cos(pkin(7));
t196 = cos(pkin(12));
t195 = sin(pkin(6));
t217 = qJD(1) * t195;
t211 = t196 * t217;
t223 = t188 * t194 + t197 * t211;
t200 = sin(qJ(3));
t203 = cos(qJ(3));
t193 = sin(pkin(12));
t212 = t193 * t217;
t222 = -t200 * t212 + t223 * t203;
t204 = qJD(3) ^ 2;
t221 = t204 / 0.2e1;
t213 = t223 * t200 + t203 * t212;
t179 = qJD(3) * pkin(9) + t213;
t181 = t197 * t188 - t194 * t211;
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t218 = t202 * t179 + t199 * t181;
t173 = qJD(4) * pkin(10) + t218;
t175 = (-pkin(4) * t202 - pkin(10) * t199 - pkin(3)) * qJD(3) - t222;
t198 = sin(qJ(5));
t201 = cos(qJ(5));
t219 = t201 * t173 + t198 * t175;
t216 = qJD(3) * t199;
t215 = t202 * qJD(3);
t214 = qJD(3) * qJD(4);
t210 = -t199 * t179 + t202 * t181;
t208 = -t198 * t173 + t201 * t175;
t172 = -qJD(4) * pkin(4) - t210;
t205 = qJD(1) ^ 2;
t189 = -qJD(5) + t215;
t185 = t198 * qJD(4) + t201 * t216;
t184 = -t201 * qJD(4) + t198 * t216;
t178 = -qJD(3) * pkin(3) - t222;
t170 = t184 * pkin(5) - t185 * qJ(6) + t172;
t169 = -t189 * qJ(6) + t219;
t168 = t189 * pkin(5) + qJD(6) - t208;
t1 = [t205 / 0.2e1, t188 ^ 2 / 0.2e1 + (t193 ^ 2 / 0.2e1 + t196 ^ 2 / 0.2e1) * t205 * t195 ^ 2, t221, t222 * qJD(3), -t213 * qJD(3), t199 ^ 2 * t221, t199 * t204 * t202, t199 * t214, t202 * t214, qJD(4) ^ 2 / 0.2e1, t210 * qJD(4) - t178 * t215, -t218 * qJD(4) + t178 * t216, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t189, t184 * t189, t189 ^ 2 / 0.2e1, t172 * t184 - t208 * t189, t172 * t185 + t219 * t189, t168 * t189 + t170 * t184, t168 * t185 - t169 * t184, -t169 * t189 - t170 * t185, t169 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1;];
T_reg  = t1;
