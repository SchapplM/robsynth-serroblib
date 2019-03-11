% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:46
% EndTime: 2019-03-09 21:49:46
% DurationCPUTime: 0.14s
% Computational Cost: add. (637->52), mult. (1443->105), div. (0->0), fcn. (1095->8), ass. (0->44)
t218 = pkin(4) + qJ(6);
t194 = sin(pkin(6));
t202 = qJD(1) ^ 2;
t217 = t194 ^ 2 * t202;
t212 = cos(pkin(6)) * qJD(1);
t192 = qJD(2) + t212;
t201 = cos(qJ(2));
t198 = sin(qJ(2));
t213 = qJD(1) * t194;
t209 = t198 * t213;
t211 = pkin(1) * t212;
t204 = -pkin(8) * t209 + t201 * t211;
t179 = -t192 * pkin(2) - t204;
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t183 = -t200 * t192 + t197 * t209;
t184 = t197 * t192 + t200 * t209;
t169 = t183 * pkin(3) - t184 * pkin(10) + t179;
t208 = t201 * t213;
t187 = -qJD(3) + t208;
t214 = pkin(8) * t208 + t198 * t211;
t180 = t192 * pkin(9) + t214;
t181 = (-pkin(2) * t201 - pkin(9) * t198 - pkin(1)) * t213;
t215 = t200 * t180 + t197 * t181;
t173 = -t187 * pkin(10) + t215;
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t216 = t196 * t169 + t199 * t173;
t210 = t201 * t217;
t207 = t199 * t169 - t196 * t173;
t206 = -t197 * t180 + t200 * t181;
t182 = qJD(4) + t183;
t166 = -qJ(5) * t182 - t216;
t205 = qJD(5) - t207;
t172 = t187 * pkin(3) - t206;
t175 = t199 * t184 - t196 * t187;
t203 = -t175 * qJ(5) + t172;
t174 = t196 * t184 + t199 * t187;
t167 = t174 * pkin(4) + t203;
t165 = -t182 * pkin(4) + t205;
t164 = t218 * t174 + t203;
t163 = -t174 * pkin(5) + qJD(6) - t166;
t162 = t175 * pkin(5) - t218 * t182 + t205;
t1 = [t202 / 0.2e1, 0, 0, t198 ^ 2 * t217 / 0.2e1, t198 * t210, t192 * t209, t192 * t208, t192 ^ 2 / 0.2e1, pkin(1) * t210 + t204 * t192, -pkin(1) * t198 * t217 - t214 * t192, t184 ^ 2 / 0.2e1, -t184 * t183, -t184 * t187, t183 * t187, t187 ^ 2 / 0.2e1, t179 * t183 - t206 * t187, t179 * t184 + t215 * t187, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t182, -t174 * t182, t182 ^ 2 / 0.2e1, t172 * t174 + t207 * t182, t172 * t175 - t216 * t182, t165 * t175 + t166 * t174, t165 * t182 - t167 * t174, -t166 * t182 - t167 * t175, t167 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t162 * t175 - t163 * t174, t163 * t182 - t164 * t175, -t162 * t182 + t164 * t174, t164 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1;];
T_reg  = t1;
