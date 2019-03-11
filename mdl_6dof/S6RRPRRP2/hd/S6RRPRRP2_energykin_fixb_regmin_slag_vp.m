% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:11
% EndTime: 2019-03-09 11:45:12
% DurationCPUTime: 0.15s
% Computational Cost: add. (560->49), mult. (1350->102), div. (0->0), fcn. (1002->8), ass. (0->42)
t201 = qJD(1) * (pkin(7) + qJ(3));
t190 = qJD(1) ^ 2;
t200 = t190 / 0.2e1;
t189 = cos(qJ(2));
t198 = t189 * t190;
t181 = qJD(2) + qJD(4);
t186 = sin(qJ(2));
t177 = qJD(2) * pkin(2) - t186 * t201;
t178 = t189 * t201;
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t168 = t183 * t177 - t182 * t178;
t175 = (t182 * t189 + t183 * t186) * qJD(1);
t161 = qJD(2) * pkin(3) - t175 * pkin(8) + t168;
t169 = t182 * t177 + t183 * t178;
t174 = (-t182 * t186 + t183 * t189) * qJD(1);
t162 = t174 * pkin(8) + t169;
t185 = sin(qJ(4));
t188 = cos(qJ(4));
t196 = t185 * t161 + t188 * t162;
t155 = t181 * pkin(9) + t196;
t166 = -t188 * t174 + t185 * t175;
t167 = t185 * t174 + t188 * t175;
t179 = qJD(3) + (-pkin(2) * t189 - pkin(1)) * qJD(1);
t170 = -t174 * pkin(3) + t179;
t157 = t166 * pkin(4) - t167 * pkin(9) + t170;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t197 = t187 * t155 + t184 * t157;
t195 = qJD(1) * qJD(2);
t194 = t186 * t195;
t193 = t189 * t195;
t192 = t188 * t161 - t185 * t162;
t191 = -t184 * t155 + t187 * t157;
t154 = -t181 * pkin(4) - t192;
t165 = qJD(5) + t166;
t164 = t187 * t167 + t184 * t181;
t163 = t184 * t167 - t187 * t181;
t152 = t163 * pkin(5) - t164 * qJ(6) + t154;
t151 = t165 * qJ(6) + t197;
t150 = -t165 * pkin(5) + qJD(6) - t191;
t1 = [t200, 0, 0, t186 ^ 2 * t200, t186 * t198, t194, t193, qJD(2) ^ 2 / 0.2e1, pkin(1) * t198 - pkin(7) * t194, -t190 * pkin(1) * t186 - pkin(7) * t193, -t168 * t175 + t169 * t174, t169 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t181, -t166 * t181, t181 ^ 2 / 0.2e1, t170 * t166 + t192 * t181, t170 * t167 - t196 * t181, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t165, -t163 * t165, t165 ^ 2 / 0.2e1, t154 * t163 + t191 * t165, t154 * t164 - t197 * t165, -t150 * t165 + t152 * t163, t150 * t164 - t151 * t163, t151 * t165 - t152 * t164, t151 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1;];
T_reg  = t1;
