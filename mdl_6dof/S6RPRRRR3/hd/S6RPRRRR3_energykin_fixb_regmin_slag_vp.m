% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:29
% EndTime: 2019-03-09 07:02:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (321->47), mult. (698->104), div. (0->0), fcn. (485->10), ass. (0->43)
t188 = qJD(1) ^ 2;
t201 = t188 / 0.2e1;
t200 = cos(qJ(5));
t183 = sin(qJ(4));
t186 = cos(qJ(4));
t184 = sin(qJ(3));
t196 = qJD(1) * t184;
t168 = t183 * qJD(3) + t186 * t196;
t187 = cos(qJ(3));
t195 = t187 * qJD(1);
t175 = -qJD(4) + t195;
t179 = sin(pkin(11));
t171 = (pkin(1) * t179 + pkin(7)) * qJD(1);
t197 = t184 * qJD(2) + t187 * t171;
t162 = qJD(3) * pkin(8) + t197;
t180 = cos(pkin(11));
t193 = -pkin(1) * t180 - pkin(2);
t163 = (-pkin(3) * t187 - pkin(8) * t184 + t193) * qJD(1);
t191 = -t183 * t162 + t186 * t163;
t150 = -t175 * pkin(4) - t168 * pkin(9) + t191;
t167 = -t186 * qJD(3) + t183 * t196;
t198 = t186 * t162 + t183 * t163;
t153 = -t167 * pkin(9) + t198;
t182 = sin(qJ(5));
t199 = t182 * t150 + t200 * t153;
t194 = qJD(1) * qJD(3);
t192 = t200 * t150 - t182 * t153;
t190 = t187 * qJD(2) - t184 * t171;
t173 = -qJD(5) + t175;
t161 = -qJD(3) * pkin(3) - t190;
t155 = t167 * pkin(4) + t161;
t185 = cos(qJ(6));
t181 = sin(qJ(6));
t172 = t193 * qJD(1);
t169 = -qJD(6) + t173;
t157 = -t182 * t167 + t200 * t168;
t156 = t200 * t167 + t182 * t168;
t151 = t156 * pkin(5) + t155;
t149 = -t181 * t156 + t185 * t157;
t148 = t185 * t156 + t181 * t157;
t145 = -t156 * pkin(10) + t199;
t144 = -t173 * pkin(5) - t157 * pkin(10) + t192;
t1 = [t201, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t179 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t188, t184 ^ 2 * t201, t184 * t188 * t187, t184 * t194, t187 * t194, qJD(3) ^ 2 / 0.2e1, t190 * qJD(3) - t172 * t195, -t197 * qJD(3) + t172 * t196, t168 ^ 2 / 0.2e1, -t168 * t167, -t168 * t175, t167 * t175, t175 ^ 2 / 0.2e1, t161 * t167 - t191 * t175, t161 * t168 + t198 * t175, t157 ^ 2 / 0.2e1, -t157 * t156, -t157 * t173, t156 * t173, t173 ^ 2 / 0.2e1, t155 * t156 - t192 * t173, t155 * t157 + t199 * t173, t149 ^ 2 / 0.2e1, -t149 * t148, -t149 * t169, t148 * t169, t169 ^ 2 / 0.2e1 -(t185 * t144 - t181 * t145) * t169 + t151 * t148 (t181 * t144 + t185 * t145) * t169 + t151 * t149;];
T_reg  = t1;
