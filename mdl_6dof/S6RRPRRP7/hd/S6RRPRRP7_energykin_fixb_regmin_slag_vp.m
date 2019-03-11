% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:05
% EndTime: 2019-03-09 12:20:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (336->49), mult. (678->102), div. (0->0), fcn. (407->6), ass. (0->38)
t182 = qJD(1) ^ 2;
t194 = t182 / 0.2e1;
t181 = cos(qJ(2));
t193 = t181 * t182;
t189 = qJD(1) * t181;
t178 = sin(qJ(2));
t190 = qJD(1) * t178;
t163 = -qJD(1) * pkin(1) - pkin(2) * t189 - qJ(3) * t190;
t155 = pkin(3) * t189 - t163;
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t160 = (t177 * t178 + t180 * t181) * qJD(1);
t161 = (-t177 * t181 + t178 * t180) * qJD(1);
t148 = t160 * pkin(4) - t161 * pkin(9) + t155;
t172 = qJD(2) - qJD(4);
t188 = pkin(7) * t190 + qJD(3);
t156 = -pkin(8) * t190 + (-pkin(2) - pkin(3)) * qJD(2) + t188;
t165 = pkin(7) * t189 + qJD(2) * qJ(3);
t162 = -pkin(8) * t189 + t165;
t191 = t177 * t156 + t180 * t162;
t151 = -t172 * pkin(9) + t191;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t192 = t176 * t148 + t179 * t151;
t187 = qJD(1) * qJD(2);
t186 = t178 * t187;
t185 = t181 * t187;
t184 = t180 * t156 - t177 * t162;
t183 = t179 * t148 - t176 * t151;
t150 = t172 * pkin(4) - t184;
t164 = -qJD(2) * pkin(2) + t188;
t159 = qJD(5) + t160;
t153 = t179 * t161 - t176 * t172;
t152 = t176 * t161 + t179 * t172;
t146 = t152 * pkin(5) - t153 * qJ(6) + t150;
t145 = t159 * qJ(6) + t192;
t144 = -t159 * pkin(5) + qJD(6) - t183;
t1 = [t194, 0, 0, t178 ^ 2 * t194, t178 * t193, t186, t185, qJD(2) ^ 2 / 0.2e1, pkin(1) * t193 - pkin(7) * t186, -t182 * pkin(1) * t178 - pkin(7) * t185, -t164 * qJD(2) - t163 * t189 (t164 * t178 + t165 * t181) * qJD(1), t165 * qJD(2) - t163 * t190, t165 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t161 ^ 2 / 0.2e1, -t161 * t160, -t161 * t172, t160 * t172, t172 ^ 2 / 0.2e1, t155 * t160 - t184 * t172, t155 * t161 + t191 * t172, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t159, -t152 * t159, t159 ^ 2 / 0.2e1, t150 * t152 + t183 * t159, t150 * t153 - t192 * t159, -t144 * t159 + t146 * t152, t144 * t153 - t145 * t152, t145 * t159 - t146 * t153, t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
