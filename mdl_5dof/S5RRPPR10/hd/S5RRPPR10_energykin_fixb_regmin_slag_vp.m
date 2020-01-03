% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:40
% EndTime: 2019-12-31 19:44:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (179->40), mult. (427->84), div. (0->0), fcn. (250->6), ass. (0->32)
t148 = qJD(1) ^ 2;
t157 = t148 / 0.2e1;
t147 = cos(qJ(2));
t156 = t147 * t148;
t145 = sin(qJ(2));
t131 = (-pkin(2) * t147 - qJ(3) * t145 - pkin(1)) * qJD(1);
t154 = t147 * qJD(1);
t137 = pkin(6) * t154 + qJD(2) * qJ(3);
t142 = sin(pkin(8));
t143 = cos(pkin(8));
t128 = t142 * t131 + t143 * t137;
t155 = qJD(1) * t145;
t153 = qJD(1) * qJD(2);
t152 = t145 * t153;
t151 = t147 * t153;
t127 = t143 * t131 - t142 * t137;
t150 = qJD(2) * pkin(2) - pkin(6) * t155 - qJD(3);
t124 = -qJ(4) * t154 + t128;
t123 = pkin(3) * t154 + qJD(4) - t127;
t133 = t142 * qJD(2) + t143 * t155;
t149 = t133 * qJ(4) + t150;
t146 = cos(qJ(5));
t144 = sin(qJ(5));
t138 = qJD(5) + t154;
t132 = -t143 * qJD(2) + t142 * t155;
t126 = t144 * t132 + t146 * t133;
t125 = -t146 * t132 + t144 * t133;
t122 = t132 * pkin(3) - t149;
t121 = t132 * pkin(7) + t124;
t120 = (-pkin(3) - pkin(4)) * t132 + t149;
t119 = pkin(4) * t154 - t133 * pkin(7) + t123;
t1 = [t157, 0, 0, t145 ^ 2 * t157, t145 * t156, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t156 - pkin(6) * t152, -t148 * pkin(1) * t145 - pkin(6) * t151, -t127 * t154 - t132 * t150, t128 * t154 - t133 * t150, -t127 * t133 - t128 * t132, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1, t122 * t132 + t123 * t154, t123 * t133 - t124 * t132, -t122 * t133 - t124 * t154, t124 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * t138, -t125 * t138, t138 ^ 2 / 0.2e1, (t146 * t119 - t144 * t121) * t138 + t120 * t125, -(t144 * t119 + t146 * t121) * t138 + t120 * t126;];
T_reg = t1;
