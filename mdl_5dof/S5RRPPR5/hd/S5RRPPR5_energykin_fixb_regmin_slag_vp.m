% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR5
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
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:57
% EndTime: 2019-12-31 19:29:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (163->37), mult. (420->79), div. (0->0), fcn. (261->6), ass. (0->35)
t138 = sin(pkin(8));
t139 = cos(pkin(8));
t141 = sin(qJ(2));
t143 = cos(qJ(2));
t128 = (t138 * t143 + t139 * t141) * qJD(1);
t133 = -(pkin(2) * t143 + pkin(1)) * qJD(1) + qJD(3);
t157 = -t128 * qJ(4) + t133;
t156 = -pkin(3) - pkin(4);
t144 = qJD(1) ^ 2;
t155 = t144 / 0.2e1;
t154 = pkin(6) + qJ(3);
t152 = t143 * t144;
t151 = qJD(1) * t141;
t131 = qJD(2) * pkin(2) - t154 * t151;
t150 = qJD(1) * t143;
t132 = t154 * t150;
t124 = t138 * t131 + t139 * t132;
t149 = qJD(1) * qJD(2);
t122 = qJD(2) * qJ(4) + t124;
t147 = t141 * t149;
t146 = t143 * t149;
t123 = t139 * t131 - t138 * t132;
t145 = qJD(4) - t123;
t142 = cos(qJ(5));
t140 = sin(qJ(5));
t135 = qJD(2) - qJD(5);
t127 = t138 * t151 - t139 * t150;
t121 = t140 * t127 + t142 * t128;
t120 = -t142 * t127 + t140 * t128;
t119 = -qJD(2) * pkin(3) + t145;
t118 = t127 * pkin(3) + t157;
t117 = t127 * pkin(7) + t122;
t116 = -t128 * pkin(7) + t156 * qJD(2) + t145;
t115 = t156 * t127 - t157;
t1 = [t155, 0, 0, t141 ^ 2 * t155, t141 * t152, t147, t146, qJD(2) ^ 2 / 0.2e1, pkin(1) * t152 - pkin(6) * t147, -t144 * pkin(1) * t141 - pkin(6) * t146, -t123 * t128 - t124 * t127, t124 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1, -t119 * qJD(2) + t118 * t127, t119 * t128 - t122 * t127, t122 * qJD(2) - t118 * t128, t122 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1, t121 ^ 2 / 0.2e1, -t121 * t120, -t121 * t135, t120 * t135, t135 ^ 2 / 0.2e1, t115 * t120 - (t142 * t116 - t140 * t117) * t135, t115 * t121 + (t140 * t116 + t142 * t117) * t135;];
T_reg = t1;
