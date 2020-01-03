% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:58
% EndTime: 2019-12-31 19:47:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (170->40), mult. (381->86), div. (0->0), fcn. (204->6), ass. (0->34)
t145 = qJD(1) ^ 2;
t155 = t145 / 0.2e1;
t154 = -pkin(2) - qJ(4);
t144 = cos(qJ(2));
t153 = t144 * t145;
t142 = sin(qJ(2));
t146 = -qJ(3) * t142 - pkin(1);
t125 = (t154 * t144 + t146) * qJD(1);
t151 = t142 * qJD(1);
t150 = pkin(6) * t151 + qJD(3);
t126 = pkin(3) * t151 + t154 * qJD(2) + t150;
t139 = sin(pkin(8));
t140 = cos(pkin(8));
t118 = t140 * t125 + t139 * t126;
t152 = qJD(1) * t144;
t133 = -pkin(6) * t152 - qJD(2) * qJ(3);
t149 = qJD(1) * qJD(2);
t148 = t142 * t149;
t147 = t144 * t149;
t117 = -t139 * t125 + t140 * t126;
t128 = pkin(3) * t152 + qJD(4) - t133;
t143 = cos(qJ(5));
t141 = sin(qJ(5));
t134 = qJD(5) + t151;
t132 = -qJD(2) * pkin(2) + t150;
t131 = t140 * qJD(2) - t139 * t152;
t130 = t139 * qJD(2) + t140 * t152;
t129 = (-pkin(2) * t144 + t146) * qJD(1);
t121 = t130 * pkin(4) + t128;
t120 = -t141 * t130 + t143 * t131;
t119 = t143 * t130 + t141 * t131;
t116 = -t130 * pkin(7) + t118;
t115 = pkin(4) * t151 - t131 * pkin(7) + t117;
t1 = [t155, 0, 0, t142 ^ 2 * t155, t142 * t153, t148, t147, qJD(2) ^ 2 / 0.2e1, pkin(1) * t153 - pkin(6) * t148, -t145 * pkin(1) * t142 - pkin(6) * t147, (t132 * t142 - t133 * t144) * qJD(1), t132 * qJD(2) + t129 * t152, -t133 * qJD(2) - t129 * t151, t129 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t117 * t151 + t128 * t130, -t118 * t151 + t128 * t131, -t117 * t131 - t118 * t130, t118 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1, t120 ^ 2 / 0.2e1, -t120 * t119, t120 * t134, -t119 * t134, t134 ^ 2 / 0.2e1, (t143 * t115 - t141 * t116) * t134 + t121 * t119, -(t141 * t115 + t143 * t116) * t134 + t121 * t120;];
T_reg = t1;
