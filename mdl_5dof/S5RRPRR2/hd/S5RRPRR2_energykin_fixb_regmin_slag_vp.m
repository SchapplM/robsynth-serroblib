% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:26
% EndTime: 2019-12-05 18:28:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (244->39), mult. (650->86), div. (0->0), fcn. (482->8), ass. (0->37)
t159 = qJD(1) * (pkin(6) + qJ(3));
t149 = qJD(1) ^ 2;
t158 = t149 / 0.2e1;
t157 = cos(qJ(4));
t148 = cos(qJ(2));
t155 = t148 * t149;
t146 = sin(qJ(2));
t137 = qJD(2) * pkin(2) - t146 * t159;
t138 = t148 * t159;
t142 = sin(pkin(9));
t143 = cos(pkin(9));
t128 = t143 * t137 - t142 * t138;
t135 = (t142 * t148 + t143 * t146) * qJD(1);
t123 = qJD(2) * pkin(3) - t135 * pkin(7) + t128;
t129 = t142 * t137 + t143 * t138;
t134 = (-t142 * t146 + t143 * t148) * qJD(1);
t124 = t134 * pkin(7) + t129;
t145 = sin(qJ(4));
t154 = t145 * t123 + t157 * t124;
t153 = qJD(1) * qJD(2);
t141 = qJD(2) + qJD(4);
t152 = t146 * t153;
t151 = t148 * t153;
t150 = t157 * t123 - t145 * t124;
t139 = qJD(3) + (-pkin(2) * t148 - pkin(1)) * qJD(1);
t130 = -t134 * pkin(3) + t139;
t147 = cos(qJ(5));
t144 = sin(qJ(5));
t140 = qJD(5) + t141;
t127 = t145 * t134 + t157 * t135;
t126 = -t157 * t134 + t145 * t135;
t119 = t126 * pkin(4) + t130;
t118 = -t144 * t126 + t147 * t127;
t117 = t147 * t126 + t144 * t127;
t116 = -t126 * pkin(8) + t154;
t115 = t141 * pkin(4) - t127 * pkin(8) + t150;
t1 = [t158, 0, 0, t146 ^ 2 * t158, t146 * t155, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t155 - pkin(6) * t152, -t149 * pkin(1) * t146 - pkin(6) * t151, -t128 * t135 + t129 * t134, t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t141, -t126 * t141, t141 ^ 2 / 0.2e1, t130 * t126 + t150 * t141, t130 * t127 - t154 * t141, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * t140, -t117 * t140, t140 ^ 2 / 0.2e1, t119 * t117 + (t147 * t115 - t144 * t116) * t140, t119 * t118 - (t144 * t115 + t147 * t116) * t140;];
T_reg = t1;
