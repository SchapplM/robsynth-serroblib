% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:12
% EndTime: 2019-12-05 18:12:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (233->42), mult. (644->88), div. (0->0), fcn. (488->8), ass. (0->37)
t152 = cos(qJ(4));
t151 = pkin(6) + qJ(2);
t136 = sin(pkin(9));
t137 = cos(pkin(9));
t140 = sin(qJ(3));
t142 = cos(qJ(3));
t126 = (t136 * t142 + t137 * t140) * qJD(1);
t148 = qJD(1) * t136;
t127 = t151 * t148;
t147 = qJD(1) * t137;
t128 = t151 * t147;
t145 = -t142 * t127 - t140 * t128;
t115 = qJD(3) * pkin(3) - t126 * pkin(7) + t145;
t125 = t140 * t148 - t142 * t147;
t149 = -t140 * t127 + t142 * t128;
t116 = -t125 * pkin(7) + t149;
t139 = sin(qJ(4));
t150 = t139 * t115 + t152 * t116;
t135 = qJD(3) + qJD(4);
t146 = t152 * t115 - t139 * t116;
t129 = qJD(2) + (-pkin(2) * t137 - pkin(1)) * qJD(1);
t120 = t125 * pkin(3) + t129;
t143 = qJD(1) ^ 2;
t141 = cos(qJ(5));
t138 = sin(qJ(5));
t134 = t137 ^ 2;
t133 = t136 ^ 2;
t132 = qJD(5) + t135;
t131 = -qJD(1) * pkin(1) + qJD(2);
t119 = -t139 * t125 + t152 * t126;
t118 = t152 * t125 + t139 * t126;
t111 = t118 * pkin(4) + t120;
t110 = -t138 * t118 + t141 * t119;
t109 = t141 * t118 + t138 * t119;
t108 = -t118 * pkin(8) + t150;
t107 = t135 * pkin(4) - t119 * pkin(8) + t146;
t1 = [t143 / 0.2e1, 0, 0, -t131 * t147, t131 * t148, (t133 + t134) * t143 * qJ(2), t131 ^ 2 / 0.2e1 + (t134 / 0.2e1 + t133 / 0.2e1) * qJ(2) ^ 2 * t143, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * qJD(3), -t125 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t145 + t129 * t125, -t149 * qJD(3) + t129 * t126, t119 ^ 2 / 0.2e1, -t119 * t118, t119 * t135, -t118 * t135, t135 ^ 2 / 0.2e1, t120 * t118 + t135 * t146, t120 * t119 - t150 * t135, t110 ^ 2 / 0.2e1, -t110 * t109, t110 * t132, -t109 * t132, t132 ^ 2 / 0.2e1, t111 * t109 + (t141 * t107 - t138 * t108) * t132, t111 * t110 - (t138 * t107 + t141 * t108) * t132;];
T_reg = t1;
