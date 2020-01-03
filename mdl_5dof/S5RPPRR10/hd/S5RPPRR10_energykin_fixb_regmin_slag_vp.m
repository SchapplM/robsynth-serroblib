% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (116->37), mult. (314->74), div. (0->0), fcn. (193->6), ass. (0->33)
t140 = qJD(1) ^ 2;
t148 = t140 / 0.2e1;
t147 = qJ(2) * t140;
t133 = sin(pkin(8));
t145 = qJD(1) * t133;
t122 = qJ(2) * t145 + qJD(3);
t120 = -pkin(6) * t145 + t122;
t134 = cos(pkin(8));
t144 = qJD(1) * t134;
t121 = (-pkin(6) + qJ(2)) * t144;
t137 = sin(qJ(4));
t139 = cos(qJ(4));
t146 = t137 * t120 + t139 * t121;
t129 = -qJD(1) * pkin(1) + qJD(2);
t143 = qJ(2) ^ 2 * t148;
t142 = t139 * t120 - t137 * t121;
t116 = -pkin(2) * t144 - qJ(3) * t145 + t129;
t112 = pkin(3) * t144 - t116;
t138 = cos(qJ(5));
t136 = sin(qJ(5));
t132 = qJD(4) + qJD(5);
t131 = t134 ^ 2;
t130 = t133 ^ 2;
t126 = t131 * t147;
t123 = t131 * t143;
t119 = (t133 * t139 - t134 * t137) * qJD(1);
t118 = (t133 * t137 + t134 * t139) * qJD(1);
t111 = -t136 * t118 + t138 * t119;
t110 = t138 * t118 + t136 * t119;
t109 = t118 * pkin(4) + t112;
t108 = -t118 * pkin(7) + t146;
t107 = qJD(4) * pkin(4) - t119 * pkin(7) + t142;
t1 = [t148, 0, 0, -t129 * t144, t129 * t145, t130 * t147 + t126, t123 + t130 * t143 + t129 ^ 2 / 0.2e1, -t116 * t144, t122 * t145 + t126, -t116 * t145, t123 + t116 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1, t119 ^ 2 / 0.2e1, -t119 * t118, t119 * qJD(4), -t118 * qJD(4), qJD(4) ^ 2 / 0.2e1, t142 * qJD(4) + t112 * t118, -t146 * qJD(4) + t112 * t119, t111 ^ 2 / 0.2e1, -t111 * t110, t111 * t132, -t110 * t132, t132 ^ 2 / 0.2e1, t109 * t110 + (t138 * t107 - t136 * t108) * t132, t109 * t111 - (t136 * t107 + t138 * t108) * t132;];
T_reg = t1;
