% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:10
% EndTime: 2019-12-31 18:22:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (153->35), mult. (361->82), div. (0->0), fcn. (211->8), ass. (0->32)
t136 = qJD(1) ^ 2;
t145 = t136 / 0.2e1;
t144 = cos(pkin(9));
t130 = sin(pkin(8));
t123 = (pkin(1) * t130 + pkin(6)) * qJD(1);
t133 = sin(qJ(3));
t135 = cos(qJ(3));
t143 = t133 * qJD(2) + t135 * t123;
t116 = qJD(3) * qJ(4) + t143;
t131 = cos(pkin(8));
t139 = -pkin(1) * t131 - pkin(2);
t117 = (-pkin(3) * t135 - qJ(4) * t133 + t139) * qJD(1);
t129 = sin(pkin(9));
t108 = t144 * t116 + t129 * t117;
t142 = qJD(1) * t133;
t141 = t135 * qJD(1);
t140 = qJD(1) * qJD(3);
t107 = -t129 * t116 + t144 * t117;
t138 = t135 * qJD(2) - t133 * t123;
t115 = -qJD(3) * pkin(3) + qJD(4) - t138;
t134 = cos(qJ(5));
t132 = sin(qJ(5));
t125 = -qJD(5) + t141;
t124 = t139 * qJD(1);
t120 = t129 * qJD(3) + t144 * t142;
t119 = -t144 * qJD(3) + t129 * t142;
t111 = -t132 * t119 + t134 * t120;
t110 = t134 * t119 + t132 * t120;
t109 = t119 * pkin(4) + t115;
t106 = -t119 * pkin(7) + t108;
t105 = -pkin(4) * t141 - t120 * pkin(7) + t107;
t1 = [t145, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t130 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t136, t133 ^ 2 * t145, t133 * t136 * t135, t133 * t140, t135 * t140, qJD(3) ^ 2 / 0.2e1, t138 * qJD(3) - t124 * t141, -t143 * qJD(3) + t124 * t142, -t107 * t141 + t115 * t119, t108 * t141 + t115 * t120, -t107 * t120 - t108 * t119, t108 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1, t111 ^ 2 / 0.2e1, -t111 * t110, -t111 * t125, t110 * t125, t125 ^ 2 / 0.2e1, -(t134 * t105 - t132 * t106) * t125 + t109 * t110, (t132 * t105 + t134 * t106) * t125 + t109 * t111;];
T_reg = t1;
