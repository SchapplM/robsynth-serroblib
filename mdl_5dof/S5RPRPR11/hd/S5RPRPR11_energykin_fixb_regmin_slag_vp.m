% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:54
% EndTime: 2019-12-31 18:27:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (154->40), mult. (416->81), div. (0->0), fcn. (267->6), ass. (0->33)
t129 = sin(pkin(8));
t130 = cos(pkin(8));
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t118 = (t129 * t134 + t130 * t132) * qJD(1);
t121 = -(pkin(2) * t130 + pkin(1)) * qJD(1) + qJD(2);
t146 = -t118 * qJ(4) + t121;
t145 = -pkin(3) - pkin(4);
t144 = pkin(6) + qJ(2);
t141 = qJD(1) * t129;
t119 = t144 * t141;
t140 = qJD(1) * t130;
t120 = t144 * t140;
t142 = -t132 * t119 + t134 * t120;
t112 = qJD(3) * qJ(4) + t142;
t138 = -t134 * t119 - t132 * t120;
t137 = qJD(4) - t138;
t135 = qJD(1) ^ 2;
t133 = cos(qJ(5));
t131 = sin(qJ(5));
t126 = qJD(3) - qJD(5);
t125 = t130 ^ 2;
t124 = t129 ^ 2;
t123 = -qJD(1) * pkin(1) + qJD(2);
t117 = t132 * t141 - t134 * t140;
t111 = -qJD(3) * pkin(3) + t137;
t110 = t131 * t117 + t133 * t118;
t109 = -t133 * t117 + t131 * t118;
t108 = t117 * pkin(3) + t146;
t107 = t117 * pkin(7) + t112;
t106 = -t118 * pkin(7) + t145 * qJD(3) + t137;
t105 = t145 * t117 - t146;
t1 = [t135 / 0.2e1, 0, 0, -t123 * t140, t123 * t141, (t124 + t125) * t135 * qJ(2), t123 ^ 2 / 0.2e1 + (t125 / 0.2e1 + t124 / 0.2e1) * qJ(2) ^ 2 * t135, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * qJD(3), -t117 * qJD(3), qJD(3) ^ 2 / 0.2e1, t138 * qJD(3) + t121 * t117, -t142 * qJD(3) + t121 * t118, -t111 * qJD(3) + t108 * t117, t111 * t118 - t112 * t117, t112 * qJD(3) - t108 * t118, t112 ^ 2 / 0.2e1 + t108 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1, t110 ^ 2 / 0.2e1, -t110 * t109, -t110 * t126, t109 * t126, t126 ^ 2 / 0.2e1, t105 * t109 - (t133 * t106 - t131 * t107) * t126, t105 * t110 + (t131 * t106 + t133 * t107) * t126;];
T_reg = t1;
