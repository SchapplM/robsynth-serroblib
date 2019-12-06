% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:48
% EndTime: 2019-12-05 15:14:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->26), mult. (162->66), div. (0->0), fcn. (107->8), ass. (0->27)
t111 = qJD(3) ^ 2;
t121 = t111 / 0.2e1;
t103 = sin(pkin(9));
t104 = cos(pkin(9));
t117 = qJD(1) * cos(qJ(3));
t118 = qJD(1) * sin(qJ(3));
t120 = t103 * t117 + t104 * t118;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t93 = qJD(3) * pkin(6) + t120;
t119 = t106 * qJD(2) + t109 * t93;
t116 = qJD(3) * t106;
t115 = qJD(3) * t109;
t114 = qJD(3) * qJD(4);
t113 = -t103 * t118 + t104 * t117;
t112 = qJD(1) ^ 2;
t108 = cos(qJ(5));
t105 = sin(qJ(5));
t102 = qJD(4) + qJD(5);
t101 = t109 * qJD(2);
t95 = (t105 * t109 + t106 * t108) * qJD(3);
t94 = t105 * t116 - t108 * t115;
t92 = -qJD(3) * pkin(3) - t113;
t90 = (-pkin(4) * t109 - pkin(3)) * qJD(3) - t113;
t89 = pkin(7) * t115 + t119;
t88 = qJD(4) * pkin(4) + t101 + (-pkin(7) * qJD(3) - t93) * t106;
t1 = [t112 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t103 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1) * t112, t121, t113 * qJD(3), -t120 * qJD(3), t106 ^ 2 * t121, t106 * t111 * t109, t106 * t114, t109 * t114, qJD(4) ^ 2 / 0.2e1, (-t106 * t93 + t101) * qJD(4) - t92 * t115, -t119 * qJD(4) + t92 * t116, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * t102, -t94 * t102, t102 ^ 2 / 0.2e1, (-t105 * t89 + t108 * t88) * t102 + t90 * t94, -(t105 * t88 + t108 * t89) * t102 + t90 * t95;];
T_reg = t1;
