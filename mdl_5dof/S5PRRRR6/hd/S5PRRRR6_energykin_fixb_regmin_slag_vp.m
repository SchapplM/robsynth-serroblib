% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (101->23), mult. (164->62), div. (0->0), fcn. (101->8), ass. (0->30)
t105 = qJD(2) + qJD(3);
t103 = t105 ^ 2;
t123 = t103 / 0.2e1;
t113 = cos(qJ(2));
t100 = qJD(2) * pkin(2) + t113 * qJD(1);
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t109 = sin(qJ(2));
t119 = qJD(1) * t109;
t122 = t108 * t100 + t112 * t119;
t107 = sin(qJ(4));
t121 = t105 * t107;
t111 = cos(qJ(4));
t120 = t105 * t111;
t118 = qJD(4) * t107;
t117 = qJD(4) * t111;
t116 = qJD(1) * qJD(2);
t95 = t105 * pkin(7) + t122;
t115 = pkin(8) * t105 + t95;
t114 = t112 * t100 - t108 * t119;
t110 = cos(qJ(5));
t106 = sin(qJ(5));
t104 = qJD(4) + qJD(5);
t97 = (t106 * t111 + t107 * t110) * t105;
t96 = t106 * t121 - t110 * t120;
t94 = -t105 * pkin(3) - t114;
t93 = (-pkin(4) * t111 - pkin(3)) * t105 - t114;
t92 = t115 * t111;
t91 = qJD(4) * pkin(4) - t115 * t107;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t113 * t116, -t109 * t116, t123, t114 * t105, -t122 * t105, t107 ^ 2 * t123, t107 * t103 * t111, t105 * t118, t105 * t117, qJD(4) ^ 2 / 0.2e1, -t95 * t118 - t94 * t120, -t95 * t117 + t94 * t121, t97 ^ 2 / 0.2e1, -t97 * t96, t97 * t104, -t96 * t104, t104 ^ 2 / 0.2e1, (-t106 * t92 + t110 * t91) * t104 + t93 * t96, -(t106 * t91 + t110 * t92) * t104 + t93 * t97;];
T_reg = t1;
