% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:09
% EndTime: 2019-12-05 15:09:09
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->22), mult. (158->58), div. (0->0), fcn. (95->6), ass. (0->21)
t99 = qJD(3) ^ 2;
t110 = t99 / 0.2e1;
t106 = qJD(1) * cos(qJ(3));
t107 = qJD(1) * sin(qJ(3));
t93 = sin(pkin(8));
t94 = cos(pkin(8));
t108 = t93 * t106 + t94 * t107;
t88 = qJD(3) * pkin(6) + t108;
t95 = sin(qJ(4));
t97 = cos(qJ(4));
t109 = t95 * qJD(2) + t97 * t88;
t102 = t94 * t106 - t93 * t107;
t85 = (-pkin(4) * t97 - qJ(5) * t95 - pkin(3)) * qJD(3) - t102;
t105 = qJD(3) * t85;
t104 = (-qJD(3) * pkin(3) - t102) * qJD(3);
t103 = qJD(3) * qJD(4);
t101 = t97 * qJD(2) - t95 * t88;
t100 = qJD(1) ^ 2;
t84 = qJD(4) * qJ(5) + t109;
t83 = -qJD(4) * pkin(4) + qJD(5) - t101;
t1 = [t100 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t93 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1) * t100, t110, t102 * qJD(3), -t108 * qJD(3), t95 ^ 2 * t110, t95 * t99 * t97, t95 * t103, t97 * t103, qJD(4) ^ 2 / 0.2e1, t101 * qJD(4) - t97 * t104, -t109 * qJD(4) + t95 * t104, -t83 * qJD(4) - t97 * t105, (t83 * t95 + t84 * t97) * qJD(3), t84 * qJD(4) - t95 * t105, t84 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1;];
T_reg = t1;
