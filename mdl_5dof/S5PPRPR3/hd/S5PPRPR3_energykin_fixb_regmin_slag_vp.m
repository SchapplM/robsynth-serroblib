% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x13]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:21
% EndTime: 2019-12-05 15:05:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (42->17), mult. (112->50), div. (0->0), fcn. (70->8), ass. (0->22)
t102 = qJD(3) ^ 2;
t108 = t102 / 0.2e1;
t101 = cos(qJ(3));
t95 = sin(pkin(8));
t107 = qJD(1) * t95;
t99 = sin(qJ(3));
t104 = t101 * qJD(2) - t99 * t107;
t90 = qJD(3) * pkin(3) + t104;
t91 = t99 * qJD(2) + t101 * t107;
t94 = sin(pkin(9));
t96 = cos(pkin(9));
t87 = t94 * t90 + t96 * t91;
t86 = t96 * t90 - t94 * t91;
t106 = (-qJD(3) * pkin(4) - t86) * qJD(3);
t105 = qJD(3) * qJD(5);
t103 = qJD(1) ^ 2;
t100 = cos(qJ(5));
t98 = sin(qJ(5));
t97 = cos(pkin(8));
t92 = -t97 * qJD(1) + qJD(4);
t85 = qJD(3) * pkin(6) + t87;
t1 = [t103 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t95 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1) * t103, t108, t104 * qJD(3), -t91 * qJD(3), t87 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1, t98 ^ 2 * t108, t98 * t102 * t100, t98 * t105, t100 * t105, qJD(5) ^ 2 / 0.2e1, (t100 * t92 - t98 * t85) * qJD(5) - t100 * t106, -(t100 * t85 + t98 * t92) * qJD(5) + t98 * t106;];
T_reg = t1;
