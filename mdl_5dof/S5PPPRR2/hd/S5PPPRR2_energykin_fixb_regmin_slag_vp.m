% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% T_reg [1x13]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->16), mult. (98->49), div. (0->0), fcn. (62->8), ass. (0->21)
t101 = qJD(4) ^ 2;
t108 = t101 / 0.2e1;
t100 = cos(qJ(4));
t94 = sin(pkin(8));
t106 = qJD(1) * t94;
t93 = sin(pkin(9));
t95 = cos(pkin(9));
t90 = t93 * qJD(2) + t95 * t106;
t96 = cos(pkin(8));
t92 = -t96 * qJD(1) + qJD(3);
t98 = sin(qJ(4));
t107 = t100 * t90 + t98 * t92;
t103 = t100 * t92 - t98 * t90;
t105 = (-qJD(4) * pkin(4) - t103) * qJD(4);
t104 = qJD(4) * qJD(5);
t102 = qJD(1) ^ 2;
t99 = cos(qJ(5));
t97 = sin(qJ(5));
t88 = -t95 * qJD(2) + t93 * t106;
t86 = qJD(4) * pkin(6) + t107;
t1 = [t102 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t94 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1) * t102, t90 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1, t108, t103 * qJD(4), -t107 * qJD(4), t97 ^ 2 * t108, t97 * t101 * t99, t97 * t104, t99 * t104, qJD(5) ^ 2 / 0.2e1, (-t97 * t86 + t99 * t88) * qJD(5) - t99 * t105, -(t99 * t86 + t97 * t88) * qJD(5) + t97 * t105;];
T_reg = t1;
