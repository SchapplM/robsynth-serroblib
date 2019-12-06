% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:06
% EndTime: 2019-12-05 14:58:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (36->16), mult. (106->49), div. (0->0), fcn. (70->8), ass. (0->21)
t90 = qJD(4) ^ 2;
t97 = t90 / 0.2e1;
t82 = sin(pkin(9));
t84 = cos(pkin(9));
t83 = sin(pkin(8));
t95 = qJD(1) * t83;
t79 = t84 * qJD(2) - t82 * t95;
t80 = t82 * qJD(2) + t84 * t95;
t87 = sin(qJ(4));
t89 = cos(qJ(4));
t96 = t87 * t79 + t89 * t80;
t92 = t89 * t79 - t87 * t80;
t94 = (-qJD(4) * pkin(4) - t92) * qJD(4);
t93 = qJD(4) * qJD(5);
t91 = qJD(1) ^ 2;
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t85 = cos(pkin(8));
t81 = -t85 * qJD(1) + qJD(3);
t76 = qJD(4) * pkin(6) + t96;
t1 = [t91 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t83 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1) * t91, t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1, t97, t92 * qJD(4), -t96 * qJD(4), t86 ^ 2 * t97, t86 * t90 * t88, t86 * t93, t88 * t93, qJD(5) ^ 2 / 0.2e1, (-t86 * t76 + t88 * t81) * qJD(5) - t88 * t94, -(t88 * t76 + t86 * t81) * qJD(5) + t86 * t94;];
T_reg = t1;
