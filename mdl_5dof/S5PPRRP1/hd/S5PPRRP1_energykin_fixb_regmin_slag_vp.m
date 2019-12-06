% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP1
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
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:13
% EndTime: 2019-12-05 15:07:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (42->21), mult. (122->54), div. (0->0), fcn. (72->6), ass. (0->21)
t92 = qJD(3) ^ 2;
t102 = t92 / 0.2e1;
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t98 = qJD(1) * cos(qJ(3));
t99 = qJD(1) * sin(qJ(3));
t100 = t86 * t98 + t87 * t99;
t80 = qJD(3) * pkin(6) + t100;
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t101 = t88 * qJD(2) + t90 * t80;
t94 = -t86 * t99 + t87 * t98;
t97 = qJD(3) * (-qJD(3) * pkin(3) - t94);
t96 = qJ(5) * qJD(3);
t95 = qJD(3) * qJD(4);
t93 = qJD(1) ^ 2;
t85 = t90 * qJD(2);
t77 = qJD(5) + (-pkin(4) * t90 - pkin(3)) * qJD(3) - t94;
t76 = t90 * t96 + t101;
t75 = qJD(4) * pkin(4) + t85 + (-t80 - t96) * t88;
t1 = [t93 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t86 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1) * t93, t102, t94 * qJD(3), -t100 * qJD(3), t88 ^ 2 * t102, t88 * t92 * t90, t88 * t95, t90 * t95, qJD(4) ^ 2 / 0.2e1, (-t88 * t80 + t85) * qJD(4) - t90 * t97, -t101 * qJD(4) + t88 * t97, (-t75 * t88 + t76 * t90) * qJD(3), t76 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1;];
T_reg = t1;
