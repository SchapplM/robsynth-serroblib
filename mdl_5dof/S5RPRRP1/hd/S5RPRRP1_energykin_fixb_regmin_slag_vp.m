% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:57
% EndTime: 2019-12-05 17:59:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (96->28), mult. (209->63), div. (0->0), fcn. (104->4), ass. (0->24)
t94 = qJD(1) ^ 2;
t101 = t94 / 0.2e1;
t93 = cos(qJ(3));
t85 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t95 = -pkin(7) * qJD(1) + t85;
t80 = qJD(3) * pkin(3) + t95 * t93;
t91 = sin(qJ(3));
t81 = t95 * t91;
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t100 = t90 * t80 + t92 * t81;
t84 = (pkin(3) * t91 + qJ(2)) * qJD(1);
t99 = t94 * qJ(2);
t98 = qJD(3) * t85;
t97 = qJD(1) * qJD(3);
t96 = t92 * t80 - t90 * t81;
t88 = qJD(3) + qJD(4);
t87 = -pkin(1) * qJD(1) + qJD(2);
t83 = (-t90 * t91 + t92 * t93) * qJD(1);
t82 = (t90 * t93 + t91 * t92) * qJD(1);
t76 = t82 * pkin(4) + qJD(5) + t84;
t75 = -t82 * qJ(5) + t100;
t74 = t88 * pkin(4) - t83 * qJ(5) + t96;
t1 = [t101, 0, 0, t87 * qJD(1), t99, qJ(2) ^ 2 * t101 + t87 ^ 2 / 0.2e1, t93 ^ 2 * t101, -t93 * t94 * t91, t93 * t97, -t91 * t97, qJD(3) ^ 2 / 0.2e1, t91 * t99 + t93 * t98, -t91 * t98 + t93 * t99, t83 ^ 2 / 0.2e1, -t83 * t82, t83 * t88, -t82 * t88, t88 ^ 2 / 0.2e1, t84 * t82 + t96 * t88, -t100 * t88 + t84 * t83, -t74 * t83 - t75 * t82, t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1;];
T_reg = t1;
