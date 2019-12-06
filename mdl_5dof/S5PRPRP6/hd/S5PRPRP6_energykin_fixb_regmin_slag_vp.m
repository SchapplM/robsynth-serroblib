% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:18
% EndTime: 2019-12-05 15:41:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (62->22), mult. (117->53), div. (0->0), fcn. (44->4), ass. (0->20)
t90 = qJD(2) ^ 2;
t98 = t90 / 0.2e1;
t86 = sin(qJ(4));
t88 = cos(qJ(4));
t87 = sin(qJ(2));
t94 = t87 * qJD(1);
t80 = t94 + (pkin(4) * t86 - qJ(5) * t88 + qJ(3)) * qJD(2);
t97 = qJD(2) * t80;
t89 = cos(qJ(2));
t91 = -t89 * qJD(1) + qJD(3);
t82 = (-pkin(2) - pkin(6)) * qJD(2) + t91;
t96 = qJD(4) * t82;
t84 = qJD(2) * qJ(3) + t94;
t95 = t84 * qJD(2);
t93 = qJD(1) * qJD(2);
t92 = qJD(2) * qJD(4);
t83 = -qJD(2) * pkin(2) + t91;
t81 = qJD(4) * qJ(5) + t86 * t82;
t79 = -qJD(4) * pkin(4) - t88 * t82 + qJD(5);
t1 = [qJD(1) ^ 2 / 0.2e1, t98, t89 * t93, -t87 * t93, t83 * qJD(2), t95, t84 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1, t88 ^ 2 * t98, -t88 * t90 * t86, t88 * t92, -t86 * t92, qJD(4) ^ 2 / 0.2e1, t86 * t95 + t88 * t96, -t86 * t96 + t88 * t95, -t79 * qJD(4) + t86 * t97, (t79 * t88 - t81 * t86) * qJD(2), t81 * qJD(4) - t88 * t97, t81 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1;];
T_reg = t1;
