% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:39
% EndTime: 2019-12-05 15:28:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (100->30), mult. (247->63), div. (0->0), fcn. (149->4), ass. (0->22)
t91 = cos(pkin(8));
t89 = t91 * qJD(1);
t90 = sin(pkin(8));
t97 = qJD(2) * t90;
t78 = t89 + (-pkin(6) - qJ(3)) * t97;
t95 = qJ(3) * qJD(2);
t83 = t90 * qJD(1) + t91 * t95;
t96 = qJD(2) * t91;
t79 = pkin(6) * t96 + t83;
t92 = sin(qJ(4));
t93 = cos(qJ(4));
t98 = t92 * t78 + t93 * t79;
t94 = t93 * t78 - t92 * t79;
t84 = qJD(3) + (-pkin(3) * t91 - pkin(2)) * qJD(2);
t87 = -qJD(2) * pkin(2) + qJD(3);
t82 = -t90 * t95 + t89;
t81 = (t90 * t93 + t91 * t92) * qJD(2);
t80 = t92 * t97 - t93 * t96;
t75 = qJD(4) * qJ(5) + t98;
t74 = t80 * pkin(4) - t81 * qJ(5) + t84;
t73 = -qJD(4) * pkin(4) + qJD(5) - t94;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, -t87 * t96, t87 * t97, (-t82 * t90 + t83 * t91) * qJD(2), t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1, t81 ^ 2 / 0.2e1, -t81 * t80, t81 * qJD(4), -t80 * qJD(4), qJD(4) ^ 2 / 0.2e1, t94 * qJD(4) + t84 * t80, -t98 * qJD(4) + t84 * t81, -t73 * qJD(4) + t74 * t80, t73 * t81 - t75 * t80, t75 * qJD(4) - t74 * t81, t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1;];
T_reg = t1;
