% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:34
% EndTime: 2019-12-05 16:17:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (81->21), mult. (125->53), div. (0->0), fcn. (61->6), ass. (0->22)
t85 = qJD(2) + qJD(3);
t83 = t85 ^ 2;
t86 = sin(pkin(9));
t100 = t83 * t86 ^ 2;
t99 = t85 * t86;
t87 = cos(pkin(9));
t98 = t87 * t85;
t97 = pkin(2) * qJD(2);
t88 = sin(qJ(5));
t96 = t88 * t99;
t90 = cos(qJ(5));
t95 = t90 * t99;
t94 = sin(qJ(3)) * t97;
t93 = cos(qJ(3)) * t97;
t92 = qJD(4) - t93;
t80 = -qJD(5) + t98;
t79 = t85 * qJ(4) + t94;
t78 = -t85 * pkin(3) + t92;
t77 = t86 * qJD(1) + t87 * t79;
t75 = -t87 * qJD(1) + t86 * t79;
t74 = (-pkin(4) * t87 - pkin(7) * t86 - pkin(3)) * t85 + t92;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t83 / 0.2e1, t85 * t93, -t85 * t94, -t78 * t98, t78 * t99, (t75 * t86 + t77 * t87) * t85, t77 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1, t90 ^ 2 * t100 / 0.2e1, -t90 * t88 * t100, -t80 * t95, t80 * t96, t80 ^ 2 / 0.2e1, -(t90 * t74 - t88 * t77) * t80 + t75 * t96, (t88 * t74 + t90 * t77) * t80 + t75 * t95;];
T_reg = t1;
