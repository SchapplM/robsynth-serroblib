% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:54
% EndTime: 2019-12-05 15:26:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->19), mult. (105->45), div. (0->0), fcn. (53->6), ass. (0->21)
t86 = qJD(2) ^ 2;
t92 = t86 / 0.2e1;
t83 = sin(qJ(2));
t91 = qJD(1) * t83;
t85 = cos(qJ(2));
t77 = qJD(2) * pkin(2) + t85 * qJD(1);
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t76 = t80 * t77 + t81 * t91;
t73 = qJD(2) * qJ(4) + t76;
t90 = t73 * qJD(2);
t89 = qJD(1) * qJD(2);
t88 = qJD(2) * qJD(5);
t75 = t81 * t77 - t80 * t91;
t87 = qJD(4) - t75;
t84 = cos(qJ(5));
t82 = sin(qJ(5));
t79 = qJD(3) ^ 2 / 0.2e1;
t72 = -qJD(2) * pkin(3) + t87;
t71 = (-pkin(3) - pkin(6)) * qJD(2) + t87;
t1 = [qJD(1) ^ 2 / 0.2e1, t92, t85 * t89, -t83 * t89, t76 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1 + t79, t72 * qJD(2), t90, t79 + t73 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1, t84 ^ 2 * t92, -t84 * t86 * t82, t84 * t88, -t82 * t88, qJD(5) ^ 2 / 0.2e1, (-t82 * qJD(3) + t84 * t71) * qJD(5) + t82 * t90, -(t84 * qJD(3) + t82 * t71) * qJD(5) + t84 * t90;];
T_reg = t1;
