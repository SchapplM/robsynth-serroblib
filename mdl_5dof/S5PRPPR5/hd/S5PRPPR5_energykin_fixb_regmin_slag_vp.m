% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (62->21), mult. (113->50), div. (0->0), fcn. (49->6), ass. (0->19)
t90 = qJD(2) ^ 2;
t95 = t90 / 0.2e1;
t89 = cos(qJ(2));
t91 = -t89 * qJD(1) + qJD(3);
t80 = (-pkin(2) - pkin(3)) * qJD(2) + t91;
t87 = sin(qJ(2));
t83 = qJD(2) * qJ(3) + t87 * qJD(1);
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t78 = t84 * t80 + t85 * t83;
t77 = t85 * t80 - t84 * t83;
t94 = (qJD(2) * pkin(4) - t77) * qJD(2);
t93 = qJD(1) * qJD(2);
t92 = qJD(2) * qJD(5);
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t82 = -qJD(2) * pkin(2) + t91;
t76 = -qJD(2) * pkin(6) + t78;
t1 = [qJD(1) ^ 2 / 0.2e1, t95, t89 * t93, -t87 * t93, -t82 * qJD(2), t83 * qJD(2), t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1, -t77 * qJD(2), t78 * qJD(2), t78 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t86 ^ 2 * t95, t86 * t90 * t88, -t86 * t92, -t88 * t92, qJD(5) ^ 2 / 0.2e1, (t88 * qJD(4) - t86 * t76) * qJD(5) + t88 * t94, -(t86 * qJD(4) + t88 * t76) * qJD(5) - t86 * t94;];
T_reg = t1;
