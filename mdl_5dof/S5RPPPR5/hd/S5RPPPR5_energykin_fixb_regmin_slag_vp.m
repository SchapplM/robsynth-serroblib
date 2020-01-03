% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:29
% EndTime: 2019-12-31 17:46:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (105->29), mult. (204->67), div. (0->0), fcn. (92->6), ass. (0->25)
t98 = qJD(1) ^ 2;
t101 = t98 / 0.2e1;
t85 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t92 = sin(pkin(7));
t94 = cos(pkin(7));
t99 = qJ(2) * qJD(1);
t81 = t92 * t85 + t94 * t99;
t79 = -qJD(1) * qJ(4) + t81;
t91 = sin(pkin(8));
t93 = cos(pkin(8));
t75 = t91 * qJD(3) + t93 * t79;
t100 = qJD(1) * t93;
t80 = t94 * t85 - t92 * t99;
t78 = qJD(1) * pkin(3) + qJD(4) - t80;
t97 = cos(qJ(5));
t96 = sin(qJ(5));
t90 = t93 * qJD(3);
t88 = -qJD(1) * pkin(1) + qJD(2);
t83 = (-t91 * t97 - t93 * t96) * qJD(1);
t82 = (t91 * t96 - t93 * t97) * qJD(1);
t76 = pkin(4) * t100 + t78;
t74 = -t91 * t79 + t90;
t73 = -pkin(6) * t100 + t75;
t72 = t90 + (pkin(6) * qJD(1) - t79) * t91;
t1 = [t101, 0, 0, -t88 * qJD(1), t98 * qJ(2), qJ(2) ^ 2 * t101 + t88 ^ 2 / 0.2e1, -t80 * qJD(1), t81 * qJD(1), t81 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t78 * t100, -t78 * t91 * qJD(1), (t74 * t91 - t75 * t93) * qJD(1), t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1, t83 ^ 2 / 0.2e1, t83 * t82, t83 * qJD(5), t82 * qJD(5), qJD(5) ^ 2 / 0.2e1, -t76 * t82 + (t97 * t72 - t96 * t73) * qJD(5), t76 * t83 - (t96 * t72 + t97 * t73) * qJD(5);];
T_reg = t1;
