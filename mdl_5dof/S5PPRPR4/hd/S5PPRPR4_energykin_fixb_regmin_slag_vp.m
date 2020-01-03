% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:23
% EndTime: 2019-12-31 17:32:23
% DurationCPUTime: 0.04s
% Computational Cost: add. (49->24), mult. (123->54), div. (0->0), fcn. (73->6), ass. (0->22)
t80 = sin(pkin(8));
t90 = qJD(3) * t80;
t81 = cos(pkin(8));
t89 = qJD(3) * t81;
t88 = t81 * qJD(1);
t87 = qJD(2) * qJD(3);
t83 = sin(qJ(3));
t77 = qJD(3) * qJ(4) + t83 * qJD(2);
t71 = -t80 * qJD(1) + t81 * t77;
t85 = cos(qJ(3));
t86 = -t85 * qJD(2) + qJD(4);
t84 = cos(qJ(5));
t82 = sin(qJ(5));
t79 = qJD(1) ^ 2 / 0.2e1;
t76 = -qJD(3) * pkin(3) + t86;
t74 = (-pkin(4) * t81 - pkin(3)) * qJD(3) + t86;
t73 = (t80 * t84 + t81 * t82) * qJD(3);
t72 = t82 * t90 - t84 * t89;
t70 = -t80 * t77 - t88;
t69 = pkin(6) * t89 + t71;
t68 = -t88 + (-pkin(6) * qJD(3) - t77) * t80;
t1 = [t79, t79 + qJD(2) ^ 2 / 0.2e1, qJD(3) ^ 2 / 0.2e1, t85 * t87, -t83 * t87, -t76 * t89, t76 * t90, (-t70 * t80 + t71 * t81) * qJD(3), t71 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1, t73 ^ 2 / 0.2e1, -t73 * t72, t73 * qJD(5), -t72 * qJD(5), qJD(5) ^ 2 / 0.2e1, t74 * t72 + (t84 * t68 - t82 * t69) * qJD(5), t74 * t73 - (t82 * t68 + t84 * t69) * qJD(5);];
T_reg = t1;
