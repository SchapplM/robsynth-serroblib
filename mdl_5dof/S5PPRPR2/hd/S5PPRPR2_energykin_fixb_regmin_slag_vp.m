% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->18), mult. (97->44), div. (0->0), fcn. (53->6), ass. (0->19)
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t75 = sin(qJ(3));
t77 = cos(qJ(3));
t86 = (t72 * t77 + t73 * t75) * qJD(1);
t78 = qJD(3) ^ 2;
t85 = t78 / 0.2e1;
t68 = qJD(3) * qJ(4) + t86;
t84 = t68 * qJD(3);
t83 = qJD(3) * qJD(5);
t82 = (-t72 * t75 + t73 * t77) * qJD(1);
t80 = qJD(4) - t82;
t79 = qJD(1) ^ 2;
t76 = cos(qJ(5));
t74 = sin(qJ(5));
t71 = qJD(2) ^ 2 / 0.2e1;
t67 = -qJD(3) * pkin(3) + t80;
t66 = (-pkin(3) - pkin(6)) * qJD(3) + t80;
t1 = [t79 / 0.2e1, t71 + (t72 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1) * t79, t85, t82 * qJD(3), -qJD(3) * t86, t67 * qJD(3), t84, t71 + t68 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1, t76 ^ 2 * t85, -t76 * t78 * t74, t76 * t83, -t74 * t83, qJD(5) ^ 2 / 0.2e1, (-t74 * qJD(2) + t76 * t66) * qJD(5) + t74 * t84, -(t76 * qJD(2) + t74 * t66) * qJD(5) + t76 * t84;];
T_reg = t1;
