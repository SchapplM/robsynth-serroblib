% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:14
% EndTime: 2019-12-05 17:38:14
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->24), mult. (127->59), div. (0->0), fcn. (51->4), ass. (0->22)
t78 = qJD(1) ^ 2;
t84 = t78 / 0.2e1;
t83 = -pkin(1) - qJ(3);
t70 = qJD(1) * qJ(2) + qJD(3);
t67 = -qJD(1) * pkin(6) + t70;
t82 = qJD(4) * t67;
t68 = -t83 * qJD(1) - qJD(2);
t81 = t68 * qJD(1);
t80 = qJD(1) * qJD(4);
t79 = -pkin(7) * qJD(1) + t67;
t77 = cos(qJ(4));
t76 = cos(qJ(5));
t75 = sin(qJ(4));
t74 = sin(qJ(5));
t72 = qJD(4) + qJD(5);
t71 = -qJD(1) * pkin(1) + qJD(2);
t66 = -qJD(2) + (pkin(4) * t75 - t83) * qJD(1);
t65 = (-t74 * t75 + t76 * t77) * qJD(1);
t64 = (t74 * t77 + t75 * t76) * qJD(1);
t63 = t79 * t75;
t62 = qJD(4) * pkin(4) + t79 * t77;
t1 = [t84, 0, 0, t71 * qJD(1), t78 * qJ(2), qJ(2) ^ 2 * t84 + t71 ^ 2 / 0.2e1, t70 * qJD(1), t81, t68 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1, t77 ^ 2 * t84, -t77 * t78 * t75, t77 * t80, -t75 * t80, qJD(4) ^ 2 / 0.2e1, t75 * t81 + t77 * t82, -t75 * t82 + t77 * t81, t65 ^ 2 / 0.2e1, -t65 * t64, t65 * t72, -t64 * t72, t72 ^ 2 / 0.2e1, t66 * t64 + (t76 * t62 - t74 * t63) * t72, t66 * t65 - (t74 * t62 + t76 * t63) * t72;];
T_reg = t1;
