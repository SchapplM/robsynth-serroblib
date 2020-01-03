% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:06
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->18), mult. (110->53), div. (0->0), fcn. (63->6), ass. (0->24)
t82 = qJD(2) ^ 2;
t90 = t82 / 0.2e1;
t77 = sin(qJ(3));
t89 = qJD(2) * t77;
t80 = cos(qJ(3));
t88 = qJD(2) * t80;
t78 = sin(qJ(2));
t72 = qJD(2) * pkin(5) + t78 * qJD(1);
t87 = qJD(3) * t72;
t81 = cos(qJ(2));
t86 = t81 * qJD(1);
t85 = qJD(1) * qJD(2);
t84 = qJD(2) * qJD(3);
t83 = pkin(6) * qJD(2) + t72;
t79 = cos(qJ(4));
t76 = sin(qJ(4));
t75 = qJD(3) + qJD(4);
t73 = -qJD(2) * pkin(2) - t86;
t71 = -t86 + (-pkin(3) * t80 - pkin(2)) * qJD(2);
t70 = (t76 * t80 + t77 * t79) * qJD(2);
t69 = t76 * t89 - t79 * t88;
t68 = t83 * t80;
t67 = qJD(3) * pkin(3) - t83 * t77;
t1 = [qJD(1) ^ 2 / 0.2e1, t90, t81 * t85, -t78 * t85, t77 ^ 2 * t90, t77 * t82 * t80, t77 * t84, t80 * t84, qJD(3) ^ 2 / 0.2e1, -t73 * t88 - t77 * t87, t73 * t89 - t80 * t87, t70 ^ 2 / 0.2e1, -t70 * t69, t70 * t75, -t69 * t75, t75 ^ 2 / 0.2e1, (t79 * t67 - t76 * t68) * t75 + t71 * t69, -(t76 * t67 + t79 * t68) * t75 + t71 * t70;];
T_reg = t1;
