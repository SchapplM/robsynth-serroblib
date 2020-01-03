% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->20), mult. (121->49), div. (0->0), fcn. (41->6), ass. (0->21)
t75 = -qJD(1) + qJD(4);
t74 = t75 ^ 2;
t90 = t74 / 0.2e1;
t78 = cos(pkin(8));
t86 = -pkin(1) * t78 - pkin(2);
t70 = qJD(3) + (-pkin(3) + t86) * qJD(1);
t77 = sin(pkin(8));
t73 = (pkin(1) * t77 + qJ(3)) * qJD(1);
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t85 = t82 * t70 - t80 * t73;
t89 = (-t75 * pkin(4) - t85) * t75;
t88 = t80 * t70 + t82 * t73;
t87 = qJD(5) * t75;
t83 = qJD(1) ^ 2;
t81 = cos(qJ(5));
t79 = sin(qJ(5));
t76 = qJD(2) ^ 2 / 0.2e1;
t72 = t86 * qJD(1) + qJD(3);
t68 = t75 * pkin(7) + t88;
t1 = [t83 / 0.2e1, 0, 0, t76 + (t77 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t83, -t72 * qJD(1), t73 * qJD(1), t73 ^ 2 / 0.2e1 + t76 + t72 ^ 2 / 0.2e1, t90, t85 * t75, -t88 * t75, t79 ^ 2 * t90, t79 * t74 * t81, t79 * t87, t81 * t87, qJD(5) ^ 2 / 0.2e1, -t81 * t89 + (-t81 * qJD(2) - t79 * t68) * qJD(5), t79 * t89 - (-t79 * qJD(2) + t81 * t68) * qJD(5);];
T_reg = t1;
