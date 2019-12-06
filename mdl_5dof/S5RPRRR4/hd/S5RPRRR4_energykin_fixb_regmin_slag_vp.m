% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:34
% EndTime: 2019-12-05 18:14:34
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->18), mult. (167->51), div. (0->0), fcn. (81->8), ass. (0->24)
t77 = qJD(1) + qJD(3);
t76 = qJD(4) + t77;
t75 = t76 ^ 2;
t94 = t75 / 0.2e1;
t79 = cos(pkin(9));
t74 = (pkin(1) * t79 + pkin(2)) * qJD(1);
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t78 = sin(pkin(9));
t90 = pkin(1) * qJD(1) * t78;
t88 = t85 * t74 - t82 * t90;
t70 = t77 * pkin(3) + t88;
t72 = t82 * t74 + t85 * t90;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t89 = t84 * t70 - t81 * t72;
t93 = (-t76 * pkin(4) - t89) * t76;
t92 = t81 * t70 + t84 * t72;
t91 = qJD(5) * t76;
t86 = qJD(1) ^ 2;
t83 = cos(qJ(5));
t80 = sin(qJ(5));
t68 = t76 * pkin(8) + t92;
t1 = [t86 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t78 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t86, t77 ^ 2 / 0.2e1, t88 * t77, -t72 * t77, t94, t89 * t76, -t92 * t76, t80 ^ 2 * t94, t80 * t75 * t83, t80 * t91, t83 * t91, qJD(5) ^ 2 / 0.2e1, -t83 * t93 + (qJD(2) * t83 - t68 * t80) * qJD(5), t80 * t93 - (qJD(2) * t80 + t68 * t83) * qJD(5);];
T_reg = t1;
