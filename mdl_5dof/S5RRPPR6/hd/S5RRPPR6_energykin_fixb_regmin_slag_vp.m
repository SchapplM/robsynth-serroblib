% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:31
% EndTime: 2021-01-15 19:47:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (277->42), mult. (689->91), div. (0->0), fcn. (462->8), ass. (0->38)
t76 = qJD(1) ^ 2;
t85 = t76 / 0.2e1;
t75 = cos(qJ(2));
t84 = t75 * t76;
t83 = pkin(6) + qJ(3);
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t80 = qJD(1) * t75;
t73 = sin(qJ(2));
t81 = qJD(1) * t73;
t60 = t70 * t81 - t71 * t80;
t61 = (t70 * t75 + t71 * t73) * qJD(1);
t66 = qJD(3) + (-pkin(2) * t75 - pkin(1)) * qJD(1);
t49 = t60 * pkin(3) - t61 * qJ(4) + t66;
t64 = qJD(2) * pkin(2) - t83 * t81;
t65 = t83 * t80;
t54 = t70 * t64 + t71 * t65;
t52 = qJD(2) * qJ(4) + t54;
t69 = sin(pkin(9));
t82 = cos(pkin(9));
t43 = t69 * t49 + t82 * t52;
t79 = qJD(1) * qJD(2);
t78 = t73 * t79;
t77 = t75 * t79;
t42 = t82 * t49 - t69 * t52;
t53 = t71 * t64 - t70 * t65;
t51 = -qJD(2) * pkin(3) + qJD(4) - t53;
t74 = cos(qJ(5));
t72 = sin(qJ(5));
t59 = qJD(5) + t60;
t57 = t69 * qJD(2) + t82 * t61;
t56 = -t82 * qJD(2) + t69 * t61;
t46 = -t72 * t56 + t74 * t57;
t45 = t74 * t56 + t72 * t57;
t44 = t56 * pkin(4) + t51;
t41 = -t56 * pkin(7) + t43;
t40 = t60 * pkin(4) - t57 * pkin(7) + t42;
t1 = [t85, 0, 0, t73 ^ 2 * t85, t73 * t84, t78, t77, qJD(2) ^ 2 / 0.2e1, pkin(1) * t84 - pkin(6) * t78, -t76 * pkin(1) * t73 - pkin(6) * t77, t53 * qJD(2) + t66 * t60, -t54 * qJD(2) + t66 * t61, -t53 * t61 - t54 * t60, t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1, t42 * t60 + t51 * t56, -t43 * t60 + t51 * t57, -t42 * t57 - t43 * t56, t43 ^ 2 / 0.2e1 + t42 ^ 2 / 0.2e1 + t51 ^ 2 / 0.2e1, t46 ^ 2 / 0.2e1, -t46 * t45, t46 * t59, -t45 * t59, t59 ^ 2 / 0.2e1, (t74 * t40 - t72 * t41) * t59 + t44 * t45, -(t72 * t40 + t74 * t41) * t59 + t44 * t46;];
T_reg = t1;
