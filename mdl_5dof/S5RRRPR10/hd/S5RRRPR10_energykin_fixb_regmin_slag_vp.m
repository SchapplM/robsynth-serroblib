% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:15
% EndTime: 2021-01-15 23:41:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (435->45), mult. (1073->97), div. (0->0), fcn. (819->10), ass. (0->42)
t89 = cos(qJ(3));
t69 = sin(pkin(5));
t77 = qJD(1) ^ 2;
t88 = t69 ^ 2 * t77;
t84 = cos(pkin(5)) * qJD(1);
t66 = qJD(2) + t84;
t73 = sin(qJ(3));
t74 = sin(qJ(2));
t85 = qJD(1) * t69;
t81 = t74 * t85;
t58 = t73 * t66 + t89 * t81;
t76 = cos(qJ(2));
t80 = t76 * t85;
t61 = -qJD(3) + t80;
t82 = pkin(1) * t84;
t86 = pkin(7) * t80 + t74 * t82;
t54 = t66 * pkin(8) + t86;
t56 = (-pkin(2) * t76 - pkin(8) * t74 - pkin(1)) * t85;
t79 = -t73 * t54 + t89 * t56;
t41 = -t61 * pkin(3) - t58 * qJ(4) + t79;
t57 = -t89 * t66 + t73 * t81;
t87 = t89 * t54 + t73 * t56;
t43 = -t57 * qJ(4) + t87;
t68 = sin(pkin(10));
t70 = cos(pkin(10));
t38 = t68 * t41 + t70 * t43;
t83 = t76 * t88;
t48 = t70 * t57 + t68 * t58;
t37 = t70 * t41 - t68 * t43;
t78 = -pkin(7) * t81 + t76 * t82;
t53 = -t66 * pkin(2) - t78;
t47 = t57 * pkin(3) + qJD(4) + t53;
t75 = cos(qJ(5));
t72 = sin(qJ(5));
t49 = -t68 * t57 + t70 * t58;
t46 = qJD(5) + t48;
t45 = t75 * t49 - t72 * t61;
t44 = t72 * t49 + t75 * t61;
t39 = t48 * pkin(4) - t49 * pkin(9) + t47;
t36 = -t61 * pkin(9) + t38;
t35 = t61 * pkin(4) - t37;
t1 = [t77 / 0.2e1, 0, 0, t74 ^ 2 * t88 / 0.2e1, t74 * t83, t66 * t81, t66 * t80, t66 ^ 2 / 0.2e1, pkin(1) * t83 + t78 * t66, -pkin(1) * t74 * t88 - t86 * t66, t58 ^ 2 / 0.2e1, -t58 * t57, -t58 * t61, t57 * t61, t61 ^ 2 / 0.2e1, t53 * t57 - t79 * t61, t53 * t58 + t87 * t61, -t37 * t61 + t47 * t48, t38 * t61 + t47 * t49, -t37 * t49 - t38 * t48, t38 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1, t45 ^ 2 / 0.2e1, -t45 * t44, t45 * t46, -t44 * t46, t46 ^ 2 / 0.2e1, (-t72 * t36 + t75 * t39) * t46 + t35 * t44, -(t75 * t36 + t72 * t39) * t46 + t35 * t45;];
T_reg = t1;
