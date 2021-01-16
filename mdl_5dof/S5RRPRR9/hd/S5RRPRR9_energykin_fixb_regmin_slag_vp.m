% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:31
% EndTime: 2021-01-15 21:47:31
% DurationCPUTime: 0.25s
% Computational Cost: add. (237->41), mult. (598->90), div. (0->0), fcn. (418->8), ass. (0->39)
t78 = qJD(1) ^ 2;
t89 = t78 / 0.2e1;
t88 = cos(qJ(4));
t77 = cos(qJ(2));
t87 = t77 * t78;
t86 = pkin(6) + qJ(3);
t71 = sin(pkin(9));
t72 = cos(pkin(9));
t83 = qJD(1) * t77;
t75 = sin(qJ(2));
t84 = qJD(1) * t75;
t62 = t71 * t84 - t72 * t83;
t63 = (t71 * t77 + t72 * t75) * qJD(1);
t68 = qJD(3) + (-pkin(2) * t77 - pkin(1)) * qJD(1);
t50 = t62 * pkin(3) - t63 * pkin(7) + t68;
t66 = qJD(2) * pkin(2) - t86 * t84;
t67 = t86 * t83;
t55 = t71 * t66 + t72 * t67;
t53 = qJD(2) * pkin(7) + t55;
t74 = sin(qJ(4));
t85 = t74 * t50 + t88 * t53;
t82 = qJD(1) * qJD(2);
t81 = t75 * t82;
t80 = t77 * t82;
t79 = t88 * t50 - t74 * t53;
t54 = t72 * t66 - t71 * t67;
t61 = qJD(4) + t62;
t52 = -qJD(2) * pkin(3) - t54;
t76 = cos(qJ(5));
t73 = sin(qJ(5));
t59 = qJD(5) + t61;
t58 = t74 * qJD(2) + t88 * t63;
t57 = -t88 * qJD(2) + t74 * t63;
t47 = -t73 * t57 + t76 * t58;
t46 = t76 * t57 + t73 * t58;
t45 = t57 * pkin(4) + t52;
t44 = -t57 * pkin(8) + t85;
t43 = t61 * pkin(4) - t58 * pkin(8) + t79;
t1 = [t89, 0, 0, t75 ^ 2 * t89, t75 * t87, t81, t80, qJD(2) ^ 2 / 0.2e1, pkin(1) * t87 - pkin(6) * t81, -t78 * pkin(1) * t75 - pkin(6) * t80, t54 * qJD(2) + t68 * t62, -t55 * qJD(2) + t68 * t63, -t54 * t63 - t55 * t62, t55 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1, t58 ^ 2 / 0.2e1, -t58 * t57, t58 * t61, -t57 * t61, t61 ^ 2 / 0.2e1, t52 * t57 + t79 * t61, t52 * t58 - t85 * t61, t47 ^ 2 / 0.2e1, -t47 * t46, t47 * t59, -t46 * t59, t59 ^ 2 / 0.2e1, (t76 * t43 - t73 * t44) * t59 + t45 * t46, -(t73 * t43 + t76 * t44) * t59 + t45 * t47;];
T_reg = t1;
