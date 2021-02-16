% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:54
% EndTime: 2021-01-15 23:23:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (306->42), mult. (695->91), div. (0->0), fcn. (476->8), ass. (0->38)
t62 = qJD(1) ^ 2;
t73 = t62 / 0.2e1;
t72 = cos(qJ(3));
t61 = cos(qJ(2));
t71 = t61 * t62;
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t68 = qJD(1) * t59;
t46 = t58 * qJD(2) + t72 * t68;
t67 = t61 * qJD(1);
t52 = -qJD(3) + t67;
t44 = (-pkin(2) * t61 - pkin(7) * t59 - pkin(1)) * qJD(1);
t49 = pkin(6) * t67 + qJD(2) * pkin(7);
t63 = t72 * t44 - t58 * t49;
t34 = -t52 * pkin(3) - t46 * qJ(4) + t63;
t45 = -t72 * qJD(2) + t58 * t68;
t70 = t58 * t44 + t72 * t49;
t36 = -t45 * qJ(4) + t70;
t56 = sin(pkin(9));
t69 = cos(pkin(9));
t28 = t56 * t34 + t69 * t36;
t66 = qJD(1) * qJD(2);
t65 = t59 * t66;
t64 = t61 * t66;
t27 = t69 * t34 - t56 * t36;
t48 = -qJD(2) * pkin(2) + pkin(6) * t68;
t40 = t45 * pkin(3) + qJD(4) + t48;
t60 = cos(qJ(5));
t57 = sin(qJ(5));
t50 = -qJD(5) + t52;
t39 = -t56 * t45 + t69 * t46;
t38 = t69 * t45 + t56 * t46;
t31 = t38 * pkin(4) + t40;
t30 = -t57 * t38 + t60 * t39;
t29 = t60 * t38 + t57 * t39;
t26 = -t38 * pkin(8) + t28;
t25 = -t52 * pkin(4) - t39 * pkin(8) + t27;
t1 = [t73, 0, 0, t59 ^ 2 * t73, t59 * t71, t65, t64, qJD(2) ^ 2 / 0.2e1, pkin(1) * t71 - pkin(6) * t65, -t62 * pkin(1) * t59 - pkin(6) * t64, t46 ^ 2 / 0.2e1, -t46 * t45, -t46 * t52, t45 * t52, t52 ^ 2 / 0.2e1, t48 * t45 - t63 * t52, t48 * t46 + t70 * t52, -t27 * t52 + t40 * t38, t28 * t52 + t40 * t39, -t27 * t39 - t28 * t38, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1, t30 ^ 2 / 0.2e1, -t30 * t29, -t30 * t50, t29 * t50, t50 ^ 2 / 0.2e1, -(t60 * t25 - t57 * t26) * t50 + t31 * t29, (t57 * t25 + t60 * t26) * t50 + t31 * t30;];
T_reg = t1;
