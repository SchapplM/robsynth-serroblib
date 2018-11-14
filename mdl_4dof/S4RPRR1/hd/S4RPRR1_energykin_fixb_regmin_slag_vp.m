% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:34
% EndTime: 2018-11-14 13:50:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (31->12), mult. (72->31), div. (0->0), fcn. (30->6), ass. (0->15)
t52 = qJD(1) + qJD(3);
t53 = sin(pkin(7));
t62 = pkin(1) * qJD(1) * t53;
t54 = cos(pkin(7));
t50 = (pkin(1) * t54 + pkin(2)) * qJD(1);
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t61 = t58 * t50 - t56 * t62;
t59 = qJD(1) ^ 2;
t57 = cos(qJ(4));
t55 = sin(qJ(4));
t51 = qJD(4) + t52;
t48 = t56 * t50 + t58 * t62;
t47 = t52 * pkin(3) + t61;
t1 = [t59 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t53 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t59, t52 ^ 2 / 0.2e1, t61 * t52, -t48 * t52, t51 ^ 2 / 0.2e1 (t57 * t47 - t55 * t48) * t51 -(t55 * t47 + t57 * t48) * t51;];
T_reg  = t1;
