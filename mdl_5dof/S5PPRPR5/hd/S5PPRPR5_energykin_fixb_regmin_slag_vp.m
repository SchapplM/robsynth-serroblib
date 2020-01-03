% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:29
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->15), mult. (61->36), div. (0->0), fcn. (21->4), ass. (0->15)
t61 = qJD(3) ^ 2;
t66 = t61 / 0.2e1;
t58 = sin(qJ(3));
t54 = qJD(3) * qJ(4) + t58 * qJD(2);
t65 = t54 * qJD(3);
t64 = qJD(2) * qJD(3);
t63 = qJD(3) * qJD(5);
t60 = cos(qJ(3));
t62 = -t60 * qJD(2) + qJD(4);
t59 = cos(qJ(5));
t57 = sin(qJ(5));
t56 = qJD(1) ^ 2 / 0.2e1;
t53 = -qJD(3) * pkin(3) + t62;
t52 = (-pkin(3) - pkin(6)) * qJD(3) + t62;
t1 = [t56, t56 + qJD(2) ^ 2 / 0.2e1, t66, t60 * t64, -t58 * t64, t53 * qJD(3), t65, t56 + t54 ^ 2 / 0.2e1 + t53 ^ 2 / 0.2e1, t59 ^ 2 * t66, -t59 * t61 * t57, t59 * t63, -t57 * t63, qJD(5) ^ 2 / 0.2e1, t57 * t65 + (t57 * qJD(1) + t59 * t52) * qJD(5), t59 * t65 - (-t59 * qJD(1) + t57 * t52) * qJD(5);];
T_reg = t1;
