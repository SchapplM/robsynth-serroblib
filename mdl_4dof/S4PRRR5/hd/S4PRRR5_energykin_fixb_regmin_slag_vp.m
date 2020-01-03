% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR5
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
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:45
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (33->10), mult. (67->35), div. (0->0), fcn. (35->6), ass. (0->19)
t63 = qJD(2) + qJD(3);
t62 = t63 ^ 2;
t77 = t62 / 0.2e1;
t69 = cos(qJ(2));
t60 = qJD(2) * pkin(2) + t69 * qJD(1);
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t66 = sin(qJ(2));
t74 = qJD(1) * t66;
t70 = t68 * t60 - t65 * t74;
t76 = (-t63 * pkin(3) - t70) * t63;
t75 = t65 * t60 + t68 * t74;
t64 = sin(qJ(4));
t73 = qJD(4) * t64;
t67 = cos(qJ(4));
t72 = qJD(4) * t67;
t71 = qJD(1) * qJD(2);
t58 = t63 * pkin(6) + t75;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t69 * t71, -t66 * t71, t77, t70 * t63, -t75 * t63, t64 ^ 2 * t77, t64 * t62 * t67, t63 * t73, t63 * t72, qJD(4) ^ 2 / 0.2e1, -t58 * t73 - t67 * t76, -t58 * t72 + t64 * t76;];
T_reg = t1;
