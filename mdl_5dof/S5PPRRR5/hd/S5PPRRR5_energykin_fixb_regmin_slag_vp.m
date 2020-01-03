% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:44
% EndTime: 2019-12-31 17:35:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->13), mult. (73->39), div. (0->0), fcn. (37->6), ass. (0->19)
t68 = qJD(3) + qJD(4);
t67 = t68 ^ 2;
t82 = t67 / 0.2e1;
t75 = cos(qJ(3));
t65 = qJD(3) * pkin(3) + t75 * qJD(2);
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t72 = sin(qJ(3));
t79 = qJD(2) * t72;
t76 = t74 * t65 - t71 * t79;
t81 = (-t68 * pkin(4) - t76) * t68;
t80 = t71 * t65 + t74 * t79;
t78 = qJD(5) * t68;
t77 = qJD(2) * qJD(3);
t73 = cos(qJ(5));
t70 = sin(qJ(5));
t69 = qJD(1) ^ 2 / 0.2e1;
t63 = t68 * pkin(7) + t80;
t1 = [t69, t69 + qJD(2) ^ 2 / 0.2e1, qJD(3) ^ 2 / 0.2e1, t75 * t77, -t72 * t77, t82, t76 * t68, -t80 * t68, t70 ^ 2 * t82, t70 * t67 * t73, t70 * t78, t73 * t78, qJD(5) ^ 2 / 0.2e1, -t73 * t81 + (-t73 * qJD(1) - t70 * t63) * qJD(5), t70 * t81 - (-t70 * qJD(1) + t73 * t63) * qJD(5);];
T_reg = t1;
