% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:27
% EndTime: 2018-11-14 13:57:27
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->12), mult. (52->26), div. (0->0), fcn. (26->4), ass. (0->11)
t34 = sin(pkin(5));
t35 = cos(pkin(5));
t36 = sin(qJ(3));
t37 = cos(qJ(3));
t41 = (t34 * t36 - t35 * t37) * qJD(1);
t40 = (t34 * t37 + t35 * t36) * qJD(1);
t38 = qJD(1) ^ 2;
t33 = qJD(2) ^ 2 / 0.2e1;
t30 = qJD(3) * qJ(4) + t40;
t29 = -qJD(3) * pkin(3) + qJD(4) + t41;
t1 = [t38 / 0.2e1, t33 + (t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * t38, qJD(3) ^ 2 / 0.2e1, -qJD(3) * t41, -t40 * qJD(3), -t29 * qJD(3), t30 * qJD(3), t30 ^ 2 / 0.2e1 + t33 + t29 ^ 2 / 0.2e1;];
T_reg  = t1;
