% Calculate kinetic energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:46
% EndTime: 2019-12-31 16:20:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (73->39), mult. (184->70), div. (0->0), fcn. (100->4), ass. (0->17)
t62 = sin(pkin(7));
t63 = cos(pkin(7));
t67 = qJD(2) * qJ(3);
t56 = t62 * qJD(1) + t63 * t67;
t65 = cos(qJ(4));
t64 = sin(qJ(4));
t61 = t63 * qJD(1);
t59 = -qJD(2) * pkin(2) + qJD(3);
t57 = qJD(3) + (-pkin(3) * t63 - pkin(2)) * qJD(2);
t55 = -t62 * t67 + t61;
t54 = (t62 * t65 + t63 * t64) * qJD(2);
t53 = (-t62 * t64 + t63 * t65) * qJD(2);
t52 = t63 * qJD(2) * pkin(5) + t56;
t51 = t61 + (-pkin(5) - qJ(3)) * t62 * qJD(2);
t50 = t64 * t51 + t65 * t52;
t49 = t65 * t51 - t64 * t52;
t1 = m(5) * (t49 ^ 2 + t50 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t57 * mrSges(5,2) - t49 * mrSges(5,3) + Ifges(5,1) * t54 / 0.2e1) * t54 + (-t57 * mrSges(5,1) + t50 * mrSges(5,3) + Ifges(5,4) * t54 + Ifges(5,2) * t53 / 0.2e1) * t53 + (t49 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,5) * t54 + Ifges(5,6) * t53 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t59 * (-mrSges(4,1) * t63 + mrSges(4,2) * t62) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t63 ^ 2 / 0.2e1 + (Ifges(4,4) * t63 + Ifges(4,1) * t62 / 0.2e1) * t62) * qJD(2) + (-t55 * t62 + t56 * t63) * mrSges(4,3)) * qJD(2);
T = t1;
