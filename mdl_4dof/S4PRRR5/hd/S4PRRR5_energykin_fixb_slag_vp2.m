% Calculate kinetic energy for
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:33
% EndTime: 2019-12-31 16:33:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (75->31), mult. (134->57), div. (0->0), fcn. (60->6), ass. (0->16)
t57 = qJD(2) + qJD(3);
t67 = t57 / 0.2e1;
t63 = cos(qJ(2));
t55 = qJD(2) * pkin(2) + t63 * qJD(1);
t59 = sin(qJ(3));
t62 = cos(qJ(3));
t60 = sin(qJ(2));
t65 = qJD(1) * t60;
t53 = t59 * t55 + t62 * t65;
t51 = t57 * pkin(6) + t53;
t66 = t51 * mrSges(5,3);
t52 = t62 * t55 - t59 * t65;
t61 = cos(qJ(4));
t58 = sin(qJ(4));
t50 = -t57 * pkin(3) - t52;
t1 = m(4) * (t52 ^ 2 + t53 ^ 2) / 0.2e1 + m(5) * (t50 ^ 2 + (t58 ^ 2 + t61 ^ 2) * t51 ^ 2) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t60 ^ 2 + t63 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(5,3) * qJD(4) / 0.2e1 + (-t58 * mrSges(5,1) - t61 * mrSges(5,2)) * t51) * qJD(4) + (Ifges(3,3) * qJD(2) / 0.2e1 + (t63 * mrSges(3,1) - t60 * mrSges(3,2)) * qJD(1)) * qJD(2) + (-t53 * mrSges(4,2) + t52 * mrSges(4,1) + Ifges(4,3) * t67 + (-t50 * mrSges(5,1) + Ifges(5,6) * qJD(4) + (Ifges(5,2) * t67 + t66) * t61) * t61 + (Ifges(5,4) * t61 * t57 + t50 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t67 + t66) * t58) * t58) * t57;
T = t1;
