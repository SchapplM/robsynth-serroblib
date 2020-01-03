% Calculate kinetic energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:34
% EndTime: 2019-12-31 19:52:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (148->52), mult. (197->77), div. (0->0), fcn. (56->4), ass. (0->17)
t60 = qJD(1) + qJD(2);
t62 = sin(qJ(2));
t69 = pkin(1) * qJD(1);
t68 = t62 * t69;
t57 = t60 * qJ(3) + t68;
t71 = t57 ^ 2;
t64 = cos(qJ(2));
t67 = -t64 * t69 + qJD(3);
t55 = (-pkin(2) - pkin(7)) * t60 + t67;
t70 = t55 * mrSges(5,3);
t63 = cos(qJ(4));
t61 = sin(qJ(4));
t56 = -t60 * pkin(2) + t67;
t53 = qJD(4) * qJ(5) + t61 * t55;
t52 = -qJD(4) * pkin(4) - t63 * t55 + qJD(5);
t51 = t68 + (pkin(4) * t61 - qJ(5) * t63 + qJ(3)) * t60;
t1 = m(5) * (t71 + (t61 ^ 2 + t63 ^ 2) * t55 ^ 2) / 0.2e1 + m(6) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t71) / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t60 ^ 2 + (Ifges(2,3) / 0.2e1 + m(3) * (t62 ^ 2 + t64 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (-t52 * mrSges(6,1) + t53 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4) + (t63 * mrSges(5,1) - t61 * mrSges(5,2)) * t55) * qJD(4) + (t56 * mrSges(4,2) + t57 * mrSges(4,3) + (mrSges(3,1) * t64 - mrSges(3,2) * t62) * t69 + (t57 * mrSges(5,2) + t52 * mrSges(6,2) - t51 * mrSges(6,3) + (-t70 + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t60) * t63 + (Ifges(6,4) + Ifges(5,5)) * qJD(4)) * t63 + (t57 * mrSges(5,1) + t51 * mrSges(6,1) - t53 * mrSges(6,2) + (-t70 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t60) * t61 + (-Ifges(5,4) + Ifges(6,5)) * t60 * t63 + (-Ifges(5,6) + Ifges(6,6)) * qJD(4)) * t61) * t60;
T = t1;
