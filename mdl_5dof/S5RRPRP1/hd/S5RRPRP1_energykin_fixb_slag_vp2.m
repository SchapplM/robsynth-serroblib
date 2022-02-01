% Calculate kinetic energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:29
% EndTime: 2022-01-20 10:19:29
% DurationCPUTime: 0.18s
% Computational Cost: add. (173->56), mult. (274->82), div. (0->0), fcn. (114->6), ass. (0->22)
t72 = qJD(1) + qJD(2);
t78 = cos(qJ(2));
t83 = pkin(1) * qJD(1);
t66 = pkin(2) * t72 + t78 * t83;
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t76 = sin(qJ(2));
t81 = t76 * t83;
t64 = t73 * t66 + t74 * t81;
t62 = pkin(7) * t72 + t64;
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t59 = t75 * qJD(3) + t77 * t62;
t82 = qJ(5) * t72;
t63 = t74 * t66 - t73 * t81;
t70 = t77 * qJD(3);
t61 = -pkin(3) * t72 - t63;
t58 = -t62 * t75 + t70;
t57 = qJD(5) + (-pkin(4) * t77 - pkin(3)) * t72 - t63;
t56 = t77 * t82 + t59;
t55 = qJD(4) * pkin(4) + t70 + (-t62 - t82) * t75;
t1 = m(5) * (t58 ^ 2 + t59 ^ 2 + t61 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t76 ^ 2 + t78 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t58 * mrSges(5,1) + t55 * mrSges(6,1) - t59 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4)) * qJD(4) + (t63 * mrSges(4,1) - t64 * mrSges(4,2) + (mrSges(3,1) * t78 - mrSges(3,2) * t76) * t83 + (-t61 * mrSges(5,1) - t57 * mrSges(6,1) + t59 * mrSges(5,3) + t56 * mrSges(6,3) + (Ifges(5,6) + Ifges(6,6)) * qJD(4)) * t77 + (t61 * mrSges(5,2) + t57 * mrSges(6,2) - t58 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,5) + Ifges(6,5)) * qJD(4)) * t75 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t77 ^ 2 + ((Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t75 + (Ifges(5,4) + Ifges(6,4)) * t77) * t75) * t72) * t72;
T = t1;
