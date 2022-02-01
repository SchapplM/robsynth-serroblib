% Calculate kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:09
% EndTime: 2022-01-23 09:12:09
% DurationCPUTime: 0.27s
% Computational Cost: add. (153->63), mult. (355->93), div. (0->0), fcn. (178->6), ass. (0->25)
t72 = sin(pkin(7));
t68 = (pkin(1) * t72 + qJ(3)) * qJD(1);
t73 = cos(pkin(8));
t70 = t73 * qJD(2);
t71 = sin(pkin(8));
t64 = t71 * t68 - t70;
t84 = t64 ^ 2;
t83 = m(3) / 0.2e1;
t74 = cos(pkin(7));
t80 = -pkin(1) * t74 - pkin(2);
t62 = qJD(3) + (-pkin(3) * t73 - pkin(6) * t71 + t80) * qJD(1);
t66 = t71 * qJD(2) + t73 * t68;
t75 = sin(qJ(4));
t76 = cos(qJ(4));
t58 = t75 * t62 + t76 * t66;
t82 = t71 * qJD(1);
t81 = t73 * qJD(1);
t79 = qJ(5) * t82;
t57 = t76 * t62 - t75 * t66;
t69 = qJD(4) - t81;
t67 = t80 * qJD(1) + qJD(3);
t59 = qJD(5) - t70 + (pkin(4) * qJD(1) * t75 + t68) * t71;
t56 = -t75 * t79 + t58;
t55 = t69 * pkin(4) - t76 * t79 + t57;
t1 = qJD(2) ^ 2 * t83 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t84) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t84) / 0.2e1 + (t57 * mrSges(5,1) + t55 * mrSges(6,1) - t58 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t69) * t69 + ((-t67 * mrSges(4,1) + t66 * mrSges(4,3) + Ifges(4,2) * t81 / 0.2e1) * t73 + (t67 * mrSges(4,2) + t64 * mrSges(4,3) + (t64 * mrSges(5,2) + t59 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t76 * t82 + (Ifges(5,5) + Ifges(6,5)) * t69) * t76 + (t64 * mrSges(5,1) + t59 * mrSges(6,1) - t58 * mrSges(5,3) - t56 * mrSges(6,3) + (-Ifges(5,6) - Ifges(6,6)) * t69 + ((Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t75 + (-Ifges(5,4) - Ifges(6,4)) * t76) * t82) * t75) * t71 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t74 * mrSges(3,1) - t72 * mrSges(3,2) + (t72 ^ 2 + t74 ^ 2) * t83 * pkin(1)) * pkin(1) + (Ifges(4,4) * t73 + Ifges(4,1) * t71 / 0.2e1) * t71) * qJD(1)) * qJD(1);
T = t1;
