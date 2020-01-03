% Calculate kinetic energy for
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:52
% EndTime: 2019-12-31 19:50:52
% DurationCPUTime: 0.21s
% Computational Cost: add. (265->60), mult. (390->94), div. (0->0), fcn. (206->6), ass. (0->25)
t76 = cos(pkin(8));
t87 = t76 ^ 2;
t86 = cos(qJ(4));
t75 = sin(pkin(8));
t74 = qJD(1) + qJD(2);
t78 = sin(qJ(2));
t84 = pkin(1) * qJD(1);
t70 = t74 * qJ(3) + t78 * t84;
t83 = pkin(7) * t74 + t70;
t63 = t83 * t75;
t64 = t83 * t76;
t77 = sin(qJ(4));
t60 = -t77 * t63 + t86 * t64;
t85 = t75 ^ 2 + t87;
t79 = cos(qJ(2));
t82 = -t79 * t84 + qJD(3);
t59 = -t63 * t86 - t77 * t64;
t65 = (-pkin(3) * t76 - pkin(2)) * t74 + t82;
t68 = -t74 * pkin(2) + t82;
t67 = (t75 * t86 + t76 * t77) * t74;
t66 = (t75 * t77 - t76 * t86) * t74;
t58 = qJD(4) * qJ(5) + t60;
t57 = -qJD(4) * pkin(4) + qJD(5) - t59;
t56 = t66 * pkin(4) - t67 * qJ(5) + t65;
t1 = m(4) * (t70 ^ 2 * t85 + t68 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t65 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t78 ^ 2 + t79 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t65 * mrSges(5,2) + t57 * mrSges(6,2) - t59 * mrSges(5,3) - t56 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t67) * t67 + (t65 * mrSges(5,1) + t56 * mrSges(6,1) - t58 * mrSges(6,2) - t60 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t66 + (-Ifges(5,4) + Ifges(6,5)) * t67) * t66 + (t59 * mrSges(5,1) - t57 * mrSges(6,1) - t60 * mrSges(5,2) + t58 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(4) + (Ifges(6,4) + Ifges(5,5)) * t67 + (-Ifges(5,6) + Ifges(6,6)) * t66) * qJD(4) + (t68 * (-mrSges(4,1) * t76 + mrSges(4,2) * t75) + (Ifges(4,2) * t87 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t76 + Ifges(4,1) * t75 / 0.2e1) * t75) * t74 + (mrSges(3,1) * t79 - mrSges(3,2) * t78) * t84 + t85 * t70 * mrSges(4,3)) * t74;
T = t1;
