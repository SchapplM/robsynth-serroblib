% Calculate kinetic energy for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:02
% EndTime: 2019-12-31 17:49:02
% DurationCPUTime: 0.27s
% Computational Cost: add. (171->65), mult. (409->97), div. (0->0), fcn. (228->6), ass. (0->25)
t85 = m(3) / 0.2e1;
t84 = cos(qJ(4));
t76 = sin(pkin(7));
t71 = (pkin(1) * t76 + qJ(3)) * qJD(1);
t77 = cos(pkin(8));
t74 = t77 * qJD(2);
t75 = sin(pkin(8));
t62 = t74 + (-pkin(6) * qJD(1) - t71) * t75;
t65 = t75 * qJD(2) + t77 * t71;
t83 = qJD(1) * t77;
t63 = pkin(6) * t83 + t65;
t79 = sin(qJ(4));
t58 = t79 * t62 + t84 * t63;
t78 = cos(pkin(7));
t82 = -pkin(1) * t78 - pkin(2);
t57 = t84 * t62 - t79 * t63;
t66 = qJD(3) + (-pkin(3) * t77 + t82) * qJD(1);
t70 = t82 * qJD(1) + qJD(3);
t68 = (t84 * t75 + t77 * t79) * qJD(1);
t67 = t79 * t75 * qJD(1) - t84 * t83;
t64 = -t75 * t71 + t74;
t59 = t67 * pkin(4) - t68 * qJ(5) + t66;
t56 = qJD(4) * qJ(5) + t58;
t55 = -qJD(4) * pkin(4) + qJD(5) - t57;
t1 = m(4) * (t64 ^ 2 + t65 ^ 2 + t70 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t85 + m(5) * (t57 ^ 2 + t58 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (t66 * mrSges(5,2) + t55 * mrSges(6,2) - t57 * mrSges(5,3) - t59 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t68) * t68 + (t66 * mrSges(5,1) + t59 * mrSges(6,1) - t56 * mrSges(6,2) - t58 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t67 + (-Ifges(5,4) + Ifges(6,5)) * t68) * t67 + (t57 * mrSges(5,1) - t55 * mrSges(6,1) - t58 * mrSges(5,2) + t56 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4) + (Ifges(6,4) + Ifges(5,5)) * t68 + (-Ifges(5,6) + Ifges(6,6)) * t67) * qJD(4) + (t70 * (-mrSges(4,1) * t77 + mrSges(4,2) * t75) + (-t64 * t75 + t65 * t77) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t78 * mrSges(3,1) - t76 * mrSges(3,2) + (t76 ^ 2 + t78 ^ 2) * t85 * pkin(1)) * pkin(1) + Ifges(4,2) * t77 ^ 2 / 0.2e1 + (Ifges(4,4) * t77 + Ifges(4,1) * t75 / 0.2e1) * t75) * qJD(1)) * qJD(1);
T = t1;
