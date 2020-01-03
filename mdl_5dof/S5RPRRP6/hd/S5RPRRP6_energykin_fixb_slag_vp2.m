% Calculate kinetic energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:06
% EndTime: 2019-12-31 18:42:06
% DurationCPUTime: 0.32s
% Computational Cost: add. (192->69), mult. (408->104), div. (0->0), fcn. (210->6), ass. (0->25)
t86 = m(3) / 0.2e1;
t76 = sin(pkin(8));
t72 = (pkin(1) * t76 + pkin(6)) * qJD(1);
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t67 = t79 * qJD(2) + t81 * t72;
t64 = qJD(3) * pkin(7) + t67;
t77 = cos(pkin(8));
t84 = -pkin(1) * t77 - pkin(2);
t65 = (-pkin(3) * t81 - pkin(7) * t79 + t84) * qJD(1);
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t58 = t80 * t64 + t78 * t65;
t85 = qJD(1) * t79;
t57 = -t78 * t64 + t80 * t65;
t66 = t81 * qJD(2) - t79 * t72;
t63 = -qJD(3) * pkin(3) - t66;
t74 = -t81 * qJD(1) + qJD(4);
t73 = t84 * qJD(1);
t71 = t78 * qJD(3) + t80 * t85;
t70 = t80 * qJD(3) - t78 * t85;
t59 = -t70 * pkin(4) + qJD(5) + t63;
t56 = t70 * qJ(5) + t58;
t55 = t74 * pkin(4) - t71 * qJ(5) + t57;
t1 = qJD(2) ^ 2 * t86 + m(4) * (t66 ^ 2 + t67 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t63 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t57 * mrSges(5,1) + t55 * mrSges(6,1) - t58 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t74) * t74 + (t63 * mrSges(5,2) + t59 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t71 + (Ifges(5,5) + Ifges(6,5)) * t74) * t71 + (-t63 * mrSges(5,1) - t59 * mrSges(6,1) + t58 * mrSges(5,3) + t56 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t70 + (Ifges(5,6) + Ifges(6,6)) * t74 + (Ifges(5,4) + Ifges(6,4)) * t71) * t70 + (t73 * (-mrSges(4,1) * t81 + mrSges(4,2) * t79) + (-t66 * t79 + t67 * t81) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t79 + Ifges(4,6) * t81) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t77 * mrSges(3,1) - t76 * mrSges(3,2) + (t76 ^ 2 + t77 ^ 2) * t86 * pkin(1)) * pkin(1) + Ifges(4,2) * t81 ^ 2 / 0.2e1 + (Ifges(4,4) * t81 + Ifges(4,1) * t79 / 0.2e1) * t79) * qJD(1)) * qJD(1);
T = t1;
