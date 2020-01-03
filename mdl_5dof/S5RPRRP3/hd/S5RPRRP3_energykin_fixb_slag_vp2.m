% Calculate kinetic energy for
% S5RPRRP3
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:50
% EndTime: 2020-01-03 11:46:50
% DurationCPUTime: 0.30s
% Computational Cost: add. (196->70), mult. (436->105), div. (0->0), fcn. (238->6), ass. (0->26)
t87 = m(3) / 0.2e1;
t77 = sin(pkin(8));
t72 = (pkin(1) * t77 + pkin(6)) * qJD(1);
t82 = cos(qJ(3));
t75 = t82 * qJD(2);
t80 = sin(qJ(3));
t86 = pkin(7) * qJD(1);
t64 = qJD(3) * pkin(3) + t75 + (-t72 - t86) * t80;
t67 = t80 * qJD(2) + t82 * t72;
t65 = t82 * t86 + t67;
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t59 = t79 * t64 + t81 * t65;
t78 = cos(pkin(8));
t85 = -pkin(1) * t78 - pkin(2);
t58 = t81 * t64 - t79 * t65;
t70 = (-pkin(3) * t82 + t85) * qJD(1);
t76 = qJD(3) + qJD(4);
t73 = t85 * qJD(1);
t69 = (t79 * t82 + t80 * t81) * qJD(1);
t68 = (-t79 * t80 + t81 * t82) * qJD(1);
t66 = -t80 * t72 + t75;
t60 = -t68 * pkin(4) + qJD(5) + t70;
t57 = t68 * qJ(5) + t59;
t56 = t76 * pkin(4) - t69 * qJ(5) + t58;
t1 = qJD(2) ^ 2 * t87 + m(4) * (t66 ^ 2 + t67 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) / 0.2e1 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t58 * mrSges(5,1) + t56 * mrSges(6,1) - t59 * mrSges(5,2) - t57 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t76) * t76 + (t70 * mrSges(5,2) + t60 * mrSges(6,2) - t58 * mrSges(5,3) - t56 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t69 + (Ifges(5,5) + Ifges(6,5)) * t76) * t69 + (-t70 * mrSges(5,1) - t60 * mrSges(6,1) + t59 * mrSges(5,3) + t57 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t68 + (Ifges(5,6) + Ifges(6,6)) * t76 + (Ifges(5,4) + Ifges(6,4)) * t69) * t68 + (t73 * (-mrSges(4,1) * t82 + mrSges(4,2) * t80) + (-t66 * t80 + t67 * t82) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t80 + Ifges(4,6) * t82) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t78 * mrSges(3,1) - t77 * mrSges(3,2) + (t77 ^ 2 + t78 ^ 2) * t87 * pkin(1)) * pkin(1) + Ifges(4,2) * t82 ^ 2 / 0.2e1 + (Ifges(4,4) * t82 + Ifges(4,1) * t80 / 0.2e1) * t80) * qJD(1)) * qJD(1);
T = t1;
