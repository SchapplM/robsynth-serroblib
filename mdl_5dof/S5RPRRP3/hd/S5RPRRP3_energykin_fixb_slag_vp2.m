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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:06
% EndTime: 2019-12-05 18:03:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (196->70), mult. (436->105), div. (0->0), fcn. (238->6), ass. (0->26)
t85 = m(3) / 0.2e1;
t75 = sin(pkin(8));
t70 = (pkin(1) * t75 + pkin(6)) * qJD(1);
t80 = cos(qJ(3));
t73 = t80 * qJD(2);
t78 = sin(qJ(3));
t84 = pkin(7) * qJD(1);
t62 = qJD(3) * pkin(3) + t73 + (-t70 - t84) * t78;
t65 = t78 * qJD(2) + t80 * t70;
t63 = t80 * t84 + t65;
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t57 = t77 * t62 + t79 * t63;
t76 = cos(pkin(8));
t83 = -pkin(1) * t76 - pkin(2);
t56 = t79 * t62 - t63 * t77;
t68 = (-pkin(3) * t80 + t83) * qJD(1);
t74 = qJD(3) + qJD(4);
t71 = t83 * qJD(1);
t67 = (t77 * t80 + t78 * t79) * qJD(1);
t66 = (-t77 * t78 + t79 * t80) * qJD(1);
t64 = -t70 * t78 + t73;
t58 = -pkin(4) * t66 + qJD(5) + t68;
t55 = qJ(5) * t66 + t57;
t54 = pkin(4) * t74 - qJ(5) * t67 + t56;
t1 = qJD(2) ^ 2 * t85 + m(5) * (t56 ^ 2 + t57 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t58 ^ 2) / 0.2e1 + (t64 * mrSges(4,1) - t65 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t56 * mrSges(5,1) + t54 * mrSges(6,1) - t57 * mrSges(5,2) - t55 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t74) * t74 + (t68 * mrSges(5,2) + t58 * mrSges(6,2) - t56 * mrSges(5,3) - t54 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t67 + (Ifges(5,5) + Ifges(6,5)) * t74) * t67 + (-t68 * mrSges(5,1) - t58 * mrSges(6,1) + t57 * mrSges(5,3) + t55 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t66 + (Ifges(5,6) + Ifges(6,6)) * t74 + (Ifges(5,4) + Ifges(6,4)) * t67) * t66 + (t71 * (-mrSges(4,1) * t80 + mrSges(4,2) * t78) + (-t64 * t78 + t65 * t80) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t78 + Ifges(4,6) * t80) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (mrSges(3,1) * t76 - mrSges(3,2) * t75 + (t75 ^ 2 + t76 ^ 2) * t85 * pkin(1)) * pkin(1) + Ifges(4,2) * t80 ^ 2 / 0.2e1 + (Ifges(4,4) * t80 + Ifges(4,1) * t78 / 0.2e1) * t78) * qJD(1)) * qJD(1);
T = t1;
