% Calculate kinetic energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:49:43
% EndTime: 2019-12-05 18:49:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (332->74), mult. (734->130), div. (0->0), fcn. (546->8), ass. (0->30)
t75 = qJD(2) + qJD(3);
t83 = cos(qJ(3));
t89 = pkin(2) * qJD(2);
t70 = pkin(3) * t75 + t83 * t89;
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t79 = sin(qJ(3));
t88 = t79 * t89;
t66 = t78 * t70 + t82 * t88;
t84 = cos(qJ(2));
t71 = (pkin(2) * t84 + pkin(1)) * qJD(1);
t80 = sin(qJ(2));
t67 = (t79 * t80 - t83 * t84) * qJD(1);
t64 = -pkin(3) * t67 + t71;
t68 = (-t79 * t84 - t80 * t83) * qJD(1);
t60 = t67 * t82 - t68 * t78;
t65 = t70 * t82 - t78 * t88;
t81 = cos(qJ(5));
t77 = sin(qJ(5));
t74 = qJD(4) + t75;
t63 = pkin(6) * t74 + t66;
t62 = -pkin(4) * t74 - t65;
t61 = t67 * t78 + t68 * t82;
t59 = qJD(5) - t60;
t58 = t61 * t81 + t74 * t77;
t57 = -t61 * t77 + t74 * t81;
t56 = -pkin(4) * t60 - pkin(6) * t61 + t64;
t55 = t56 * t77 + t63 * t81;
t54 = t56 * t81 - t63 * t77;
t1 = m(4) * (t71 ^ 2 + (t79 ^ 2 + t83 ^ 2) * pkin(2) ^ 2 * qJD(2) ^ 2) / 0.2e1 + Ifges(4,3) * t75 ^ 2 / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t62 ^ 2) / 0.2e1 + (t65 * mrSges(5,1) - t66 * mrSges(5,2) + Ifges(5,3) * t74 / 0.2e1) * t74 + (t54 * mrSges(6,1) - t55 * mrSges(6,2) + Ifges(6,3) * t59 / 0.2e1) * t59 + (t71 * mrSges(4,2) + Ifges(4,5) * t75 + Ifges(4,1) * t68 / 0.2e1) * t68 + (t64 * mrSges(5,2) - t65 * mrSges(5,3) + Ifges(5,5) * t74 + Ifges(5,1) * t61 / 0.2e1) * t61 + (t62 * mrSges(6,2) - t54 * mrSges(6,3) + Ifges(6,5) * t59 + Ifges(6,1) * t58 / 0.2e1) * t58 + (-t71 * mrSges(4,1) + Ifges(4,4) * t68 + Ifges(4,6) * t75 + Ifges(4,2) * t67 / 0.2e1) * t67 + (-t64 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,4) * t61 + Ifges(5,6) * t74 + Ifges(5,2) * t60 / 0.2e1) * t60 + (-t62 * mrSges(6,1) + t55 * mrSges(6,3) + Ifges(6,4) * t58 + Ifges(6,6) * t59 + Ifges(6,2) * t57 / 0.2e1) * t57 + (Ifges(3,3) * qJD(2) / 0.2e1 + (-Ifges(3,5) * t80 - Ifges(3,6) * t84) * qJD(1) + (t79 * (-mrSges(4,2) * t75 + mrSges(4,3) * t67) + t83 * (mrSges(4,1) * t75 - mrSges(4,3) * t68)) * pkin(2)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * pkin(1) ^ 2 / 0.2e1 + (pkin(1) * mrSges(3,1) + Ifges(3,2) * t84 / 0.2e1) * t84 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t84 + Ifges(3,1) * t80 / 0.2e1) * t80) * qJD(1) ^ 2;
T = t1;
