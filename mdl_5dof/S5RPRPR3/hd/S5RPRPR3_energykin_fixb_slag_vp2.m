% Calculate kinetic energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:11
% EndTime: 2019-12-05 17:51:11
% DurationCPUTime: 0.18s
% Computational Cost: add. (186->52), mult. (328->85), div. (0->0), fcn. (160->8), ass. (0->27)
t79 = cos(pkin(8));
t70 = (pkin(1) * t79 + pkin(2)) * qJD(1);
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t77 = sin(pkin(8));
t88 = pkin(1) * qJD(1) * t77;
t68 = t81 * t70 + t83 * t88;
t75 = qJD(1) + qJD(3);
t66 = t75 * qJ(4) + t68;
t76 = sin(pkin(9));
t78 = cos(pkin(9));
t62 = -t78 * qJD(2) + t76 * t66;
t91 = t62 ^ 2;
t90 = m(3) / 0.2e1;
t89 = t78 * t75;
t67 = t83 * t70 - t81 * t88;
t87 = qJD(4) - t67;
t84 = qJD(2) ^ 2;
t82 = cos(qJ(5));
t80 = sin(qJ(5));
t71 = qJD(5) - t89;
t65 = -t75 * pkin(3) + t87;
t64 = t76 * qJD(2) + t78 * t66;
t61 = (-pkin(4) * t78 - pkin(7) * t76 - pkin(3)) * t75 + t87;
t60 = t80 * t61 + t82 * t64;
t59 = t82 * t61 - t80 * t64;
t1 = m(4) * (t67 ^ 2 + t68 ^ 2 + t84) / 0.2e1 + t84 * t90 + m(5) * (t64 ^ 2 + t65 ^ 2 + t91) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t91) / 0.2e1 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t79 * mrSges(3,1) - t77 * mrSges(3,2) + (t77 ^ 2 + t79 ^ 2) * t90 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (-t68 * mrSges(4,2) + t67 * mrSges(4,1) + Ifges(4,3) * t75 / 0.2e1 + (-t65 * mrSges(5,1) + t64 * mrSges(5,3) + Ifges(5,2) * t89 / 0.2e1) * t78 + (t65 * mrSges(5,2) + (Ifges(5,4) * t78 + (Ifges(6,1) * t82 ^ 2 / 0.2e1 + Ifges(5,1) / 0.2e1 + (-Ifges(6,4) * t82 + Ifges(6,2) * t80 / 0.2e1) * t80) * t76) * t75 + (mrSges(6,1) * t80 + mrSges(6,2) * t82 + mrSges(5,3)) * t62 + (-t59 * t82 - t60 * t80) * mrSges(6,3) + t71 * (Ifges(6,5) * t82 - Ifges(6,6) * t80)) * t76) * t75;
T = t1;
