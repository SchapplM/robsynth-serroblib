% Calculate kinetic energy for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:38
% EndTime: 2019-12-05 17:07:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (173->54), mult. (255->92), div. (0->0), fcn. (122->6), ass. (0->26)
t78 = qJD(2) + qJD(3);
t83 = cos(qJ(4));
t90 = t78 * t83;
t81 = sin(qJ(3));
t89 = qJD(2) * pkin(2);
t73 = t78 * pkin(7) + t81 * t89;
t80 = sin(qJ(4));
t68 = t80 * qJD(1) + t83 * t73;
t84 = cos(qJ(3));
t88 = t84 * t89;
t86 = qJD(1) ^ 2;
t85 = qJD(2) ^ 2;
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t77 = qJD(4) + qJD(5);
t76 = t83 * qJD(1);
t74 = -t78 * pkin(3) - t88;
t71 = -t88 + (-pkin(4) * t83 - pkin(3)) * t78;
t70 = (t79 * t83 + t80 * t82) * t78;
t69 = (-t79 * t80 + t82 * t83) * t78;
t67 = -t80 * t73 + t76;
t66 = pkin(8) * t90 + t68;
t65 = qJD(4) * pkin(4) + t76 + (-pkin(8) * t78 - t73) * t80;
t64 = t79 * t65 + t82 * t66;
t63 = t82 * t65 - t79 * t66;
t1 = m(5) * (t67 ^ 2 + t68 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t86 + (t81 ^ 2 + t84 ^ 2) * pkin(2) ^ 2 * t85) / 0.2e1 + t85 * Ifges(3,3) / 0.2e1 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t77 / 0.2e1) * t77 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t71 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t77 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t71 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t77 + Ifges(6,2) * t69 / 0.2e1) * t69 + (Ifges(4,3) * t78 / 0.2e1 + (mrSges(4,1) * t84 - mrSges(4,2) * t81) * t89 + (-t74 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t90 / 0.2e1) * t83 + (t74 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t83 + Ifges(5,1) * t80 / 0.2e1) * t78) * t80) * t78 + (m(2) + m(3)) * t86 / 0.2e1;
T = t1;
