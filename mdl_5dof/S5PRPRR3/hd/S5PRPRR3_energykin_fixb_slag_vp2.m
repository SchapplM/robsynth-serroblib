% Calculate kinetic energy for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:33
% EndTime: 2019-12-05 15:46:33
% DurationCPUTime: 0.19s
% Computational Cost: add. (151->59), mult. (322->98), div. (0->0), fcn. (186->8), ass. (0->27)
t90 = cos(qJ(2));
t77 = qJD(2) * pkin(2) + qJD(1) * t90;
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t87 = sin(qJ(2));
t94 = qJD(1) * t87;
t73 = t83 * t77 + t84 * t94;
t71 = qJD(2) * pkin(6) + t73;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t67 = t86 * qJD(3) + t89 * t71;
t93 = t89 * qJD(2);
t72 = t77 * t84 - t83 * t94;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t82 = qJD(4) + qJD(5);
t81 = t89 * qJD(3);
t75 = (t85 * t89 + t86 * t88) * qJD(2);
t74 = (-t85 * t86 + t88 * t89) * qJD(2);
t70 = -qJD(2) * pkin(3) - t72;
t68 = (-pkin(4) * t89 - pkin(3)) * qJD(2) - t72;
t66 = -t71 * t86 + t81;
t65 = pkin(7) * t93 + t67;
t64 = qJD(4) * pkin(4) + t81 + (-pkin(7) * qJD(2) - t71) * t86;
t63 = t64 * t85 + t65 * t88;
t62 = t64 * t88 - t65 * t85;
t1 = m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + (m(3) * (t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,6) * t82 + Ifges(6,2) * t74 / 0.2e1) * t74 + (t72 * mrSges(4,1) - t73 * mrSges(4,2) + (t90 * mrSges(3,1) - t87 * mrSges(3,2)) * qJD(1) + (-t70 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t93 / 0.2e1) * t89 + (t70 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t86 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t89 + Ifges(5,1) * t86 / 0.2e1) * t86) * qJD(2)) * qJD(2);
T = t1;
