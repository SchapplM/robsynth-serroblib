% Calculate kinetic energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:31
% DurationCPUTime: 0.31s
% Computational Cost: add. (238->63), mult. (509->104), div. (0->0), fcn. (304->6), ass. (0->29)
t84 = cos(pkin(8));
t93 = t84 ^ 2;
t92 = m(3) / 0.2e1;
t83 = sin(pkin(8));
t75 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t90 = -pkin(6) * qJD(1) + t75;
t69 = t90 * t83;
t70 = t90 * t84;
t86 = sin(qJ(4));
t88 = cos(qJ(4));
t62 = t88 * t69 + t86 * t70;
t91 = t83 ^ 2 + t93;
t77 = qJD(1) * qJ(2) + qJD(3);
t73 = t83 * qJD(1) * pkin(3) + t77;
t61 = -t86 * t69 + t88 * t70;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t81 = qJD(4) + qJD(5);
t78 = -qJD(1) * pkin(1) + qJD(2);
t72 = (-t83 * t86 + t84 * t88) * qJD(1);
t71 = (-t83 * t88 - t84 * t86) * qJD(1);
t65 = -t71 * pkin(4) + t73;
t64 = t85 * t71 + t87 * t72;
t63 = t87 * t71 - t85 * t72;
t60 = t71 * pkin(7) + t62;
t59 = qJD(4) * pkin(4) - t72 * pkin(7) + t61;
t58 = t85 * t59 + t87 * t60;
t57 = t87 * t59 - t85 * t60;
t1 = t78 ^ 2 * t92 + m(5) * (t61 ^ 2 + t62 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t91 * t75 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t65 ^ 2) / 0.2e1 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t73 * mrSges(5,2) - t61 * mrSges(5,3) + Ifges(5,1) * t72 / 0.2e1) * t72 + (-t73 * mrSges(5,1) + t62 * mrSges(5,3) + Ifges(5,4) * t72 + Ifges(5,2) * t71 / 0.2e1) * t71 + (t65 * mrSges(6,2) - t57 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t64 / 0.2e1) * t64 + (-t65 * mrSges(6,1) + t58 * mrSges(6,3) + Ifges(6,4) * t64 + Ifges(6,6) * t81 + Ifges(6,2) * t63 / 0.2e1) * t63 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,5) * t72 + Ifges(5,6) * t71 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t77 * (mrSges(4,1) * t83 + mrSges(4,2) * t84) + t78 * mrSges(3,2) - t91 * t75 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t92 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t93 / 0.2e1 + (-Ifges(4,4) * t84 + Ifges(4,2) * t83 / 0.2e1) * t83) * qJD(1)) * qJD(1);
T = t1;
