% Calculate kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:01
% EndTime: 2019-12-05 18:20:01
% DurationCPUTime: 0.20s
% Computational Cost: add. (208->51), mult. (327->86), div. (0->0), fcn. (160->8), ass. (0->26)
t76 = qJD(1) + qJD(2);
t84 = cos(qJ(2));
t89 = qJD(1) * pkin(1);
t70 = t76 * pkin(2) + t84 * t89;
t78 = sin(pkin(8));
t80 = cos(pkin(8));
t82 = sin(qJ(2));
t88 = t82 * t89;
t68 = t78 * t70 + t80 * t88;
t66 = t76 * qJ(4) + t68;
t77 = sin(pkin(9));
t79 = cos(pkin(9));
t62 = -t79 * qJD(3) + t77 * t66;
t91 = t62 ^ 2;
t90 = t79 * t76;
t67 = t80 * t70 - t78 * t88;
t87 = qJD(4) - t67;
t83 = cos(qJ(5));
t81 = sin(qJ(5));
t71 = qJD(5) - t90;
t65 = -t76 * pkin(3) + t87;
t64 = t77 * qJD(3) + t79 * t66;
t61 = (-pkin(4) * t79 - pkin(7) * t77 - pkin(3)) * t76 + t87;
t60 = t81 * t61 + t83 * t64;
t59 = t83 * t61 - t81 * t64;
t1 = m(6) * (t59 ^ 2 + t60 ^ 2 + t91) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t91) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t82 ^ 2 + t84 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,3) * t71 / 0.2e1) * t71 + (t67 * mrSges(4,1) - t68 * mrSges(4,2) + (mrSges(3,1) * t84 - mrSges(3,2) * t82) * t89 + (-t65 * mrSges(5,1) + t64 * mrSges(5,3) + Ifges(5,2) * t90 / 0.2e1) * t79 + (t65 * mrSges(5,2) + (mrSges(6,1) * t81 + mrSges(6,2) * t83 + mrSges(5,3)) * t62 + (-t59 * t83 - t60 * t81) * mrSges(6,3) + t71 * (Ifges(6,5) * t83 - Ifges(6,6) * t81)) * t77 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t79 + (Ifges(6,1) * t83 ^ 2 / 0.2e1 + Ifges(5,1) / 0.2e1 + (-Ifges(6,4) * t83 + Ifges(6,2) * t81 / 0.2e1) * t81) * t77) * t77) * t76) * t76;
T = t1;
