% Calculate kinetic energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:43
% EndTime: 2019-12-05 15:21:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (159->60), mult. (397->101), div. (0->0), fcn. (238->6), ass. (0->26)
t82 = sin(pkin(8));
t84 = cos(pkin(8));
t73 = qJD(3) + (-pkin(3) * t84 - qJ(4) * t82 - pkin(2)) * qJD(2);
t89 = qJ(3) * qJD(2);
t77 = qJD(1) * t82 + t84 * t89;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t67 = t81 * t73 + t83 * t77;
t91 = qJD(2) * t82;
t90 = qJD(2) * t84;
t88 = t81 * t91;
t66 = t83 * t73 - t77 * t81;
t76 = qJD(1) * t84 - t82 * t89;
t75 = qJD(4) - t76;
t86 = cos(qJ(5));
t85 = sin(qJ(5));
t80 = -qJD(2) * pkin(2) + qJD(3);
t78 = qJD(5) - t90;
t70 = (-t81 * t85 + t83 * t86) * t91;
t69 = (-t81 * t86 - t83 * t85) * t91;
t68 = pkin(4) * t88 + t75;
t65 = -pkin(6) * t88 + t67;
t64 = (-pkin(6) * t82 * t83 - pkin(4) * t84) * qJD(2) + t66;
t63 = t64 * t85 + t65 * t86;
t62 = t64 * t86 - t65 * t85;
t1 = m(4) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t75 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t78 / 0.2e1) * t78 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t78 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t78 + Ifges(6,2) * t69 / 0.2e1) * t69 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t75 * (mrSges(5,1) * t81 + mrSges(5,2) * t83) + t80 * mrSges(4,2) - t76 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) * t83 ^ 2 / 0.2e1 + (-Ifges(5,4) * t83 + Ifges(5,2) * t81 / 0.2e1) * t81) * t91 + (-t66 * t83 - t67 * t81) * mrSges(5,3)) * t82 + (t67 * mrSges(5,2) - t66 * mrSges(5,1) - t80 * mrSges(4,1) + t77 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t90 + (-Ifges(5,5) * t83 + Ifges(5,6) * t81 + Ifges(4,4)) * t91) * t84) * qJD(2);
T = t1;
