% Calculate kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:15
% EndTime: 2020-01-03 11:22:16
% DurationCPUTime: 0.44s
% Computational Cost: add. (285->79), mult. (764->131), div. (0->0), fcn. (500->8), ass. (0->36)
t90 = sin(pkin(7));
t92 = cos(pkin(8));
t105 = t90 * t92;
t88 = sin(pkin(9));
t91 = cos(pkin(9));
t93 = cos(pkin(7));
t77 = (t88 * t105 + t91 * t93) * qJD(1);
t106 = m(3) / 0.2e1;
t102 = qJD(1) * t93;
t79 = qJD(2) + (-pkin(2) * t93 - qJ(3) * t90 - pkin(1)) * qJD(1);
t89 = sin(pkin(8));
t101 = qJD(1) * qJ(2);
t99 = t93 * t101;
t72 = t89 * t79 + t92 * t99;
t68 = -qJ(4) * t102 + t72;
t103 = qJD(1) * t90;
t82 = t90 * t101 + qJD(3);
t74 = (pkin(3) * t89 - qJ(4) * t92) * t103 + t82;
t65 = t91 * t68 + t88 * t74;
t100 = t89 * t103;
t71 = t79 * t92 - t89 * t99;
t64 = -t68 * t88 + t74 * t91;
t67 = pkin(3) * t102 + qJD(4) - t71;
t95 = cos(qJ(5));
t94 = sin(qJ(5));
t85 = -qJD(1) * pkin(1) + qJD(2);
t78 = (t91 * t105 - t88 * t93) * qJD(1);
t75 = qJD(5) + t77;
t70 = t94 * t100 + t78 * t95;
t69 = t95 * t100 - t78 * t94;
t63 = pkin(6) * t100 + t65;
t62 = -pkin(4) * t100 - t64;
t61 = pkin(4) * t77 - pkin(6) * t78 + t67;
t60 = t61 * t94 + t63 * t95;
t59 = t61 * t95 - t63 * t94;
t1 = t85 ^ 2 * t106 + m(4) * (t71 ^ 2 + t72 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t62 ^ 2) / 0.2e1 + (t67 * mrSges(5,2) - t64 * mrSges(5,3) + Ifges(5,1) * t78 / 0.2e1) * t78 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,3) * t75 / 0.2e1) * t75 - (-t67 * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t78 - Ifges(5,2) * t77 / 0.2e1) * t77 + (t62 * mrSges(6,2) - t59 * mrSges(6,3) + Ifges(6,5) * t75 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t62 * mrSges(6,1) + t60 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t75 + Ifges(6,2) * t69 / 0.2e1) * t69 + ((-t85 * mrSges(3,1) - t71 * mrSges(4,1) + t72 * mrSges(4,2) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t102) * t93 + (t85 * mrSges(3,2) + (t82 * mrSges(4,2) - t71 * mrSges(4,3)) * t92 + (t82 * mrSges(4,1) + t64 * mrSges(5,1) - t65 * mrSges(5,2) - t72 * mrSges(4,3) + Ifges(5,5) * t78 - Ifges(5,6) * t77) * t89) * t90 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t106 + mrSges(3,3)) * (t90 ^ 2 + t93 ^ 2) * qJ(2) + ((Ifges(4,1) * t92 ^ 2 / 0.2e1 + Ifges(3,1) / 0.2e1) * t90 + (-Ifges(4,5) * t92 + Ifges(3,4)) * t93 + (Ifges(4,6) * t93 + (-Ifges(4,4) * t92 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t89) * t90) * t89) * t90) * qJD(1)) * qJD(1);
T = t1;
