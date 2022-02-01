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
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:00
% EndTime: 2022-01-23 08:59:00
% DurationCPUTime: 0.43s
% Computational Cost: add. (285->79), mult. (764->131), div. (0->0), fcn. (500->8), ass. (0->36)
t93 = sin(pkin(7));
t95 = cos(pkin(8));
t108 = t95 * t93;
t91 = sin(pkin(9));
t94 = cos(pkin(9));
t96 = cos(pkin(7));
t80 = (t91 * t108 + t94 * t96) * qJD(1);
t109 = m(3) / 0.2e1;
t105 = t96 * qJD(1);
t104 = qJ(2) * qJD(1);
t102 = t96 * t104;
t82 = qJD(2) + (-pkin(2) * t96 - qJ(3) * t93 - pkin(1)) * qJD(1);
t92 = sin(pkin(8));
t75 = t95 * t102 + t92 * t82;
t71 = -qJ(4) * t105 + t75;
t106 = t93 * qJD(1);
t85 = t93 * t104 + qJD(3);
t77 = (pkin(3) * t92 - qJ(4) * t95) * t106 + t85;
t68 = t94 * t71 + t91 * t77;
t103 = t92 * t106;
t74 = -t92 * t102 + t82 * t95;
t67 = -t71 * t91 + t77 * t94;
t70 = pkin(3) * t105 + qJD(4) - t74;
t98 = cos(qJ(5));
t97 = sin(qJ(5));
t88 = -qJD(1) * pkin(1) + qJD(2);
t81 = (t94 * t108 - t91 * t96) * qJD(1);
t78 = qJD(5) + t80;
t73 = t97 * t103 + t81 * t98;
t72 = t98 * t103 - t81 * t97;
t66 = pkin(6) * t103 + t68;
t65 = -pkin(4) * t103 - t67;
t64 = pkin(4) * t80 - pkin(6) * t81 + t70;
t63 = t64 * t97 + t66 * t98;
t62 = t64 * t98 - t66 * t97;
t1 = t88 ^ 2 * t109 + m(4) * (t74 ^ 2 + t75 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t70 ^ 2) / 0.2e1 + (t70 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t81 / 0.2e1) * t81 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t78 / 0.2e1) * t78 - (-t70 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t81 - Ifges(5,2) * t80 / 0.2e1) * t80 + (t65 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t78 + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t65 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,6) * t78 + Ifges(6,2) * t72 / 0.2e1) * t72 + ((-t88 * mrSges(3,1) - t74 * mrSges(4,1) + t75 * mrSges(4,2) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t105) * t96 + (t88 * mrSges(3,2) + (t85 * mrSges(4,2) - t74 * mrSges(4,3)) * t95 + (t85 * mrSges(4,1) + t67 * mrSges(5,1) - t68 * mrSges(5,2) - t75 * mrSges(4,3) + Ifges(5,5) * t81 - Ifges(5,6) * t80) * t92) * t93 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t109 + mrSges(3,3)) * (t93 ^ 2 + t96 ^ 2) * qJ(2) + ((Ifges(3,1) / 0.2e1 + Ifges(4,1) * t95 ^ 2 / 0.2e1) * t93 + (-Ifges(4,5) * t95 + Ifges(3,4)) * t96 + (Ifges(4,6) * t96 + (-Ifges(4,4) * t95 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t92) * t93) * t92) * t93) * qJD(1)) * qJD(1);
T = t1;
