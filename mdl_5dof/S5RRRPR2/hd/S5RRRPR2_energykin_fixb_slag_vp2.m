% Calculate kinetic energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:31
% EndTime: 2022-01-20 11:30:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (207->43), mult. (302->76), div. (0->0), fcn. (140->8), ass. (0->23)
t76 = qJD(1) + qJD(2);
t84 = cos(qJ(2));
t88 = pkin(1) * qJD(1);
t73 = t76 * pkin(2) + t84 * t88;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t81 = sin(qJ(2));
t87 = t81 * t88;
t70 = t83 * t73 - t80 * t87;
t75 = qJD(3) + t76;
t68 = t75 * pkin(3) + t70;
t71 = t80 * t73 + t83 * t87;
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t66 = t77 * t68 + t78 * t71;
t65 = t78 * t68 - t77 * t71;
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t64 = t75 * pkin(8) + t66;
t63 = -t75 * pkin(4) - t65;
t62 = t79 * qJD(4) + t82 * t64;
t61 = t82 * qJD(4) - t79 * t64;
t1 = m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t81 ^ 2 + t84 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t76 / 0.2e1 + (mrSges(3,1) * t84 - mrSges(3,2) * t81) * t88) * t76 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t70 * mrSges(4,1) + t65 * mrSges(5,1) - t71 * mrSges(4,2) - t66 * mrSges(5,2) + (-t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t82 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t79 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) * t82 ^ 2 / 0.2e1 + (Ifges(6,4) * t82 + Ifges(6,1) * t79 / 0.2e1) * t79) * t75) * t75;
T = t1;
