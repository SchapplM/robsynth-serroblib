% Calculate kinetic energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:19
% EndTime: 2020-01-03 12:00:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (195->44), mult. (303->75), div. (0->0), fcn. (140->8), ass. (0->25)
t76 = qJD(1) + qJD(2);
t74 = qJD(4) + t76;
t90 = t74 / 0.2e1;
t84 = cos(qJ(2));
t89 = qJD(1) * pkin(1);
t73 = t76 * pkin(2) + t84 * t89;
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t81 = sin(qJ(2));
t88 = t81 * t89;
t70 = t78 * t73 - t77 * t88;
t68 = t76 * pkin(3) + t70;
t71 = t77 * t73 + t78 * t88;
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t66 = t80 * t68 + t83 * t71;
t65 = t83 * t68 - t80 * t71;
t85 = qJD(3) ^ 2;
t82 = cos(qJ(5));
t79 = sin(qJ(5));
t64 = t74 * pkin(8) + t66;
t63 = -t74 * pkin(4) - t65;
t62 = t79 * qJD(3) + t82 * t64;
t61 = t82 * qJD(3) - t79 * t64;
t1 = m(4) * (t70 ^ 2 + t71 ^ 2 + t85) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t85) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t81 ^ 2 + t84 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t66 * mrSges(5,2) + t65 * mrSges(5,1) + Ifges(5,3) * t90 + (Ifges(6,2) * t82 * t90 - t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t82 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t82 + Ifges(6,1) * t79 / 0.2e1) * t74) * t79) * t74 + (t70 * mrSges(4,1) - t71 * mrSges(4,2) + (mrSges(3,1) * t84 - mrSges(3,2) * t81) * t89 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t76) * t76;
T = t1;
