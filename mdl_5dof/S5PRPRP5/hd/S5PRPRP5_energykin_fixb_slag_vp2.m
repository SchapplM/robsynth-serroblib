% Calculate kinetic energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:30
% EndTime: 2019-12-05 15:37:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (150->59), mult. (358->91), div. (0->0), fcn. (206->6), ass. (0->23)
t74 = cos(pkin(8));
t83 = t74 ^ 2;
t82 = cos(qJ(4));
t73 = sin(pkin(8));
t76 = sin(qJ(2));
t69 = qJD(2) * qJ(3) + qJD(1) * t76;
t80 = pkin(6) * qJD(2) + t69;
t62 = t80 * t73;
t63 = t80 * t74;
t75 = sin(qJ(4));
t59 = -t75 * t62 + t82 * t63;
t81 = t73 ^ 2 + t83;
t77 = cos(qJ(2));
t79 = -t77 * qJD(1) + qJD(3);
t58 = -t62 * t82 - t75 * t63;
t66 = (-pkin(3) * t74 - pkin(2)) * qJD(2) + t79;
t67 = -qJD(2) * pkin(2) + t79;
t65 = (t73 * t82 + t74 * t75) * qJD(2);
t64 = (t73 * t75 - t74 * t82) * qJD(2);
t57 = qJD(4) * qJ(5) + t59;
t56 = -qJD(4) * pkin(4) + qJD(5) - t58;
t55 = t64 * pkin(4) - t65 * qJ(5) + t66;
t1 = m(4) * (t69 ^ 2 * t81 + t67 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + (m(3) * (t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t66 * mrSges(5,2) + t56 * mrSges(6,2) - t58 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t65) * t65 + (t66 * mrSges(5,1) + t55 * mrSges(6,1) - t57 * mrSges(6,2) - t59 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t64 + (-Ifges(5,4) + Ifges(6,5)) * t65) * t64 + (t58 * mrSges(5,1) - t56 * mrSges(6,1) - t59 * mrSges(5,2) + t57 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(4) + (Ifges(6,4) + Ifges(5,5)) * t65 + (-Ifges(5,6) + Ifges(6,6)) * t64) * qJD(4) + (t67 * (-mrSges(4,1) * t74 + mrSges(4,2) * t73) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t83 / 0.2e1 + (Ifges(4,4) * t74 + Ifges(4,1) * t73 / 0.2e1) * t73) * qJD(2) + (mrSges(3,1) * t77 - mrSges(3,2) * t76) * qJD(1) + t81 * t69 * mrSges(4,3)) * qJD(2);
T = t1;
