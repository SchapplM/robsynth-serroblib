% Calculate kinetic energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:38
% EndTime: 2019-12-05 16:08:38
% DurationCPUTime: 0.21s
% Computational Cost: add. (163->66), mult. (377->99), div. (0->0), fcn. (210->6), ass. (0->24)
t72 = sin(qJ(2));
t67 = qJD(2) * pkin(6) + t72 * qJD(1);
t81 = t67 * mrSges(4,3);
t80 = qJD(2) / 0.2e1;
t71 = sin(qJ(3));
t76 = qJ(4) * qJD(2) + t67;
t61 = qJD(3) * pkin(3) - t76 * t71;
t73 = cos(qJ(3));
t62 = t76 * t73;
t70 = sin(pkin(8));
t79 = cos(pkin(8));
t58 = t70 * t61 + t79 * t62;
t78 = qJD(2) * t73;
t74 = cos(qJ(2));
t77 = t74 * qJD(1);
t57 = t79 * t61 - t70 * t62;
t65 = -t77 + qJD(4) + (-pkin(3) * t73 - pkin(2)) * qJD(2);
t68 = -qJD(2) * pkin(2) - t77;
t64 = (t70 * t73 + t79 * t71) * qJD(2);
t63 = t70 * t71 * qJD(2) - t79 * t78;
t56 = qJD(3) * qJ(5) + t58;
t55 = t63 * pkin(4) - t64 * qJ(5) + t65;
t54 = -qJD(3) * pkin(4) + qJD(5) - t57;
t1 = m(4) * (t68 ^ 2 + (t71 ^ 2 + t73 ^ 2) * t67 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + (m(3) * (t72 ^ 2 + t74 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t65 * mrSges(5,2) + t54 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t64) * t64 + (t65 * mrSges(5,1) + t55 * mrSges(6,1) - t56 * mrSges(6,2) - t58 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t63 + (-Ifges(5,4) + Ifges(6,5)) * t64) * t63 + (Ifges(3,3) * t80 + (t74 * mrSges(3,1) - t72 * mrSges(3,2)) * qJD(1) + (-t68 * mrSges(4,1) + (Ifges(4,2) * t80 + t81) * t73) * t73 + (Ifges(4,4) * t78 + t68 * mrSges(4,2) + (Ifges(4,1) * t80 + t81) * t71) * t71) * qJD(2) + (t57 * mrSges(5,1) - t54 * mrSges(6,1) - t58 * mrSges(5,2) + t56 * mrSges(6,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(3) + (-t71 * mrSges(4,1) - t73 * mrSges(4,2)) * t67 + (Ifges(6,4) + Ifges(5,5)) * t64 + (-Ifges(5,6) + Ifges(6,6)) * t63 + (Ifges(4,5) * t71 + Ifges(4,6) * t73) * qJD(2)) * qJD(3);
T = t1;
