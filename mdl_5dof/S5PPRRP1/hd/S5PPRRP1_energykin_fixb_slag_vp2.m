% Calculate kinetic energy for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:35
% EndTime: 2019-12-05 15:06:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->52), mult. (212->75), div. (0->0), fcn. (112->6), ass. (0->20)
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t79 = qJD(1) * cos(qJ(3));
t80 = qJD(1) * sin(qJ(3));
t64 = t70 * t79 + t71 * t80;
t62 = qJD(3) * pkin(6) + t64;
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t58 = t72 * qJD(2) + t74 * t62;
t78 = qJ(5) * qJD(3);
t63 = -t70 * t80 + t71 * t79;
t77 = qJD(1) ^ 2;
t76 = qJD(2) ^ 2;
t69 = t74 * qJD(2);
t61 = -qJD(3) * pkin(3) - t63;
t59 = qJD(5) + (-pkin(4) * t74 - pkin(3)) * qJD(3) - t63;
t57 = -t72 * t62 + t69;
t56 = t74 * t78 + t58;
t55 = qJD(4) * pkin(4) + t69 + (-t62 - t78) * t72;
t1 = m(5) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + m(3) * (t76 + (t70 ^ 2 + t71 ^ 2) * t77) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t76) / 0.2e1 + m(2) * t77 / 0.2e1 + (t57 * mrSges(5,1) + t55 * mrSges(6,1) - t58 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4)) * qJD(4) + (-t64 * mrSges(4,2) + t63 * mrSges(4,1) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t61 * mrSges(5,1) - t59 * mrSges(6,1) + t58 * mrSges(5,3) + t56 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(3) * t74 + (Ifges(5,6) + Ifges(6,6)) * qJD(4)) * t74 + (t61 * mrSges(5,2) + t59 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,5) + Ifges(6,5)) * qJD(4) + ((Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t72 + (Ifges(5,4) + Ifges(6,4)) * t74) * qJD(3)) * t72) * qJD(3);
T = t1;
