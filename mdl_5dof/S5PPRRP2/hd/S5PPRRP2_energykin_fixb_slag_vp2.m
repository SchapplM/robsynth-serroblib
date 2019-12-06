% Calculate kinetic energy for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:29
% EndTime: 2019-12-05 15:08:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (91->51), mult. (210->74), div. (0->0), fcn. (110->6), ass. (0->18)
t69 = sin(pkin(8));
t70 = cos(pkin(8));
t77 = qJD(1) * cos(qJ(3));
t78 = qJD(1) * sin(qJ(3));
t64 = t69 * t77 + t70 * t78;
t62 = qJD(3) * pkin(6) + t64;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t59 = t71 * qJD(2) + t73 * t62;
t63 = -t69 * t78 + t70 * t77;
t58 = t73 * qJD(2) - t71 * t62;
t76 = qJD(1) ^ 2;
t75 = qJD(2) ^ 2;
t61 = -qJD(3) * pkin(3) - t63;
t57 = (-pkin(4) * t73 - qJ(5) * t71 - pkin(3)) * qJD(3) - t63;
t56 = qJD(4) * qJ(5) + t59;
t55 = -qJD(4) * pkin(4) + qJD(5) - t58;
t1 = m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t75 + (t69 ^ 2 + t70 ^ 2) * t76) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t75) / 0.2e1 + m(2) * t76 / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t61 ^ 2) / 0.2e1 + (t58 * mrSges(5,1) - t55 * mrSges(6,1) - t59 * mrSges(5,2) + t56 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4)) * qJD(4) + (t63 * mrSges(4,1) - t64 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t61 * mrSges(5,1) - t57 * mrSges(6,1) + t56 * mrSges(6,2) + t59 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(3) * t73 + (Ifges(5,6) - Ifges(6,6)) * qJD(4)) * t73 + (t61 * mrSges(5,2) + t55 * mrSges(6,2) - t58 * mrSges(5,3) - t57 * mrSges(6,3) + (Ifges(6,4) + Ifges(5,5)) * qJD(4) + ((Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t71 + (Ifges(5,4) - Ifges(6,5)) * t73) * qJD(3)) * t71) * qJD(3);
T = t1;
