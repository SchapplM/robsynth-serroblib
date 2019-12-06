% Calculate kinetic energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.21s
% Computational Cost: add. (191->66), mult. (379->97), div. (0->0), fcn. (188->4), ass. (0->24)
t65 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t78 = t65 * mrSges(4,3);
t77 = qJD(1) / 0.2e1;
t74 = cos(qJ(3));
t76 = -pkin(7) * qJD(1) + t65;
t59 = qJD(3) * pkin(3) + t76 * t74;
t72 = sin(qJ(3));
t60 = t76 * t72;
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t54 = t71 * t59 + t73 * t60;
t70 = qJD(1) * qJ(2);
t63 = t72 * qJD(1) * pkin(3) + t70;
t53 = t73 * t59 - t71 * t60;
t75 = qJD(1) ^ 2;
t69 = t75 * qJ(2) ^ 2;
t68 = qJD(3) + qJD(4);
t67 = -qJD(1) * pkin(1) + qJD(2);
t62 = (-t71 * t72 + t73 * t74) * qJD(1);
t61 = (-t71 * t74 - t72 * t73) * qJD(1);
t55 = -t61 * pkin(4) + qJD(5) + t63;
t52 = t61 * qJ(5) + t54;
t51 = t68 * pkin(4) - t62 * qJ(5) + t53;
t1 = m(4) * (t69 + (t72 ^ 2 + t74 ^ 2) * t65 ^ 2) / 0.2e1 + m(3) * (t67 ^ 2 + t69) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t63 ^ 2) / 0.2e1 + m(6) * (t51 ^ 2 + t52 ^ 2 + t55 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t75 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t74 * mrSges(4,1) - t72 * mrSges(4,2)) * t65) * qJD(3) + (t67 * mrSges(3,2) + (mrSges(4,2) * t70 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t77 - t78) * t74) * t74 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t74) * qJD(1) + (Ifges(4,2) * t77 - t78) * t72) * t72) * qJD(1) + (t53 * mrSges(5,1) + t51 * mrSges(6,1) - t54 * mrSges(5,2) - t52 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t68) * t68 + (t63 * mrSges(5,2) + t55 * mrSges(6,2) - t53 * mrSges(5,3) - t51 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t62 + (Ifges(5,5) + Ifges(6,5)) * t68) * t62 + (-t63 * mrSges(5,1) - t55 * mrSges(6,1) + t54 * mrSges(5,3) + t52 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t61 + (Ifges(5,6) + Ifges(6,6)) * t68 + (Ifges(5,4) + Ifges(6,4)) * t62) * t61;
T = t1;
