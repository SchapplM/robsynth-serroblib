% Calculate kinetic energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:29
% EndTime: 2019-12-05 18:47:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (296->67), mult. (417->103), div. (0->0), fcn. (218->6), ass. (0->26)
t72 = qJD(1) + qJD(2);
t85 = t72 / 0.2e1;
t75 = sin(qJ(2));
t83 = qJD(1) * pkin(1);
t69 = pkin(7) * t72 + t75 * t83;
t84 = t69 * mrSges(4,3);
t74 = sin(qJ(3));
t81 = pkin(8) * t72 + t69;
t63 = qJD(3) * pkin(3) - t81 * t74;
t77 = cos(qJ(3));
t64 = t81 * t77;
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t58 = t73 * t63 + t76 * t64;
t78 = cos(qJ(2));
t82 = t78 * t83;
t57 = t76 * t63 - t64 * t73;
t67 = -t82 + (-pkin(3) * t77 - pkin(2)) * t72;
t71 = qJD(3) + qJD(4);
t70 = -pkin(2) * t72 - t82;
t66 = (t73 * t77 + t74 * t76) * t72;
t65 = (-t73 * t74 + t76 * t77) * t72;
t59 = -pkin(4) * t65 + qJD(5) + t67;
t56 = qJ(5) * t65 + t58;
t55 = pkin(4) * t71 - qJ(5) * t66 + t57;
t1 = m(4) * (t70 ^ 2 + (t74 ^ 2 + t77 ^ 2) * t69 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t75 ^ 2 + t78 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t74 * mrSges(4,1) - t77 * mrSges(4,2)) * t69) * qJD(3) + (Ifges(3,3) * t85 + (mrSges(3,1) * t78 - mrSges(3,2) * t75) * t83 + (-t70 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t85 + t84) * t77) * t77 + (Ifges(4,4) * t72 * t77 + t70 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t85 + t84) * t74) * t74) * t72 + (t57 * mrSges(5,1) + t55 * mrSges(6,1) - t58 * mrSges(5,2) - t56 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t71) * t71 + (t67 * mrSges(5,2) + t59 * mrSges(6,2) - t57 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t66 + (Ifges(5,5) + Ifges(6,5)) * t71) * t66 + (-t67 * mrSges(5,1) - t59 * mrSges(6,1) + t58 * mrSges(5,3) + t56 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t65 + (Ifges(5,6) + Ifges(6,6)) * t71 + (Ifges(5,4) + Ifges(6,4)) * t66) * t65;
T = t1;
