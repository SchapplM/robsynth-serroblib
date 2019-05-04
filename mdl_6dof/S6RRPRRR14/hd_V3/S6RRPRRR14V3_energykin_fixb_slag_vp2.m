% Calculate kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14V3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:18
% EndTime: 2019-04-12 15:03:19
% DurationCPUTime: 0.41s
% Computational Cost: add. (277->79), mult. (591->124), div. (0->0), fcn. (410->8), ass. (0->33)
t73 = sin(qJ(4));
t77 = cos(qJ(4));
t74 = sin(qJ(2));
t83 = qJD(1) * t74;
t67 = qJD(2) * t77 - t73 * t83;
t65 = t67 * qJ(3);
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t60 = -t76 * qJD(3) + t65 * t72;
t87 = t60 ^ 2;
t68 = qJD(2) * t73 + t77 * t83;
t63 = t68 * qJ(3);
t86 = t63 ^ 2;
t85 = m(4) / 0.2e1;
t84 = qJ(3) * mrSges(4,3);
t78 = cos(qJ(2));
t82 = qJD(1) * t78;
t69 = qJD(4) - t82;
t58 = -t68 * t72 + t69 * t76;
t81 = qJ(3) ^ 2;
t80 = qJD(1) ^ 2;
t79 = qJD(3) ^ 2;
t75 = cos(qJ(6));
t71 = sin(qJ(6));
t66 = qJD(5) - t67;
t62 = qJD(3) * t72 + t65 * t76;
t59 = t68 * t76 + t69 * t72;
t57 = qJD(6) - t58;
t56 = t62 * t75 + t63 * t71;
t55 = -t62 * t71 + t63 * t75;
t54 = t59 * t75 + t66 * t71;
t53 = -t59 * t71 + t66 * t75;
t1 = m(5) * (t65 ^ 2 + t79 + t86) / 0.2e1 + t80 * Ifges(2,3) / 0.2e1 + m(6) * (t62 ^ 2 + t86 + t87) / 0.2e1 + m(7) * (t55 ^ 2 + t56 ^ 2 + t87) / 0.2e1 + (t74 ^ 2 * t80 * t81 + t79) * t85 + (-t63 * mrSges(5,1) - t65 * mrSges(5,2) + Ifges(5,3) * t69 / 0.2e1) * t69 + (-t60 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t66 / 0.2e1) * t66 + (t55 * mrSges(7,1) - t56 * mrSges(7,2) + Ifges(7,3) * t57 / 0.2e1) * t57 + (qJD(3) * mrSges(5,2) + t63 * mrSges(5,3) + Ifges(5,5) * t69 + Ifges(5,1) * t68 / 0.2e1) * t68 + (t63 * mrSges(6,2) + t60 * mrSges(6,3) + Ifges(6,5) * t66 + Ifges(6,1) * t59 / 0.2e1) * t59 + (t60 * mrSges(7,2) - t55 * mrSges(7,3) + Ifges(7,5) * t57 + Ifges(7,1) * t54 / 0.2e1) * t54 + (-qJD(3) * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t68 + Ifges(5,6) * t69 + Ifges(5,2) * t67 / 0.2e1) * t67 + (-t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t59 + Ifges(6,6) * t66 + Ifges(6,2) * t58 / 0.2e1) * t58 + (-t60 * mrSges(7,1) + t56 * mrSges(7,3) + Ifges(7,4) * t54 + Ifges(7,6) * t57 + Ifges(7,2) * t53 / 0.2e1) * t53 + ((Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t78 ^ 2 * qJD(1) + (qJD(3) * mrSges(4,2) + (t84 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t83 + (mrSges(4,1) * qJ(3) + Ifges(3,4) - Ifges(4,5)) * t82) * t74) * qJD(1) + (-qJD(3) * mrSges(4,1) + (t81 * t85 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + t84) * qJD(2) + ((Ifges(4,4) + Ifges(3,5)) * t74 + (mrSges(4,2) * qJ(3) + Ifges(3,6) - Ifges(4,6)) * t78) * qJD(1)) * qJD(2);
T  = t1;
