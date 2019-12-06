% Return the minimum parameter vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t102 = (pkin(7) ^ 2);
t104 = (pkin(3) ^ 2);
t101 = (pkin(8) ^ 2);
t103 = (pkin(4) ^ 2);
t109 = (Ifges(5,2) + (t101 + t103) * m(6));
t113 = -pkin(8) * m(6) - mrSges(6,3);
t92 = (mrSges(5,3) - t113);
t106 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(7) * t92 + Ifges(6,2) + t109;
t99 = (m(5) + m(6));
t86 = Ifges(4,2) + (t102 + t104) * t99 + t106;
t87 = t102 * t99 + Ifges(4,1) + t106;
t114 = -t86 + t87;
t112 = -pkin(7) * t99 - t92;
t97 = sin(pkin(9));
t98 = cos(pkin(9));
t110 = t97 * t98;
t94 = t97 ^ 2;
t95 = t98 ^ 2;
t108 = t95 - t94;
t107 = Ifges(4,4) * t110;
t90 = pkin(3) * t99 + mrSges(4,1);
t105 = -t97 * mrSges(4,2) + t98 * t90;
t88 = t112 * pkin(3) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t94 * t87 + 0.2e1 * t107 + t95 * t86 + (2 * pkin(6) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t114 * t108 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t107; t108 * Ifges(4,4) + t114 * t110 + Ifges(3,4); -t97 * Ifges(4,6) + t98 * t88 + Ifges(3,5); t98 * Ifges(4,6) + t97 * t88 + Ifges(3,6); 0.2e1 * pkin(2) * t105 + (t104 * t99) + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t105; t98 * mrSges(4,2) + t97 * t90 + mrSges(3,2); mrSges(4,3) - t112; m(4) + t99; m(6) * t101 + Ifges(5,1) - t109; Ifges(5,4); t113 * pkin(4) + Ifges(5,5); Ifges(5,6); t103 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
