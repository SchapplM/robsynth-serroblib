% Return the minimum parameter vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t107 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t97 = pkin(8) ^ 2;
t99 = (pkin(4) ^ 2);
t86 = Ifges(5,2) + (t97 + t99) * m(6) + t107;
t87 = m(6) * t97 + Ifges(5,1) + t107;
t109 = -t86 + t87;
t94 = sin(pkin(9));
t95 = cos(pkin(9));
t108 = t94 * t95;
t91 = t94 ^ 2;
t92 = t95 ^ 2;
t106 = t92 - t91;
t105 = Ifges(5,4) * t108;
t104 = -pkin(7) * m(4) - mrSges(4,3);
t103 = -pkin(8) * m(6) - mrSges(6,3);
t102 = (mrSges(3,3) - t104);
t89 = m(6) * pkin(4) + mrSges(5,1);
t101 = -t94 * mrSges(5,2) + t95 * t89;
t100 = pkin(2) ^ 2;
t98 = pkin(7) ^ 2;
t96 = (m(3) + m(4));
t90 = t98 + t100;
t88 = t103 * pkin(4) + Ifges(5,5);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t91 * t87 + 0.2e1 * t105 + t92 * t86 + (2 * pkin(7) * mrSges(4,3)) + t90 * m(4) + (2 * pkin(6) * t102) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * t96); pkin(1) * t96 + mrSges(2,1); -pkin(6) * t96 + mrSges(2,2) - t102; Ifges(3,1) - Ifges(3,2) + (-t90 + t98) * m(4); Ifges(3,4); t104 * pkin(2) + Ifges(3,5); Ifges(3,6); t100 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); t109 * t106 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t105; t106 * Ifges(5,4) + t109 * t108 + Ifges(4,4); -t94 * Ifges(5,6) + t95 * t88 + Ifges(4,5); t95 * Ifges(5,6) + t94 * t88 + Ifges(4,6); (t99 * m(6)) + 0.2e1 * pkin(3) * t101 + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t101; t95 * mrSges(5,2) + t94 * t89 + mrSges(4,2); mrSges(5,3) - t103; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
