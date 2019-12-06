% Return the minimum parameter vector for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t77 = (pkin(7) ^ 2);
t78 = (pkin(4) ^ 2);
t83 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t68 = Ifges(5,2) + (t77 + t78) * m(6) + t83;
t69 = m(6) * t77 + Ifges(5,1) + t83;
t86 = -t68 + t69;
t85 = (pkin(6) * m(4));
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t84 = t75 * t76;
t72 = t75 ^ 2;
t73 = t76 ^ 2;
t82 = t73 - t72;
t81 = Ifges(5,4) * t84;
t80 = -pkin(7) * m(6) - mrSges(6,3);
t71 = m(6) * pkin(4) + mrSges(5,1);
t79 = -t75 * mrSges(5,2) + t76 * t71;
t70 = t80 * pkin(4) + Ifges(5,5);
t1 = [0.2e1 * t81 + t73 * t68 + t72 * t69 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t85) * pkin(6)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t85; mrSges(3,3); m(3) + m(4); t86 * t82 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t81; t82 * Ifges(5,4) + t86 * t84 + Ifges(4,4); -t75 * Ifges(5,6) + t76 * t70 + Ifges(4,5); t76 * Ifges(5,6) + t75 * t70 + Ifges(4,6); (t78 * m(6)) + 0.2e1 * pkin(3) * t79 + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t79; t76 * mrSges(5,2) + t75 * t71 + mrSges(4,2); mrSges(5,3) - t80; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
