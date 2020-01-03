% Return the minimum parameter vector for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% MPV [19x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t71 = sin(pkin(7));
t69 = t71 ^ 2;
t72 = cos(pkin(7));
t70 = t72 ^ 2;
t83 = t70 - t69;
t82 = (pkin(6) * m(4));
t81 = t71 * t72;
t74 = Ifges(5,2) + Ifges(6,3);
t77 = Ifges(5,1) + Ifges(6,1);
t80 = t74 - t77;
t76 = Ifges(5,4) - Ifges(6,5);
t79 = t76 * t81;
t78 = t72 * mrSges(5,1) - t71 * mrSges(5,2);
t75 = Ifges(5,5) + Ifges(6,4);
t73 = Ifges(5,6) - Ifges(6,6);
t1 = [0.2e1 * t79 + t69 * t77 + t70 * t74 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t82) * pkin(6)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t82; mrSges(3,3); m(3) + m(4); -t83 * t80 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t79; t83 * t76 - t80 * t81 + Ifges(4,4); -t71 * t73 + t72 * t75 + Ifges(4,5); t71 * t75 + t72 * t73 + Ifges(4,6); 0.2e1 * pkin(3) * t78 + Ifges(6,2) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t78; t71 * mrSges(5,1) + t72 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
