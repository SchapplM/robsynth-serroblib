% Return the minimum parameter vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [17x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRRPP2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t88 = sin(pkin(8));
t86 = t88 ^ 2;
t89 = cos(pkin(8));
t87 = t89 ^ 2;
t99 = t87 - t86;
t98 = t88 * t89;
t91 = Ifges(5,2) + Ifges(6,3);
t94 = Ifges(5,1) + Ifges(6,1);
t97 = t91 - t94;
t93 = Ifges(5,4) - Ifges(6,5);
t96 = t93 * t98;
t95 = t89 * mrSges(5,1) - t88 * mrSges(5,2);
t92 = Ifges(5,5) + Ifges(6,4);
t90 = Ifges(5,6) - Ifges(6,6);
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + Ifges(4,2) + t86 * t94 + 0.2e1 * t96 + t87 * t91 + (2 * pkin(6) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(6) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -pkin(6) * m(4) + mrSges(3,2) - mrSges(4,3); -t99 * t97 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t96; t99 * t93 - t97 * t98 + Ifges(4,4); -t88 * t90 + t89 * t92 + Ifges(4,5); t88 * t92 + t89 * t90 + Ifges(4,6); 0.2e1 * pkin(3) * t95 + Ifges(6,2) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t95; t88 * mrSges(5,1) + t89 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
