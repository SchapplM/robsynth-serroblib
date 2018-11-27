% Return the minimum parameter vector for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [45x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MPV = S7RRRRRRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_convert_par2_MPV_fixb: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_convert_par2_MPV_fixb: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_convert_par2_MPV_fixb: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_convert_par2_MPV_fixb: Ifges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t141 = (pkin(4) * m(8));
t131 = m(6) + m(7) + m(8);
t140 = (pkin(2) * (m(4) + m(5) + t131));
t139 = (pkin(3) * t131);
t138 = Ifges(4,2) + (2 * mrSges(4,3) + t140) * pkin(2);
t137 = Ifges(6,2) + (2 * mrSges(6,3) + t139) * pkin(3);
t136 = Ifges(8,2) + (2 * mrSges(8,3) + t141) * pkin(4);
t1 = [Ifges(2,3) + Ifges(3,2); mrSges(2,1); mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + t138; Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + t138; mrSges(3,1); mrSges(3,2) + mrSges(4,3) + t140; Ifges(4,1) + Ifges(5,2) - Ifges(4,2); Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + Ifges(5,2); mrSges(4,1); mrSges(4,2) + mrSges(5,3); Ifges(5,1) - Ifges(5,2) + t137; Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + t137; mrSges(5,1); mrSges(5,2) - mrSges(6,3) - t139; Ifges(6,1) + Ifges(7,2) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2) - mrSges(7,3); Ifges(7,1) - Ifges(7,2) + t136; Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3) + t136; mrSges(7,1); mrSges(7,2) + mrSges(8,3) + t141; Ifges(8,1) - Ifges(8,2); Ifges(8,4); Ifges(8,5); Ifges(8,6); Ifges(8,3); mrSges(8,1); mrSges(8,2);];
MPV  = t1;
