% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:14
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR4_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:14:36
% EndTime: 2019-02-22 11:14:36
% DurationCPUTime: 0.04s
% Computational Cost: add. (26->8), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
t86 = cos(pkin(6));
t83 = sin(pkin(11));
t85 = cos(pkin(11));
t87 = sin(qJ(2));
t89 = cos(qJ(2));
t91 = t89 * t83 + t87 * t85;
t79 = t91 * t86;
t80 = t87 * t83 - t89 * t85;
t88 = sin(qJ(1));
t90 = cos(qJ(1));
t93 = t90 * t79 - t88 * t80;
t92 = t88 * t79 + t90 * t80;
t84 = sin(pkin(6));
t78 = t80 * t86;
t77 = -t88 * t78 + t90 * t91;
t76 = -t90 * t78 - t88 * t91;
t1 = [t90 * t84, 0, 0, 0, 0, 0; t88 * t84, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t93, t77, 0, 0, 0, 0; t92, -t76, 0, 0, 0, 0; 0, t80 * t84, 0, 0, 0, 0; t76, -t92, 0, 0, 0, 0; t77, t93, 0, 0, 0, 0; 0, t91 * t84, 0, 0, 0, 0;];
JR_rot  = t1;
