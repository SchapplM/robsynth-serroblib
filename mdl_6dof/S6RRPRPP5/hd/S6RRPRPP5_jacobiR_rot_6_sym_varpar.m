% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiR_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:22:04
% EndTime: 2019-02-22 11:22:05
% DurationCPUTime: 0.03s
% Computational Cost: add. (15->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t69 = sin(qJ(2));
t70 = sin(qJ(1));
t80 = t70 * t69;
t71 = cos(qJ(4));
t79 = t70 * t71;
t68 = sin(qJ(4));
t72 = cos(qJ(2));
t78 = t72 * t68;
t77 = t72 * t71;
t73 = cos(qJ(1));
t76 = t73 * t69;
t75 = t73 * t71;
t74 = t73 * t72;
t67 = -t68 * t80 + t75;
t66 = t73 * t68 + t69 * t79;
t65 = t68 * t76 + t79;
t64 = t70 * t68 - t69 * t75;
t1 = [t67, t68 * t74, 0, -t64, 0, 0; t65, t70 * t78, 0, t66, 0, 0; 0, t69 * t68, 0, -t77, 0, 0; t66, -t71 * t74, 0, t65, 0, 0; t64, -t70 * t77, 0, -t67, 0, 0; 0, -t69 * t71, 0, -t78, 0, 0; t70 * t72, t76, 0, 0, 0, 0; -t74, t80, 0, 0, 0, 0; 0, -t72, 0, 0, 0, 0;];
JR_rot  = t1;
