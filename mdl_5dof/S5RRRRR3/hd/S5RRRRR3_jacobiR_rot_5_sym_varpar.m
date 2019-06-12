% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-31 10:31
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-31 10:31:13
% EndTime: 2019-05-31 10:31:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (99->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t75 = qJ(4) + qJ(5);
t71 = sin(t75);
t76 = qJ(2) + qJ(3);
t72 = sin(t76);
t87 = t72 * t71;
t73 = cos(t75);
t86 = t72 * t73;
t74 = cos(t76);
t85 = t74 * t71;
t77 = sin(qJ(1));
t84 = t77 * t72;
t83 = t77 * t73;
t69 = t77 * t74;
t78 = cos(qJ(1));
t82 = t78 * t72;
t81 = t78 * t73;
t70 = t78 * t74;
t80 = t72 * t83;
t79 = t72 * t81;
t68 = t74 * t73;
t67 = t71 * t82;
t66 = t71 * t84;
t65 = t73 * t70 + t77 * t71;
t64 = -t71 * t70 + t83;
t63 = -t73 * t69 + t78 * t71;
t62 = t71 * t69 + t81;
t1 = [t63, -t79, -t79, t64, t64; t65, -t80, -t80, -t62, -t62; 0, t68, t68, -t87, -t87; t62, t67, t67, -t65, -t65; t64, t66, t66, t63, t63; 0, -t85, -t85, -t86, -t86; -t84, t70, t70, 0, 0; t82, t69, t69, 0, 0; 0, t72, t72, 0, 0;];
JR_rot  = t1;
