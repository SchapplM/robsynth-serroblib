% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:41
% EndTime: 2019-02-26 21:36:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t73 = sin(qJ(2));
t74 = sin(qJ(1));
t81 = t74 * t73;
t72 = qJ(4) + pkin(9);
t70 = sin(t72);
t75 = cos(qJ(2));
t80 = t75 * t70;
t71 = cos(t72);
t79 = t75 * t71;
t76 = cos(qJ(1));
t78 = t76 * t73;
t77 = t76 * t75;
t69 = -t70 * t81 + t76 * t71;
t68 = t76 * t70 + t71 * t81;
t67 = t70 * t78 + t74 * t71;
t66 = -t74 * t70 + t71 * t78;
t1 = [t69, t70 * t77, 0, t66, 0, 0; t67, t74 * t80, 0, t68, 0, 0; 0, t73 * t70, 0, -t79, 0, 0; -t68, t71 * t77, 0, -t67, 0, 0; t66, t74 * t79, 0, t69, 0, 0; 0, t73 * t71, 0, t80, 0, 0; -t74 * t75, -t78, 0, 0, 0, 0; t77, -t81, 0, 0, 0, 0; 0, t75, 0, 0, 0, 0;];
JR_rot  = t1;
