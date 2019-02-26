% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:45
% EndTime: 2019-02-26 21:23:45
% DurationCPUTime: 0.05s
% Computational Cost: add. (34->17), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
t85 = cos(qJ(2));
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t81 = sin(qJ(6));
t84 = cos(qJ(6));
t89 = t79 * t81 + t80 * t84;
t87 = t89 * t85;
t82 = sin(qJ(2));
t83 = sin(qJ(1));
t94 = t83 * t82;
t86 = cos(qJ(1));
t93 = t86 * t82;
t76 = t86 * t79 + t80 * t94;
t77 = -t79 * t94 + t86 * t80;
t92 = t76 * t84 - t77 * t81;
t91 = t76 * t81 + t77 * t84;
t90 = t79 * t84 - t80 * t81;
t88 = t90 * t85;
t75 = t79 * t93 + t83 * t80;
t74 = t83 * t79 - t80 * t93;
t73 = t74 * t81 + t75 * t84;
t72 = t74 * t84 - t75 * t81;
t1 = [t91, t86 * t88, 0, 0, 0, t72; t73, t83 * t88, 0, 0, 0, -t92; 0, t90 * t82, 0, 0, 0, t87; t92, -t86 * t87, 0, 0, 0, -t73; t72, -t83 * t87, 0, 0, 0, t91; 0, -t89 * t82, 0, 0, 0, t88; t83 * t85, t93, 0, 0, 0, 0; -t86 * t85, t94, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0;];
JR_rot  = t1;
