% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:14
% EndTime: 2019-02-26 20:50:14
% DurationCPUTime: 0.03s
% Computational Cost: add. (62->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t72 = pkin(11) + qJ(5);
t68 = sin(t72);
t74 = sin(qJ(3));
t79 = t74 * t68;
t70 = cos(t72);
t78 = t74 * t70;
t75 = cos(qJ(3));
t77 = t75 * t68;
t76 = t75 * t70;
t73 = qJ(1) + pkin(10);
t71 = cos(t73);
t69 = sin(t73);
t67 = t69 * t68 + t71 * t76;
t66 = t69 * t70 - t71 * t77;
t65 = t71 * t68 - t69 * t76;
t64 = t69 * t77 + t71 * t70;
t1 = [t65, 0, -t71 * t78, 0, t66, 0; t67, 0, -t69 * t78, 0, -t64, 0; 0, 0, t76, 0, -t79, 0; t64, 0, t71 * t79, 0, -t67, 0; t66, 0, t69 * t79, 0, t65, 0; 0, 0, -t77, 0, -t78, 0; -t69 * t74, 0, t71 * t75, 0, 0, 0; t71 * t74, 0, t69 * t75, 0, 0, 0; 0, 0, t74, 0, 0, 0;];
JR_rot  = t1;
