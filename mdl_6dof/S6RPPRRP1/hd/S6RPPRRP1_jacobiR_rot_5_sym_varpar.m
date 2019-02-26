% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:10
% EndTime: 2019-02-26 20:30:10
% DurationCPUTime: 0.03s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t76 = qJ(1) + pkin(9);
t72 = sin(t76);
t77 = sin(qJ(5));
t84 = t72 * t77;
t78 = cos(qJ(5));
t83 = t72 * t78;
t75 = pkin(10) + qJ(4);
t73 = cos(t75);
t82 = t73 * t77;
t81 = t73 * t78;
t74 = cos(t76);
t80 = t74 * t77;
t79 = t74 * t78;
t71 = sin(t75);
t70 = t73 * t79 + t84;
t69 = -t73 * t80 + t83;
t68 = -t72 * t81 + t80;
t67 = t72 * t82 + t79;
t1 = [t68, 0, 0, -t71 * t79, t69, 0; t70, 0, 0, -t71 * t83, -t67, 0; 0, 0, 0, t81, -t71 * t77, 0; t67, 0, 0, t71 * t80, -t70, 0; t69, 0, 0, t71 * t84, t68, 0; 0, 0, 0, -t82, -t71 * t78, 0; -t72 * t71, 0, 0, t74 * t73, 0, 0; t74 * t71, 0, 0, t72 * t73, 0, 0; 0, 0, 0, t71, 0, 0;];
JR_rot  = t1;
